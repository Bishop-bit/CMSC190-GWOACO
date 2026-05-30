function out = ALG_DCGWO(Nodes, sink, P)

    HP = dcgwo_hyperparameters(P);

    n = numel(Nodes);
    E_init_total = n * P.E0;

    dead = zeros(1, P.rmax);
    sur  = zeros(1, P.rmax);
    re   = zeros(1, P.rmax);

    FND = NaN;
    HND = NaN;
    AND = NaN;
    fin = P.rmax;

    for r = 1:P.rmax


        for i = 1:n
            if Nodes(i).E <= 0
                Nodes(i).E     = 0;
                Nodes(i).alive = 0;
            else
                Nodes(i).alive = 1;
            end
            Nodes(i).cluster_head = -1;
            Nodes(i).relay_load   = 0;
        end

        alive_idx = find([Nodes.alive] == 1);
        na = numel(alive_idx);

        if na == 0
            AND = r - 1;
            fin = r - 1;
            dead = dead(1:fin);
            sur  = sur(1:fin);
            re   = re(1:fin);
            break;
        end

        F = -inf(1, n);
        for k = 1:na
            i = alive_idx(k);
            F(i) = HP.w1 * Nodes(i).E + HP.w2 * rho(i);
        end

        ch_idx = select_potential_CHs(alive_idx, A, F, HP);

        if isempty(ch_idx)
            [~, pos] = max(F(alive_idx));
            ch_idx = alive_idx(pos);
        end

        cluster_of = assign_nodes_to_CHs_dijkstra(Nodes, alive_idx, ch_idx, A);

        for i = alive_idx
            cid = cluster_of(i);
            if cid > 0
                Nodes(i).cluster_head = ch_idx(cid);
            end
        end
        for c = 1:numel(ch_idx)
            Nodes(ch_idx(c)).cluster_head = ch_idx(c);
        end

        for i = alive_idx
            ch = Nodes(i).cluster_head;

            if ch <= 0 || ~Nodes(ch).alive || i == ch
                continue;
            end

            d = dist_xy(Nodes(i).x, Nodes(i).y, Nodes(ch).x, Nodes(ch).y);

            Etx = tx_energy(P, P.k_bits, d);
            Nodes(i).E = max(0, Nodes(i).E - Etx);

            Nodes(ch).E = max(0, Nodes(ch).E - rx_energy(P, P.k_bits));
            Nodes(ch).E = max(0, Nodes(ch).E - P.E_agg * P.k_bits);
        end

        alive_ch = ch_idx(arrayfun(@(x) Nodes(x).alive == 1, ch_idx));

        for cc = 1:numel(alive_ch)
            src = alive_ch(cc);

            if ~Nodes(src).alive
                continue;
            end

            best_path = gwo_route_path(src, alive_ch, Nodes, sink, HP);

            if isempty(best_path)
                d = dist_xy(Nodes(src).x, Nodes(src).y, sink.x, sink.y);
                Nodes(src).E = max(0, Nodes(src).E - tx_energy(P, P.k_bits, d));
                continue;
            end

            current = src;
            for h = 1:numel(best_path)
                nxt = best_path(h);

                if nxt == 0
                    d = dist_xy(Nodes(current).x, Nodes(current).y, sink.x, sink.y);
                    Nodes(current).E = max(0, Nodes(current).E - tx_energy(P, P.k_bits, d));
                    break;
                else
                    if ~Nodes(nxt).alive
                        break;
                    end

                    d = dist_xy(Nodes(current).x, Nodes(current).y, Nodes(nxt).x, Nodes(nxt).y);

                    Nodes(current).E = max(0, Nodes(current).E - tx_energy(P, P.k_bits, d));

                    Nodes(nxt).E = max(0, Nodes(nxt).E - rx_energy(P, P.k_bits));
                    Nodes(nxt).E = max(0, Nodes(nxt).E - P.E_agg * P.k_bits);
                    Nodes(nxt).relay_load = Nodes(nxt).relay_load + 1;

                    current = nxt;
                end
            end
        end

        for i = 1:n
            if Nodes(i).E <= 0
                Nodes(i).E     = 0;
                Nodes(i).alive = 0;
            else
                Nodes(i).alive = 1;
            end
        end

        alive_after = sum([Nodes.alive]);
        dead_after  = n - alive_after;
        total_E     = sum([Nodes.E]);

        dead(r) = dead_after;
        sur(r)  = 100 * alive_after / n;
        re(r)   = 100 * total_E / E_init_total;

        if isnan(FND) && dead_after >= 1
            FND = r;
        end
        if isnan(HND) && dead_after >= ceil(n/2)
            HND = r;
        end
        if dead_after == n
            AND = r;
            fin = r;
            dead = dead(1:fin);
            sur  = sur(1:fin);
            re   = re(1:fin);
            break;
        end
    end

    if isnan(AND)
        fin = numel(dead);
        dead = dead(1:fin);
        sur  = sur(1:fin);
        re   = re(1:fin);
    end

    out = struct();
    out.dead = dead;
    out.sur  = sur;
    out.re   = re;
    out.fin  = fin;
    out.FND  = FND;
    out.HND  = HND;
    out.AND  = AND;
    out.name = 'DC-GWO';
end


function HP = dcgwo_hyperparameters(P)
    HP = struct();
    HP.w1 = 0.6;
    HP.w2 = 0.4;
    HP.cluster_radius = 20;   
    HP.min_seed_gap   = 15;   
    HP.max_CH_frac    = 0.08; 
    HP.min_CH_frac    = 0.03; 
    HP.num_wolves = 10;
    HP.max_iter   = 15;
    HP.max_hops   = 3;
    HP.route_neighbor_radius = 22;
    HP.low_energy_threshold = 0.18;  
    HP.beta_sink = 0.5;
    HP.Pt_min = 0;
    HP.Pt_max = 1;
end

function [A, rho] = build_graph_and_density(Nodes, alive_idx, R)
    n = numel(Nodes);
    A = inf(n, n);
    rho = zeros(1, n);

    for ii = 1:numel(alive_idx)
        i = alive_idx(ii);
        A(i,i) = 0;
        for jj = ii+1:numel(alive_idx)
            j = alive_idx(jj);
            d = dist_xy(Nodes(i).x, Nodes(i).y, Nodes(j).x, Nodes(j).y);

            if d <= R
                A(i,j) = d;
                A(j,i) = d;
                rho(i) = rho(i) + 1;
                rho(j) = rho(j) + 1;
            end
        end
    end
end


function ch_idx = select_potential_CHs(alive_idx, A, F, HP)
    nAlive = numel(alive_idx);
    kmin = max(1, ceil(HP.min_CH_frac * nAlive));
    kmax = max(1, ceil(HP.max_CH_frac * nAlive));

    [~, ord] = sort(F(alive_idx), 'descend');
    ranked = alive_idx(ord);

    ch_idx = [];

    for t = 1:numel(ranked)
        i = ranked(t);

        too_close = false;
        for c = 1:numel(ch_idx)
            j = ch_idx(c);
            if isfinite(A(i,j)) && A(i,j) < HP.min_seed_gap
                too_close = true;
                break;
            end
        end

        if ~too_close
            ch_idx(end+1) = i; %#ok<AGROW>
        end

        if numel(ch_idx) >= kmax
            break;
        end
    end

    if numel(ch_idx) < kmin
        need = setdiff(ranked, ch_idx, 'stable');
        take = min(kmin - numel(ch_idx), numel(need));
        ch_idx = [ch_idx, need(1:take)];
    end
end

function cluster_of = assign_nodes_to_CHs_dijkstra(Nodes, alive_idx, ch_idx, A)
    n = numel(Nodes);
    cluster_of = zeros(1, n);

    for ii = 1:numel(alive_idx)
        i = alive_idx(ii);

        best_c = 0;
        best_d = inf;

        for c = 1:numel(ch_idx)
            s = ch_idx(c);
            d = dijkstra_shortest_path_cost(A, s, i);

            if isfinite(d) && d < best_d
                best_d = d;
                best_c = c;
            end
        end

        if best_c == 0
            ds = zeros(1, numel(ch_idx));
            for c = 1:numel(ch_idx)
                s = ch_idx(c);
                ds(c) = dist_xy(Nodes(i).x, Nodes(i).y, Nodes(s).x, Nodes(s).y);
            end
            [~, best_c] = min(ds);
        end

        cluster_of(i) = best_c;
    end
end

function dist_target = dijkstra_shortest_path_cost(A, src, dst)
    n = size(A,1);
    visited = false(1,n);
    distv = inf(1,n);
    distv(src) = 0;

    while true
        unvis = find(~visited);
        if isempty(unvis)
            break;
        end

        [~, pos] = min(distv(unvis));
        u = unvis(pos);

        if isinf(distv(u))
            break;
        end

        visited(u) = true;

        if u == dst
            break;
        end

        nbrs = find(isfinite(A(u,:)) & ~visited);
        for v = nbrs
            alt = distv(u) + A(u,v);
            if alt < distv(v)
                distv(v) = alt;
            end
        end
    end

    dist_target = distv(dst);
end

function best_path = gwo_route_path(src, alive_ch, Nodes, sink, HP)

    candidates = setdiff(alive_ch, src, 'stable');

    neigh = [];
    for i = 1:numel(candidates)
        j = candidates(i);
        d = dist_xy(Nodes(src).x, Nodes(src).y, Nodes(j).x, Nodes(j).y);
        if d <= HP.route_neighbor_radius
            neigh(end+1) = j; %#ok<AGROW>
        end
    end

    if isempty(neigh)
        best_path = 0;
        return;
    end

    wolves = cell(1, HP.num_wolves);
    fit = inf(1, HP.num_wolves);

    for w = 1:HP.num_wolves
        wolves{w} = random_path(src, neigh, Nodes, sink, HP);
        fit(w) = route_fitness(wolves{w}, src, Nodes, sink, HP);
    end

    [~, ord] = sort(fit, 'ascend');
    alpha = wolves{ord(1)};
    beta  = wolves{ord(min(2,end))};
    delta = wolves{ord(min(3,end))};

    for iter = 1:HP.max_iter
        a = 2 - 2 * (iter / HP.max_iter);

        for w = 1:HP.num_wolves
            curr = wolves{w};

            trial = discrete_gwo_update(curr, alpha, beta, delta, a, src, neigh, Nodes, sink, HP);

            ftrial = route_fitness(trial, src, Nodes, sink, HP);
            if ftrial < fit(w)
                wolves{w} = trial;
                fit(w) = ftrial;
            end
        end

        [~, ord] = sort(fit, 'ascend');
        alpha = wolves{ord(1)};
        beta  = wolves{ord(min(2,end))};
        delta = wolves{ord(min(3,end))};
    end

    best_path = alpha;
end

function path = random_path(src, neigh, Nodes, sink, HP)
    current = src;
    available = neigh(:)';
    path = [];

    for h = 1:HP.max_hops
        if isempty(available) || rand < 0.35
            path(end+1) = 0; %#ok<AGROW>
            return;
        end

        j = available(randi(numel(available)));
        path(end+1) = j; %#ok<AGROW>
        current = j; %#ok<NASGU>
        available = setdiff(available, j, 'stable');
    end

    if isempty(path) || path(end) ~= 0
        path(end+1) = 0;
    end
end

function trial = discrete_gwo_update(curr, alpha, beta, delta, a, src, neigh, Nodes, sink, HP) %#ok<INUSD>
    leaders = {alpha, beta, delta};

    pick = leaders{randi(3)};
    trial = pick;

    pmut = 0.2 + 0.4 * abs(a) / 2;
    if rand < pmut
        trial = random_path(src, neigh, Nodes, sink, HP);
    end

    if ~isempty(curr) && rand < 0.5
        L = min(numel(curr), numel(trial));
        if L >= 1
            cut = randi(L);
            trial(1:cut) = curr(1:cut);
        end
    end

    if isempty(trial) || trial(end) ~= 0
        if numel(trial) >= HP.max_hops
            trial(end) = 0;
        else
            trial(end+1) = 0;
        end
    end

    seen = [];
    cleaned = [];
    for k = 1:numel(trial)
        x = trial(k);
        if x == 0
            cleaned(end+1) = 0; %#ok<AGROW>
            break;
        end
        if ~ismember(x, seen)
            cleaned(end+1) = x; %#ok<AGROW>
            seen(end+1) = x; %#ok<AGROW>
        end
    end

    if isempty(cleaned) || cleaned(end) ~= 0
        cleaned(end+1) = 0;
    end

    trial = cleaned;
end

function f = route_fitness(path, src, Nodes, sink, HP)
    current = src;
    hop_costs = [];
    hop_energies = [];

    for k = 1:numel(path)
        nxt = path(k);

        if nxt == 0
            d = dist_xy(Nodes(current).x, Nodes(current).y, sink.x, sink.y);
            Ehop = d^2;
            hop_costs(end+1) = Ehop; %#ok<AGROW>
            break;
        else
            if ~Nodes(nxt).alive
                f = inf;
                return;
            end
            d = dist_xy(Nodes(current).x, Nodes(current).y, Nodes(nxt).x, Nodes(nxt).y);
            Ehop = d^2;
            hop_costs(end+1) = Ehop; %#ok<AGROW>
            hop_energies(end+1) = Nodes(nxt).E; %#ok<AGROW>
            current = nxt;
        end
    end

    if isempty(hop_costs)
        f = inf;
        return;
    end

    if isempty(hop_energies)
        Efac = Nodes(src).E;
        Vfac = 0;
    else
        Efac = mean(hop_energies);
        Vfac = var(hop_energies);
    end

    total_cost = sum(hop_costs);

    dsink = dist_xy(Nodes(src).x, Nodes(src).y, sink.x, sink.y);
    sink_balance = HP.beta_sink / (1 + dsink);

    f = total_cost / (Efac + eps) + Vfac - sink_balance;
end

function Etx = tx_energy(P, L, d)
    d0 = sqrt(P.E_fs / P.E_mp);

    if d <= d0
        Etx = L * P.E_elec + L * P.E_fs * d^2;
    else
        Etx = L * P.E_elec + L * P.E_mp * d^4;
    end
end


function Erx = rx_energy(P, L)
    Erx = L * P.E_elec;
end


function d = dist_xy(x1, y1, x2, y2)
    d = sqrt((x1 - x2)^2 + (y1 - y2)^2);
end