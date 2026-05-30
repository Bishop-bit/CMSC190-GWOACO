function out = ALG_RLBEEP(Nodes, sink, P)

    HP = rlbeep_hyperparameters(P);

    n = numel(Nodes);
    E_init_total = n * P.E0;

    Q = zeros(n, n);

    node_mode = ones(1, n);

    restrict_counter = zeros(1, n);

    sleep_timer = zeros(1, n);

    data_cache_min =  inf(1, n);
    data_cache_max = -inf(1, n);
    first_active_run = true(1, n);

    [cluster_id, ch_idx, centroids] = bootstrap_clusters_and_CHs(Nodes, HP.K_clusters);

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
                node_mode(i)   = 0;
            else
                Nodes(i).alive = 1;
            end
            Nodes(i).cluster_head = -1;
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

        [cluster_id, ch_idx, centroids] = maintain_clusters_and_CHs( ...
            Nodes, cluster_id, ch_idx, centroids, HP);

        for i = alive_idx
            cid = cluster_id(i);
            if cid > 0 && cid <= numel(ch_idx) && ch_idx(cid) > 0 && Nodes(ch_idx(cid)).alive
                Nodes(i).cluster_head = ch_idx(cid);
            end
        end

        send_permission = false(1, n);

        for i = alive_idx

            if any(ch_idx == i)
                continue;
            end

            if node_mode(i) == 0
                continue;
            end

            sensed_value = synthetic_sensor_value(i, r, Nodes(i), sink);

            if first_active_run(i)
                data_cache_min(i) = inf;
                data_cache_max(i) = -inf;
                first_active_run(i) = false;
            end

            if sensed_value < data_cache_min(i)
                data_cache_min(i) = sensed_value;
            elseif sensed_value > data_cache_max(i)
                data_cache_max(i) = sensed_value;
            end

            cond1 = (data_cache_max(i) - sensed_value) > HP.change_threshold;
            cond2 = (sensed_value - data_cache_min(i)) < HP.change_threshold;

            if cond1 || cond2
                send_permission(i) = true;

             
                data_cache_min(i) = sensed_value;
                data_cache_max(i) = sensed_value;
            else
                send_permission(i) = false;
            end
        end

    
        for i = alive_idx

            if any(ch_idx == i)
                node_mode(i) = 1;
                restrict_counter(i) = 0;
                sleep_timer(i) = 0;
                continue;
            end

            if node_mode(i) == 1
                if ~send_permission(i)
                    restrict_counter(i) = restrict_counter(i) + 1;

                    if restrict_counter(i) >= HP.sleep_restrict_repeat_threshold
                        node_mode(i) = 0;
                        sleep_timer(i) = HP.sleep_interval_rounds;
                        restrict_counter(i) = 0;
                    end
                else
                    restrict_counter(i) = 0;
                end
            else
                sleep_timer(i) = max(0, sleep_timer(i) - 1);
                if sleep_timer(i) == 0
                    node_mode(i) = 1;
                    first_active_run(i) = true;
                end
            end
        end

        for i = alive_idx

            if any(ch_idx == i)
                continue;
            end

            if node_mode(i) == 0 || ~send_permission(i)
                continue;
            end

            ch = Nodes(i).cluster_head;
            if ch <= 0 || ~Nodes(ch).alive
                continue;
            end

            d_cm_ch = dist_xy(Nodes(i).x, Nodes(i).y, Nodes(ch).x, Nodes(ch).y);
            
            Etx = HP.tax_member_tx * tx_energy(P, P.k_bits, d_cm_ch);
            Nodes(i).E = max(0, Nodes(i).E - Etx);
            
            Erx = HP.tax_ch_rx * rx_energy(P, P.k_bits);
            Nodes(ch).E = max(0, Nodes(ch).E - Erx);
            
            Eagg = HP.tax_ch_agg * P.E_agg * P.k_bits;
            Nodes(ch).E = max(0, Nodes(ch).E - Eagg);
        end

        h_to_sink = estimate_hop_count_to_sink(Nodes, ch_idx, sink, HP.send_distance_range);

        for c = 1:numel(ch_idx)
            cur = ch_idx(c);

            if cur <= 0 || ~Nodes(cur).alive
                continue;
            end

            nbrs = candidate_forwarders(cur, ch_idx, Nodes, HP.send_distance_range);

            d_sink = dist_xy(Nodes(cur).x, Nodes(cur).y, sink.x, sink.y);
            can_send_direct = (d_sink <= HP.send_distance_range);

            if can_send_direct
                Etx = tx_energy(P, P.k_bits, d_sink);
                Nodes(cur).E = max(0, Nodes(cur).E - Etx);
                continue;
            end

            if isempty(nbrs)
                Nodes(cur).E = max(0, Nodes(cur).E - HP.tax_failed_route * tx_energy(P, P.k_bits, d_sink));
                continue;
            end

            best_q = -inf;
            best_nbr = -1;

            for kk = 1:numel(nbrs)
                nbr = nbrs(kk);

                d_cn = dist_xy(Nodes(cur).x, Nodes(cur).y, Nodes(nbr).x, Nodes(nbr).y);

                nd = d_cn / max(HP.MNDlon, HP.MNDlat);

                n_exp = nd * (HP.DFRmax - HP.DFRmin) + HP.DFRmin;

                h_nbr = max(1, h_to_sink(nbr));
                R = Nodes(nbr).E / ((d_cn^n_exp) * h_nbr + eps);

                Qnbr = max(Q(nbr, :));

                Q(cur, nbr) = (1 - HP.alpha) * Q(cur, nbr) + HP.alpha * (R + Qnbr);

                if Q(cur, nbr) > best_q
                    best_q = Q(cur, nbr);
                    best_nbr = nbr;
                end
            end

            if best_nbr > 0 && Nodes(best_nbr).alive
                d_cn = dist_xy(Nodes(cur).x, Nodes(cur).y, Nodes(best_nbr).x, Nodes(best_nbr).y);

                Etx = HP.tax_ch_forward * tx_energy(P, P.k_bits, d_cn);
                Nodes(cur).E = max(0, Nodes(cur).E - Etx);
                
                Erx = HP.tax_ch_rx * rx_energy(P, P.k_bits);
                Nodes(best_nbr).E = max(0, Nodes(best_nbr).E - Erx);
                
                Eagg = HP.tax_ch_agg * P.E_agg * P.k_bits;
                Nodes(best_nbr).E = max(0, Nodes(best_nbr).E - Eagg);
            end
        end

        for i = 1:n
            if Nodes(i).E <= 0
                Nodes(i).E     = 0;
                Nodes(i).alive = 0;
                node_mode(i)   = 0;
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
    out.name = 'RLBEEP';
end

function HP = rlbeep_hyperparameters(P)
    HP = struct();
    HP.K_clusters = 4;
    HP.send_distance_range = 18;  
    HP.alpha = 0.5;
    HP.DFRmin = 5.0;
    HP.DFRmax = 55.0;
    HP.MNDlon = 60.0;
    HP.MNDlat = 60.0;
    HP.total_epochs = 300;
    HP.tax_member_tx   = 1.3;
    HP.tax_ch_rx       = 1.5;
    HP.tax_ch_agg      = 1.5;
    HP.tax_ch_forward  = 1.50;
    HP.tax_failed_route = 1.50;
    HP.change_threshold = 0.25;
    HP.sleep_restrict_repeat_threshold = 10;
    HP.sleep_interval_rounds = 1;
    HP.reselect_CH_on_death_only = true;
    HP.packet_bits = P.k_bits;
end

function [cluster_id, ch_idx, centroids] = bootstrap_clusters_and_CHs(Nodes, K)
    n = numel(Nodes);
    XY = [[Nodes.x]' [Nodes.y]'];

    centroids = XY(randperm(n, K), :);
    cluster_id = ones(n, 1);

    for iter = 1:20
        for i = 1:n
            ds = sum((centroids - XY(i,:)).^2, 2);
            [~, cluster_id(i)] = min(ds);
        end

        for k = 1:K
            members = find(cluster_id == k);
            if ~isempty(members)
                centroids(k,:) = mean(XY(members,:), 1);
            end
        end
    end

    ch_idx = zeros(1, K);
    for k = 1:K
        members = find(cluster_id == k);
        if isempty(members)
            ch_idx(k) = 0;
            continue;
        end
        ds = zeros(1, numel(members));
        for t = 1:numel(members)
            i = members(t);
            ds(t) = dist_xy(Nodes(i).x, Nodes(i).y, centroids(k,1), centroids(k,2));
        end
        [~, pos] = min(ds);
        ch_idx(k) = members(pos);
    end
end

function [cluster_id, ch_idx, centroids] = maintain_clusters_and_CHs(Nodes, cluster_id, ch_idx, centroids, HP)
    K = numel(ch_idx);
    n = numel(Nodes);

    for k = 1:K
        members = find(cluster_id == k & [Nodes.alive]' == 1);

        if isempty(members)
            ch_idx(k) = 0;
            continue;
        end

        xk = mean([Nodes(members).x]);
        yk = mean([Nodes(members).y]);
        centroids(k,:) = [xk yk];

        if ch_idx(k) == 0 || ~Nodes(ch_idx(k)).alive
            ds = zeros(1, numel(members));
            for t = 1:numel(members)
                i = members(t);
                ds(t) = dist_xy(Nodes(i).x, Nodes(i).y, xk, yk);
            end
            [~, pos] = min(ds);
            ch_idx(k) = members(pos);
        end
    end

    alive_idx = find([Nodes.alive] == 1);
    valid_clusters = find(ch_idx > 0);

    if isempty(valid_clusters)
        return;
    end

    for i = alive_idx
        ds = inf(1, numel(valid_clusters));
        for t = 1:numel(valid_clusters)
            k = valid_clusters(t);
            ds(t) = dist_xy(Nodes(i).x, Nodes(i).y, centroids(k,1), centroids(k,2));
        end
        [~, pos] = min(ds);
        cluster_id(i) = valid_clusters(pos);
    end

    for i = 1:n
        if ~Nodes(i).alive
            cluster_id(i) = 0;
        end
    end
end

function nbrs = candidate_forwarders(cur, ch_idx, Nodes, send_range)
    nbrs = [];
    for k = 1:numel(ch_idx)
        nbr = ch_idx(k);
        if nbr <= 0 || nbr == cur || ~Nodes(nbr).alive
            continue;
        end
        d = dist_xy(Nodes(cur).x, Nodes(cur).y, Nodes(nbr).x, Nodes(nbr).y);
        if d <= send_range
            nbrs(end+1) = nbr; %#ok<AGROW>
        end
    end
end


function h = estimate_hop_count_to_sink(Nodes, ch_idx, sink, send_range)
    n = numel(Nodes);
    h = inf(1, n);

    alive_ch = ch_idx(ch_idx > 0);
    if isempty(alive_ch)
        return;
    end

    frontier = [];
    for ii = 1:numel(alive_ch)
        i = alive_ch(ii);
        d = dist_xy(Nodes(i).x, Nodes(i).y, sink.x, sink.y);
        if d <= send_range
            h(i) = 1;
            frontier(end+1) = i; %#ok<AGROW>
        end
    end

    while ~isempty(frontier)
        cur = frontier(1);
        frontier(1) = [];

        for ii = 1:numel(alive_ch)
            j = alive_ch(ii);
            if j == cur
                continue;
            end
            d = dist_xy(Nodes(j).x, Nodes(j).y, Nodes(cur).x, Nodes(cur).y);
            if d <= send_range && h(j) > h(cur) + 1
                h(j) = h(cur) + 1;
                frontier(end+1) = j; %#ok<AGROW>
            end
        end
    end

    h(isinf(h)) = 1e6;
end

function val = synthetic_sensor_value(i, r, node, sink)
    base = 20 + 2*sin(0.05*r + 0.2*i);
    spatial = 0.01 * dist_xy(node.x, node.y, sink.x, sink.y);
    noise = 0.15 * sin(0.17*r + i);
    val = base + spatial + noise;
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