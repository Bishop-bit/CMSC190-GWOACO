function out = ALG_KEAC(Nodes, sink, P)
    HP = keac_hyperparameters();

    HP.tax_member_tx  = 1.10;
    HP.tax_ch_rx      = 1.15;
    HP.tax_ch_agg     = 1.20;
    HP.tax_ch_bs      = 1.25;

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
        mx = mean([Nodes(alive_idx).x]);
        my = mean([Nodes(alive_idx).y]);
        d_to_M = zeros(1, na);
        for t = 1:na
            idx = alive_idx(t);
            d_to_M(t) = dist_xy(Nodes(idx).x, Nodes(idx).y, mx, my);
        end
        D = mean(d_to_M);
        d_to_BS = zeros(1, na);
        for t = 1:na
            idx = alive_idx(t);
            d_to_BS(t) = dist_xy(Nodes(idx).x, Nodes(idx).y, sink.x, sink.y);
        end
        dBS = mean(d_to_BS);

        k_est = sqrt(na / (2*pi)) * sqrt(P.E_fs / P.E_mp) * (P.xm / (dBS^2 + eps));

        k = round(k_est);
        k = max(1, k);
        k = min(k, na);

        Cx = zeros(k,1);
        Cy = zeros(k,1);

        for ci = 1:k
            theta = 2*pi*(ci-1)/k;
            Cx(ci) = mx + D*cos(theta);
            Cy(ci) = my + D*sin(theta);
        end

        assign = zeros(na, 1);

        for iter = 1:HP.kmeans_max_iter
            prev_assign = assign;
            for t = 1:na
                idx = alive_idx(t);

                best_c = 1;
                best_d = inf;

                for ci = 1:k
                    dij = dist_xy(Nodes(idx).x, Nodes(idx).y, Cx(ci), Cy(ci));
                    if dij < best_d
                        best_d = dij;
                        best_c = ci;
                    end
                end

                assign(t) = best_c;
            end

            for ci = 1:k
                members_mask = (assign == ci);

                if ~any(members_mask)
                    continue;
                end

                members = alive_idx(members_mask);

                Cx(ci) = mean([Nodes(members).x]);
                Cy(ci) = mean([Nodes(members).y]);
            end
            if isequal(assign, prev_assign)
                break;
            end
        end

        ch_idx = [];

        for ci = 1:k
            members = alive_idx(assign == ci);

            if isempty(members)
                continue;
            end

            Wi = inf(1, numel(members));

            for m = 1:numel(members)
                idx = members(m);
                Di = dist_xy(Nodes(idx).x, Nodes(idx).y, Cx(ci), Cy(ci));
                Ei = Nodes(idx).E;
                Wi(m) = HP.c1 * (1 / (Ei + eps)) + HP.c2 * Di;
            end

            [~, best_pos] = min(Wi);
            ch = members(best_pos);

            ch_idx(end+1) = ch; %#ok<AGROW>

            for m = 1:numel(members)
                idx = members(m);
                Nodes(idx).cluster_head = ch;
            end
            Nodes(ch).cluster_head = ch;
        end

        if isempty(ch_idx)
            ch = alive_idx(randi(na));
            ch_idx = ch;
            for t = 1:na
                idx = alive_idx(t);
                Nodes(idx).cluster_head = ch;
            end
            Nodes(ch).cluster_head = ch;
        end

        for t = 1:na
            idx = alive_idx(t);

            if any(ch_idx == idx)
                continue;
            end

            ch = Nodes(idx).cluster_head;

            if ch <= 0 || ~Nodes(ch).alive
                continue;
            end

            d_cm_ch = dist_xy(Nodes(idx).x, Nodes(idx).y, Nodes(ch).x, Nodes(ch).y);

            Etx_member = HP.tax_member_tx * tx_energy(P, P.k_bits, d_cm_ch);
            Nodes(idx).E = max(0, Nodes(idx).E - Etx_member);

            Erx_ch     = HP.tax_ch_rx * rx_energy(P, P.k_bits);
            Nodes(ch).E = max(0, Nodes(ch).E - Erx_ch);

            if HP.use_CH_agg
                Eagg       = HP.tax_ch_agg * P.E_agg * P.k_bits;
                Nodes(ch).E = max(0, Nodes(ch).E - Eagg);
            end
        end

        for c = 1:numel(ch_idx)
            ch = ch_idx(c);

            if Nodes(ch).E <= 0
                Nodes(ch).E     = 0;
                Nodes(ch).alive = 0;
                continue;
            end

            d_ch_bs = dist_xy(Nodes(ch).x, Nodes(ch).y, sink.x, sink.y);
            Etx_bs     = HP.tax_ch_bs * tx_energy(P, P.k_bits, d_ch_bs);

            Nodes(ch).E = max(0, Nodes(ch).E - Etx_bs);
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
    out.name = 'KEAC';
end

function HP = keac_hyperparameters()

    HP = struct();

    HP.kmeans_max_iter = 12;
    HP.c1 = 0.6;
    HP.c2 = 1.8;

    HP.use_CH_agg = true;
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