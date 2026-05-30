function out = ALG_LEACH(Nodes, sink, P)
    HP = leach_hyperparameters();
    n = numel(Nodes);
    E_init_total = n * P.E0;
    dead = zeros(1, P.rmax);
    sur  = zeros(1, P.rmax);
    re   = zeros(1, P.rmax);

    FND = NaN;
    HND = NaN; 
    AND = NaN;  

    epoch_len = round(1 / HP.p);

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
        num_alive = numel(alive_idx);
        if num_alive == 0
            AND = r - 1;
            fin = r - 1;
            dead = dead(1:fin);
            sur  = sur(1:fin);
            re   = re(1:fin);
            break;
        end

        ch_idx = [];

        for i = alive_idx

            if (r - Nodes(i).lastCH) >= epoch_len
                denom = 1 - HP.p * mod(r - 1, epoch_len);

                if denom <= 0
                    Tn = 1;
                else
                    Tn = HP.p / denom;
                end

                if rand <= Tn
                    ch_idx(end+1) = i; %#ok<AGROW>
                end
            end
        end

       
        if isempty(ch_idx)
            pick = alive_idx(randi(num_alive));
            ch_idx = pick;
        end

        for k = 1:numel(ch_idx)
            i = ch_idx(k);
            Nodes(i).lastCH = r;
            Nodes(i).cluster_head = i;
        end

      
        for i = alive_idx
            if any(ch_idx == i)
                continue; 
            end

            best_ch = -1;
            best_d  = inf;

            for c = 1:numel(ch_idx)
                j = ch_idx(c);

                d = euclid_dist(Nodes(i).x, Nodes(i).y, Nodes(j).x, Nodes(j).y);

                if d < best_d
                    best_d  = d;
                    best_ch = j;
                end
            end

            Nodes(i).cluster_head = best_ch;
        end

       
        members_per_ch = zeros(1, n);

        for i = alive_idx
            ch = Nodes(i).cluster_head;
            if ch > 0 && ch ~= i
                members_per_ch(ch) = members_per_ch(ch) + 1;
            end
        end

       
        for i = alive_idx
            ch = Nodes(i).cluster_head;

            if ch > 0 && ch ~= i && Nodes(ch).alive == 1

                d_to_ch = euclid_dist(Nodes(i).x, Nodes(i).y, ...
                                      Nodes(ch).x, Nodes(ch).y);

                Etx_member = tx_energy(P, P.k_bits, d_to_ch);

                Nodes(i).E = max(0, Nodes(i).E - Etx_member);

                Erx_ch = rx_energy(P, P.k_bits);
                Nodes(ch).E = max(0, Nodes(ch).E - Erx_ch);
            end
        end

        for c = 1:numel(ch_idx)
            ch = ch_idx(c);

            if Nodes(ch).E <= 0
                Nodes(ch).E     = 0;
                Nodes(ch).alive = 0;
                continue;
            end

            m = members_per_ch(ch);

            
            Eagg = m * P.k_bits * P.E_agg;

            d_to_bs = euclid_dist(Nodes(ch).x, Nodes(ch).y, sink.x, sink.y);
            Etx_bs  = tx_energy(P, P.k_bits, d_to_bs);

            Nodes(ch).E = max(0, Nodes(ch).E - Eagg - Etx_bs);
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
    out.name = 'LEACH';
end

function HP = leach_hyperparameters()
    HP = struct();
    HP.p = 0.05;
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

function d = euclid_dist(x1, y1, x2, y2)
    d = sqrt((x1 - x2)^2 + (y1 - y2)^2);
end