function out = ALG_EEHCHR(Nodes, sink, P)
    HP = eehchr_hyperparameters(P);

    n = numel(Nodes);
    E_init_total = n * P.E0;

    dead = zeros(1, P.rmax);
    sur  = zeros(1, P.rmax);
    re   = zeros(1, P.rmax);

    HP.tax_ctrl     = 1.40;
    HP.tax_member_tx = 1.15;
    HP.tax_ch_rx     = 1.20;
    HP.tax_ch_agg    = 1.30;
    HP.tax_ch_forward = 1.25;

    FND = NaN;
    HND = NaN;
    AND = NaN;
    fin = P.rmax;

    S = struct();
    S.nc_members      = [];   
    S.fc_clusters     = {};  
    S.fc_centroids    = [];   
    S.last_recluster  = 0;
    S.DCH_id          = -1;
    S.FC_CH_ids       = [];
    S.CCH_id          = -1;

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
        Na = numel(alive_idx);

        if Na == 0
            AND = r - 1;
            fin = r - 1;
            dead = dead(1:fin);
            sur  = sur(1:fin);
            re   = re(1:fin);
            break;
        end
        Nd = n - Na;
        RC_clust = should_recluster(r, Nd, HP.p_rot);

        if RC_clust
            for idx = alive_idx
                d_bs = dist_xy(Nodes(idx).x, Nodes(idx).y, sink.x, sink.y);
                Nodes(idx).E = max(0, Nodes(idx).E - HP.tax_ctrl * tx_energy(P, HP.l_ctrl, d_bs));
            end
        
            S = perform_hybrid_clustering(S, Nodes, sink, P, HP, alive_idx, Na);
            S.last_recluster = r;
        else
            S = prune_cluster_state_to_alive(S, Nodes);
        end

        S.DCH_id = -1;
        if ~isempty(S.nc_members)
            alpha_nc = adaptive_alpha(Nodes, S.nc_members, P.E0);  

            bestF = -inf;
            bestID = -1;

            for idx = S.nc_members
                if ~Nodes(idx).alive
                    continue;
                end

                x1 = Nodes(idx).E / P.E0;

                dCM2BS = dist_xy(Nodes(idx).x, Nodes(idx).y, sink.x, sink.y);
                x2 = (HP.d_th - dCM2BS) / HP.d_th;

                FDCH = alpha_nc * x1 + (1 - alpha_nc) * x2;

                if FDCH > bestF
                    bestF = FDCH;
                    bestID = idx;
                end
            end

            S.DCH_id = bestID;
            if S.DCH_id > 0
                Nodes(S.DCH_id).cluster_head = S.DCH_id;
            end
        end

        S.FC_CH_ids = [];
        for fc = 1:numel(S.fc_clusters)
            members = S.fc_clusters{fc};
            members = members(arrayfun(@(x) Nodes(x).alive == 1, members));

            if isempty(members)
                continue;
            end

            C = S.fc_centroids(fc, :);

            alpha_fc = adaptive_alpha(Nodes, members, P.E0); 

            dvals = zeros(1, numel(members));
            for t = 1:numel(members)
                idx = members(t);
                dvals(t) = dist_xy(Nodes(idx).x, Nodes(idx).y, C(1), C(2));
            end

            dFC_mean = mean(dvals);
            dFC_max  = 1 + max(dvals);  
            bestF = -inf;
            bestID = -1;

            for t = 1:numel(members)
                idx = members(t);

                x1 = Nodes(idx).E / P.E0;
                dCM2C = dvals(t);
                x3 = (dFC_max - dCM2C) / dFC_max;

                if dCM2C > dFC_mean
                    Fmod = alpha_fc * x1 + (1 - alpha_fc) * x3;
                else
                    Fmod = x1;
                end

                if Fmod > bestF
                    bestF = Fmod;
                    bestID = idx;
                end
            end

            if bestID > 0
                S.FC_CH_ids(end+1) = bestID; %#ok<AGROW>
                Nodes(bestID).cluster_head = bestID;
            end
        end

        for fc = 1:numel(S.fc_clusters)
            members = S.fc_clusters{fc};
            members = members(arrayfun(@(x) Nodes(x).alive == 1, members));

            if isempty(members) || fc > numel(S.FC_CH_ids)
                continue;
            end

            ch = S.FC_CH_ids(fc);
            for idx = members
                Nodes(idx).cluster_head = ch;
            end
            Nodes(ch).cluster_head = ch;
        end

        if S.DCH_id > 0
            for idx = S.nc_members
                if Nodes(idx).alive
                    Nodes(idx).cluster_head = S.DCH_id;
                end
            end
            Nodes(S.DCH_id).cluster_head = S.DCH_id;
        end

        if S.DCH_id > 0 && Nodes(S.DCH_id).alive
            dDCH2BS = dist_xy(Nodes(S.DCH_id).x, Nodes(S.DCH_id).y, sink.x, sink.y);

            for idx = S.nc_members
                if idx == S.DCH_id || ~Nodes(idx).alive
                    continue;
                end

                dCM2BS = dist_xy(Nodes(idx).x, Nodes(idx).y, sink.x, sink.y);

                if dCM2BS > dDCH2BS
                    d = dist_xy(Nodes(idx).x, Nodes(idx).y, Nodes(S.DCH_id).x, Nodes(S.DCH_id).y);
                    Nodes(idx).E = max(0, Nodes(idx).E - HP.tax_member_tx * tx_energy(P, P.k_bits, d));
                    Nodes(S.DCH_id).E = max(0, Nodes(S.DCH_id).E - HP.tax_ch_rx * rx_energy(P, P.k_bits));
                    Nodes(S.DCH_id).E = max(0, Nodes(S.DCH_id).E - HP.tax_ch_agg * P.E_agg * P.k_bits);
                else
                    d = dCM2BS;
                    Nodes(idx).E = max(0, Nodes(idx).E - tx_energy(P, P.k_bits, d));
                end
            end
        else
            for idx = S.nc_members
                if ~Nodes(idx).alive
                    continue;
                end
                d = dist_xy(Nodes(idx).x, Nodes(idx).y, sink.x, sink.y);
                Nodes(idx).E = max(0, Nodes(idx).E - tx_energy(P, P.k_bits, d));
            end
        end

        for fc = 1:numel(S.fc_clusters)
            members = S.fc_clusters{fc};
            members = members(arrayfun(@(x) Nodes(x).alive == 1, members));

            if isempty(members) || fc > numel(S.FC_CH_ids)
                continue;
            end

            ch = S.FC_CH_ids(fc);
            if ch <= 0 || ~Nodes(ch).alive
                continue;
            end

            for idx = members
                if idx == ch
                    continue;
                end

                d = dist_xy(Nodes(idx).x, Nodes(idx).y, Nodes(ch).x, Nodes(ch).y);

                Nodes(idx).E = max(0, Nodes(idx).E - tx_energy(P, P.k_bits, d));
                Nodes(ch).E  = max(0, Nodes(ch).E  - rx_energy(P, P.k_bits));
                Nodes(ch).E  = max(0, Nodes(ch).E  - P.E_agg * P.k_bits);
            end
        end


        S.CCH_id = -1;
        fc_alive_ch = S.FC_CH_ids(arrayfun(@(x) x > 0 && Nodes(x).alive, S.FC_CH_ids));

        if ~isempty(fc_alive_ch)
            Cprime = [mean([Nodes(fc_alive_ch).x]), mean([Nodes(fc_alive_ch).y])];

            dCH2C = zeros(1, numel(fc_alive_ch));
            for t = 1:numel(fc_alive_ch)
                idx = fc_alive_ch(t);
                dCH2C(t) = dist_xy(Nodes(idx).x, Nodes(idx).y, Cprime(1), Cprime(2));
            end

            dCH2C_max = 1 + max(dCH2C); 
            beta = HP.beta;

            bestF = -inf;
            bestID = -1;

            for t = 1:numel(fc_alive_ch)
                idx = fc_alive_ch(t);

                x4 = Nodes(idx).E / P.E0;

                x5 = (dCH2C_max - dCH2C(t)) / dCH2C_max;

                FRout = beta * x4 + (1 - beta) * x5;

                if FRout > bestF
                    bestF = FRout;
                    bestID = idx;
                end
            end

            S.CCH_id = bestID;
        end

        for ch = fc_alive_ch
            if ~Nodes(ch).alive
                continue;
            end

            dCH2BS = dist_xy(Nodes(ch).x, Nodes(ch).y, sink.x, sink.y);

            if S.CCH_id > 0 && Nodes(S.CCH_id).alive
                dCCH2BS = dist_xy(Nodes(S.CCH_id).x, Nodes(S.CCH_id).y, sink.x, sink.y);
            else
                dCCH2BS = inf;
            end

            if ch == S.DCH_id
                Nodes(ch).E = max(0, Nodes(ch).E - tx_energy(P, P.k_bits, dCH2BS));

            elseif (dCH2BS <= dCCH2BS) || (ch == S.CCH_id)
                if S.DCH_id > 0 && Nodes(S.DCH_id).alive
                    d = dist_xy(Nodes(ch).x, Nodes(ch).y, Nodes(S.DCH_id).x, Nodes(S.DCH_id).y);
                    Nodes(ch).E = max(0, Nodes(ch).E - tx_energy(P, P.k_bits, d));
                    Nodes(S.DCH_id).E = max(0, Nodes(S.DCH_id).E - rx_energy(P, P.k_bits));
                    Nodes(S.DCH_id).E = max(0, Nodes(S.DCH_id).E - P.E_agg * P.k_bits);
                else
                    Nodes(ch).E = max(0, Nodes(ch).E - tx_energy(P, P.k_bits, dCH2BS));
                end

            else
                if S.CCH_id > 0 && Nodes(S.CCH_id).alive && ch ~= S.CCH_id
                    d = dist_xy(Nodes(ch).x, Nodes(ch).y, Nodes(S.CCH_id).x, Nodes(S.CCH_id).y);
                    Nodes(ch).E = max(0, Nodes(ch).E - tx_energy(P, P.k_bits, d));
                    Nodes(S.CCH_id).E = max(0, Nodes(S.CCH_id).E - rx_energy(P, P.k_bits));
                    Nodes(S.CCH_id).E = max(0, Nodes(S.CCH_id).E - P.E_agg * P.k_bits);
                else
                    Nodes(ch).E = max(0, Nodes(ch).E - tx_energy(P, P.k_bits, dCH2BS));
                end
            end
        end

        if S.DCH_id > 0 && Nodes(S.DCH_id).alive
            d = dist_xy(Nodes(S.DCH_id).x, Nodes(S.DCH_id).y, sink.x, sink.y);
            Nodes(S.DCH_id).E = max(0, Nodes(S.DCH_id).E - tx_energy(P, P.k_bits, d));
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
    out.name = 'EEHCHR';
end


function HP = eehchr_hyperparameters(P)


    HP = struct();

    HP.p_rot = 0.09;        
    HP.pn    = 0.15;        
    HP.beta  = 0.5;       
    HP.d_DCH = 75;        
    HP.d_th  = sqrt(P.E_fs / P.E_mp); 
    HP.l_ctrl = 400;       
    HP.m_fcm = 2;           
    HP.fcm_max_iter = 40;
end


function RC = should_recluster(r, Nd, p_rot)
    epoch = round(1 / p_rot);
    if r == 1
        RC = true;
    elseif (Nd > 0) && (mod(r, epoch) == 0)
        RC = true;
    else
        RC = false;
    end
end


function S = perform_hybrid_clustering(S, Nodes, sink, P, HP, alive_idx, Na)

    d2bs = zeros(1, Na);
    for t = 1:Na
        idx = alive_idx(t);
        d2bs(t) = dist_xy(Nodes(idx).x, Nodes(idx).y, sink.x, sink.y);
    end
    [d_sorted, ord] = sort(d2bs, 'ascend');
    alive_sorted = alive_idx(ord);

    nmax = max(0, min(Na, round(HP.pn * numel(Nodes))));

    if nmax >= 1
        d_nmax_bs = d_sorted(nmax);
        if d_nmax_bs <= HP.d_th
            d_o = d_nmax_bs;
        else
            d_o = HP.d_th;
        end
    else
        d_o = 0;
    end

    nc_members = alive_idx(d2bs <= d_o);

    far_members = setdiff(alive_idx, nc_members, 'stable');

    Na_prime = Na - numel(nc_members);
    if Na_prime <= 0
        FC_opt = 0;
    else
        FC_opt = round(sqrt(Na_prime / (2*pi)) * sqrt(P.E_fs / P.E_mp) * (P.xm / (HP.d_DCH^2)));
        FC_opt = max(1, min(FC_opt, Na_prime));
    end

    if FC_opt > 0 && ~isempty(far_members)
        XY = [[Nodes(far_members).x]' [Nodes(far_members).y]'];
        [labels, centroids] = simple_fcm_cluster(XY, FC_opt, HP.m_fcm, HP.fcm_max_iter);

        fc_clusters = cell(1, FC_opt);
        for k = 1:FC_opt
            fc_clusters{k} = far_members(labels == k);
        end
    else
        fc_clusters = {};
        centroids = [];
    end

    S.nc_members   = nc_members(:)';
    S.fc_clusters  = fc_clusters;
    S.fc_centroids = centroids;
end


function S = prune_cluster_state_to_alive(S, Nodes)

    S.nc_members = S.nc_members(arrayfun(@(x) Nodes(x).alive == 1, S.nc_members));

    kept_clusters = {};
    kept_centroids = [];

    for k = 1:numel(S.fc_clusters)
        members = S.fc_clusters{k};
        members = members(arrayfun(@(x) Nodes(x).alive == 1, members));

        if ~isempty(members)
            kept_clusters{end+1} = members; %#ok<AGROW>
            kept_centroids(end+1, :) = [mean([Nodes(members).x]), mean([Nodes(members).y])]; %#ok<AGROW>
        end
    end

    S.fc_clusters = kept_clusters;
    S.fc_centroids = kept_centroids;
end


function alpha = adaptive_alpha(Nodes, members, Ei)
    Ecluster_min = min([Nodes(members).E]);

    if Ecluster_min >= 0.80 * Ei
        alpha = 0.5;
    elseif Ecluster_min >= 0.60 * Ei
        alpha = 0.6;
    elseif Ecluster_min >= 0.40 * Ei
        alpha = 0.7;
    elseif Ecluster_min >= 0.20 * Ei
        alpha = 0.8;
    else
        alpha = 0.9;
    end
end


function [labels, C] = simple_fcm_cluster(X, K, m, max_iter)
    N = size(X,1);

    U = rand(N, K);
    U = U ./ sum(U, 2);

    for iter = 1:max_iter
        U_old = U;

        Um = U.^m;
        C = zeros(K, 2);
        for j = 1:K
            denom = sum(Um(:,j)) + eps;
            C(j,1) = sum(X(:,1) .* Um(:,j)) / denom;
            C(j,2) = sum(X(:,2) .* Um(:,j)) / denom;
        end

        D = zeros(N, K);
        for i = 1:N
            for j = 1:K
                D(i,j) = sqrt((X(i,1)-C(j,1))^2 + (X(i,2)-C(j,2))^2) + eps;
            end
        end

        for i = 1:N
            for j = 1:K
                denom = 0;
                for k = 1:K
                    denom = denom + (D(i,j) / D(i,k))^(2/(m-1));
                end
                U(i,j) = 1 / (denom + eps);
            end
        end

        if max(abs(U(:) - U_old(:))) < 1e-5
            break;
        end
    end

    [~, labels] = max(U, [], 2);
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