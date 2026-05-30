function out = ALG_GWOACO(Nodes0, sink, P)
n = numel(Nodes0);

for i = 1:n
    Nodes0(i).alive        = (Nodes0(i).E > 0);
    Nodes0(i).cluster_head = -1;
    Nodes0(i).relay_load   = 0;
end

S = init_gwo_aco_state_clean(Nodes0, sink, P, n);

re   = zeros(1, P.rmax);
dead = zeros(1, P.rmax);
sur  = zeros(1, P.rmax);

fin = P.rmax;
FND = NaN;
HND = NaN;
AND = NaN;

for r = 1:P.rmax
    if S.active
        S = step_gwo_aco(S, sink, P, n);
        if S.num_alive == 0
            S.active = false;
            S.final_round = r;
            fin = r;
        end
    end

    alive = sum([S.Nodes.alive]);
    e     = sum([S.Nodes.E]);

    re(r)   = (e / (n * P.E0)) * 100;
    dead(r) = n - alive;
    sur(r)  = (alive / n) * 100;

    if isnan(FND) && dead(r) >= 1
        FND = r;
    end
    if isnan(HND) && dead(r) >= ceil(n/2)
        HND = r;
    end
    if dead(r) == n
        AND = r;
    end

    if ~S.active
        break;
    end
end

re(fin+1:end)   = [];
dead(fin+1:end) = [];
sur(fin+1:end)  = [];

if isnan(AND) && ~isempty(dead) && dead(end) == n
    AND = numel(dead);
end

out = struct();
out.dead = dead;
out.sur  = sur;
out.re   = re;
out.fin  = fin;
out.FND  = FND;
out.HND  = HND;
out.AND  = AND;
out.name = 'GWO-ACO';
end

function S = init_gwo_aco_state_clean(Nodes, sink, P, n)
S.Nodes = Nodes;
S.active = true;
S.final_round = P.rmax;
S.survival = [];
S.res_energy = 0;
S.num_wolves = 50;
S.max_iter   = 25;

S.ACO.num_ants = 15;
S.ACO.alpha    = 1;
S.ACO.beta     = 4;
S.ACO.rho      = 0.2;
S.ACO.Q        = 1;

S.pheromone = ones(n+1, n+1);
S.ch_indices = [];
S.num_alive  = n;
end

function S = step_gwo_aco(S, sink, P, n)
alive_idx = find([S.Nodes.alive] == 1);

for i = alive_idx
    S.Nodes(i).relay_load = 0;
end

S.num_alive = numel(alive_idx);
if S.num_alive == 0
    return;
end

if S.num_alive == 1
    i = alive_idx(1);
    d = sqrt((S.Nodes(i).x - sink.x)^2 + (S.Nodes(i).y - sink.y)^2);
    S.Nodes(i).E = max(0, S.Nodes(i).E - tx_energy(P, d));
    if S.Nodes(i).E <= 0
        S.Nodes(i).alive = 0;
    end
    S.num_alive = sum([S.Nodes.alive]);
    S.ch_indices = [];
    return;
end

na = S.num_alive;

d_sum = 0;
for ii = alive_idx
    d_sum = d_sum + sqrt((S.Nodes(ii).x - sink.x)^2 + (S.Nodes(ii).y - sink.y)^2);
end
dBS_alive = d_sum / na;

k_opt = sqrt(na/(2*pi)) * sqrt(P.E_fs/P.E_mp) * (P.xm/(dBS_alive^2));

k_min = ceil(0.04 * na);
k_max = ceil(0.12 * na);

k = round(k_opt);
k = max(k_min, min(k, k_max));
k = max(1, min(k, na));

wolves = zeros(S.num_wolves, k);
for w = 1:S.num_wolves
    wolves(w,:) = randperm(na, k);
end

fitness = zeros(1, S.num_wolves);
for w = 1:S.num_wolves
    real_ids = alive_idx(wolves(w,:));
    fitness(w) = gwo_fitness_aco(real_ids, S.Nodes, sink, P, n);
end

[~, si] = sort(fitness);
ap = wolves(si(1), :);
bp = wolves(si(min(2,end)), :);
dp = wolves(si(min(3,end)), :);

for iter = 1:S.max_iter
    a = 2 - 2*(iter / S.max_iter);

    for w = 1:S.num_wolves
        nw2 = zeros(1, k);

        for g = 1:k
            r1 = rand; r2 = rand;
            X1 = ap(g) - (2*a*r1 - a) * abs(2*r2*ap(g) - wolves(w,g));

            r1 = rand; r2 = rand;
            X2 = bp(g) - (2*a*r1 - a) * abs(2*r2*bp(g) - wolves(w,g));

            r1 = rand; r2 = rand;
            X3 = dp(g) - (2*a*r1 - a) * abs(2*r2*dp(g) - wolves(w,g));

            raw = round((X1 + X2 + X3) / 3);
            raw = max(1, min(na, raw));
            nw2(g) = raw;
        end

        nw2 = deduplicate_wolf(nw2, 1:na);
        wolves(w,:) = nw2;

        real_ids = alive_idx(nw2);
        nf = gwo_fitness_aco(real_ids, S.Nodes, sink, P, n);
        if nf < fitness(w)
            fitness(w) = nf;
        end
    end

    [~, si] = sort(fitness);
    ap = wolves(si(1), :);
    bp = wolves(si(min(2,end)), :);
    dp = wolves(si(min(3,end)), :);
end

chi = unique(alive_idx(ap));

[chi, S.Nodes, member_map] = assign_and_prune_CHs(S.Nodes, alive_idx, chi, sink);
S.ch_indices = chi;

if isempty(chi)
    S.num_alive = sum([S.Nodes.alive]);
    return;
end

non = setdiff(alive_idx, chi);
for i = non
    ch = S.Nodes(i).cluster_head;
    if ch > 0 && ch <= n && S.Nodes(ch).alive
        d = dist2D(S.Nodes(i), S.Nodes(ch));
        S.Nodes(i).E  = max(0, S.Nodes(i).E  - tx_energy(P, d));
        S.Nodes(ch).E = max(0, S.Nodes(ch).E - rx_energy(P));

        if S.Nodes(i).E <= 0
            S.Nodes(i).alive = 0;
        end
        if S.Nodes(ch).E <= 0
            S.Nodes(ch).alive = 0;
        end
    end
end

alive_idx = find([S.Nodes.alive] == 1);
chi = intersect(chi, alive_idx, 'stable');
S.ch_indices = chi;

if isempty(chi)
    S.num_alive = sum([S.Nodes.alive]);
    return;
end

[chi, S.Nodes, member_map] = assign_and_prune_CHs(S.Nodes, alive_idx, chi, sink);
S.ch_indices = chi;

if isempty(chi)
    S.num_alive = sum([S.Nodes.alive]);
    return;
end

S.pheromone = max(S.pheromone * (1 - S.ACO.rho), 1e-6);

cs = zeros(1, n);
for cidx = 1:numel(chi)
    c = chi(cidx);
    if isKey(member_map, c)
        cs(c) = numel(member_map(c));
    else
        cs(c) = 0;
    end
end

for ci = 1:numel(chi)
    src = chi(ci);
    if ~S.Nodes(src).alive
        continue;
    end

    d_direct = sqrt((S.Nodes(src).x - sink.x)^2 + (S.Nodes(src).y - sink.y)^2);
    direct_cost = tx_energy(P, d_direct) + P.E_agg * P.k_bits * cs(src);
    best_path = [src, n+1];
    best_cost = direct_cost;

    for ant = 1:S.ACO.num_ants
        p = aco_build_path(src, n+1, chi, S.Nodes, sink, S.pheromone, S.ACO, P, n);
        if ~isempty(p)
            c2 = aco_path_cost(p, S.Nodes, sink, P, n);
            if c2 < best_cost
                best_cost = c2;
                best_path = p;
            end
        end
    end

    if ~isempty(best_path) && numel(best_path) >= 2
        valid_nodes = best_path(best_path <= n);
        if isempty(valid_nodes)
            avg_load = 0;
        else
            avg_load = mean(arrayfun(@(x) S.Nodes(x).relay_load, valid_nodes));
        end

        dep = S.ACO.Q / (best_cost * (1 + avg_load + eps));

        for s = 1:numel(best_path)-1
            u = best_path(s);
            v = best_path(s+1);
            S.pheromone(u,v) = S.pheromone(u,v) + dep;
            S.pheromone(v,u) = S.pheromone(v,u) + dep;
        end

        for s = 1:numel(best_path)-1
            u = best_path(s);
            if u > n
                continue;
            end

            v = best_path(s+1);

            if v <= n
                S.Nodes(v).relay_load = S.Nodes(v).relay_load + 1;
            end

            ea = P.E_agg * P.k_bits * cs(u);

            if v > n
                d = sqrt((S.Nodes(u).x - sink.x)^2 + (S.Nodes(u).y - sink.y)^2);
            else
                d = dist2D(S.Nodes(u), S.Nodes(v));
            end

            S.Nodes(u).E = max(0, S.Nodes(u).E - tx_energy(P, d) - ea);

            if v <= n
                S.Nodes(v).E = max(0, S.Nodes(v).E - rx_energy(P));
            end

            if S.Nodes(u).E <= 0
                S.Nodes(u).alive = 0;
            end
            if v <= n && S.Nodes(v).E <= 0
                S.Nodes(v).alive = 0;
            end
        end
    end
end

S.num_alive = sum([S.Nodes.alive]);
end

function [chi, Nodes, member_map] = assign_and_prune_CHs(Nodes, alive_idx, chi, sink)
member_map = containers.Map('KeyType','double','ValueType','any');

if isempty(chi)
    for i = alive_idx
        Nodes(i).cluster_head = -1;
    end
    return;
end

if numel(alive_idx) == 1
    i = alive_idx(1);
    chi = i;
    Nodes(i).cluster_head = i;
    member_map(i) = [];
    return;
end

MAX_ITERS = 50;
iter_count = 0;
changed = true;

while changed && iter_count < MAX_ITERS
    changed = false;
    iter_count = iter_count + 1;

    for i = alive_idx
        Nodes(i).cluster_head = -1;
    end

    member_map = containers.Map('KeyType','double','ValueType','any');
    for c = chi
        Nodes(c).cluster_head = c;
        member_map(c) = [];
    end

    non = setdiff(alive_idx, chi);

    for i = non
        bd = inf;
        bc = chi(1);
        for c = chi
            d = dist2D(Nodes(i), Nodes(c));
            if d < bd
                bd = d;
                bc = c;
            end
        end
        Nodes(i).cluster_head = bc;
        tmp = member_map(bc);
        tmp(end+1) = i;
        member_map(bc) = tmp;
    end

    keep_mask = true(size(chi));
    zero_member_count = 0;

    for kk = 1:numel(chi)
        c = chi(kk);
        if ~isKey(member_map, c) || isempty(member_map(c))
            keep_mask(kk) = false;
            zero_member_count = zero_member_count + 1;
        end
    end

    if zero_member_count == 0
        break;
    end

    if ~any(keep_mask)
        ds = inf(1, numel(chi));
        for kk = 1:numel(chi)
            c = chi(kk);
            ds(kk) = sqrt((Nodes(c).x - sink.x)^2 + (Nodes(c).y - sink.y)^2);
        end
        [~, best_idx] = min(ds);
        keep_mask(best_idx) = true;
    end

    chi = chi(keep_mask);
    changed = true;
end

member_map = containers.Map('KeyType','double','ValueType','any');

for i = alive_idx
    Nodes(i).cluster_head = -1;
end

for c = chi
    Nodes(c).cluster_head = c;
    member_map(c) = [];
end

non = setdiff(alive_idx, chi);

for i = non
    bd = inf;
    bc = chi(1);
    for c = chi
        d = dist2D(Nodes(i), Nodes(c));
        if d < bd
            bd = d;
            bc = c;
        end
    end

    Nodes(i).cluster_head = bc;
    tmp = member_map(bc);
    tmp(end+1) = i;
    member_map(bc) = tmp;
end
end

function f = gwo_fitness_aco(cs, Nodes, sink, P, n)
w1 = 0.25;
w2 = 0.25;
w3 = 0.35;
w4 = 0.15;

cs = unique(cs);
k = numel(cs);

di  = 0;
cnt = 0;

for i = 1:n
    if ~Nodes(i).alive
        continue;
    end
    md = inf;
    for c = cs
        d = sqrt((Nodes(i).x - Nodes(c).x)^2 + (Nodes(i).y - Nodes(c).y)^2);
        if d < md
            md = d;
        end
    end
    di  = di + md;
    cnt = cnt + 1;
end

if cnt > 0
    di = di / cnt;
end

db = 0;
for c = cs
    db = db + sqrt((Nodes(c).x - sink.x)^2 + (Nodes(c).y - sink.y)^2);
end
db = db / k;

min_db = inf;
for c = cs
    d_s = sqrt((Nodes(c).x - sink.x)^2 + (Nodes(c).y - sink.y)^2);
    if d_s < min_db
        min_db = d_s;
    end
end

tE = sum(arrayfun(@(c) Nodes(c).E, cs));
if tE <= 0
    tE = 1e-9;
end

Dm = sqrt(P.xm^2 + P.ym^2);
f  = w1*(di/Dm) + w2*(db/Dm) + w3*((1/tE)*P.E0) + w4*(min_db/Dm);
end

function wolf = deduplicate_wolf(wolf, valid_pool)
[uw, ia] = unique(wolf, 'stable');
if numel(uw) == numel(wolf)
    return;
end

used = uw;
dp   = setdiff(1:numel(wolf), ia);
ca   = setdiff(valid_pool, used);

for p = dp
    if isempty(ca)
        break;
    end
    wolf(p) = ca(1);
    ca(1) = [];
    used(end+1) = wolf(p); %#ok<AGROW>
end
end

function path = aco_build_path(src, dest, chi, Nodes, sink, ph, ACO, P, n)
path = [src];
vis  = [src];
cur  = src;
mh   = numel(chi) + 2;

for hop = 1:mh
    ca = union(setdiff(chi, vis), dest);
    if isempty(ca)
        break;
    end

    pr = zeros(1, numel(ca));
    for ci = 1:numel(ca)
        nx = ca(ci);
        if nx == dest
            d   = sqrt((Nodes(cur).x - sink.x)^2 + (Nodes(cur).y - sink.y)^2);
            eta = 1 / (d + 1e-9);
        else
            d   = dist2D(Nodes(cur), Nodes(nx));
            if Nodes(nx).alive
                eta = Nodes(nx).E / (d + 1e-9);
            else
                eta = 0;
            end
        end

        tau   = ph(cur, nx);
        pr(ci)= (tau^ACO.alpha) * (eta^ACO.beta);
    end

    tot = sum(pr);
    if tot == 0
        nx = dest;
    else
        pr = pr / tot;
        cm = cumsum(pr);
        idx = find(cm >= rand(), 1, 'first');
        if isempty(idx)
            idx = numel(ca);
        end
        nx = ca(idx);
    end

    path(end+1) = nx; %#ok<AGROW>
    vis(end+1)  = nx; %#ok<AGROW>
    cur         = nx;

    if nx == dest
        break;
    end
end

if path(end) ~= dest
    path(end+1) = dest;
end
end

function cost = aco_path_cost(path, Nodes, sink, P, n)
cost = 0;

for s = 1:numel(path)-1
    u = path(s);
    v = path(s+1);

    if v > n
        d = sqrt((Nodes(u).x - sink.x)^2 + (Nodes(u).y - sink.y)^2);
    else
        d = dist2D(Nodes(u), Nodes(v));
    end

    if v <= n
        load_penalty = 0.5 * Nodes(v).relay_load;
    else
        load_penalty = 0;
    end

    cost = cost + tx_energy(P, d) + rx_energy(P) + load_penalty;
end
end


function d = dist2D(A, B)
d = sqrt((A.x - B.x)^2 + (A.y - B.y)^2);
end

function e = tx_energy(P, d)
d0 = sqrt(P.E_fs / P.E_mp);
if d < d0
    e = P.E_elec * P.k_bits + P.E_fs * P.k_bits * d^2;
else
    e = P.E_elec * P.k_bits + P.E_mp * P.k_bits * d^4;
end
end

function e = rx_energy(P)
e = P.E_elec * P.k_bits;
end