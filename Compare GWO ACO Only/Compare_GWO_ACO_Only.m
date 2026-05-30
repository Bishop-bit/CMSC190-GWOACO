function Compare_GWO_ACO_Only()
% Runs the GWO-ACO ONLY (no comparison with other protocols yet)
% for a specified number of nodes across multiple runs

% Clear workspace, close all figures, and clean command window 
% (useful for after runs, to just clear the workspace)
% Ensures no leftover variables or plots interfere with this run
clear; close all; clc;

% Seed the random number generator using the system clock
% This makes the overall simulation environment random each time
rng('shuffle');

% Define which node counts to simulate
% Can be extended to [50, 100, 200, ...] to test multiple network sizes
% Use 50:50:500 to generate from 50 to 500 with 50 increments.
node_counts = 50;

% Count how many different node sizes will be tested
% Here nN = 1 since only one value (50) is in node_counts
nN = numel(node_counts);

% Number of independent runs per node count
NUM_RUNS = 5;

% Checkpoint interval: log data to CSV every 100 rounds
% Avoids writing every single round to keep file sizes manageable
CHKPT = 100;

% Initial energy of every sensor node in Joules
% Used to override whatever Init_Network sets
E0 = 0.5;

% Maximum number of simulation rounds
RMAX_FIXED = 5000;


% OUTPUT FOLDER SETUP
% Build the full path for the output folder inside the current directory
out_dir = fullfile(pwd, 'simulation_output_gwo_aco_only');

if exist(out_dir, 'dir')
    % Folder already exists: test if it is actually writable
    % by trying to create a temporary probe file inside it
    probe = fullfile(out_dir, '.write_test');
    [ftest, ~] = fopen(probe, 'wt');

    if ftest == -1
        % Could not write: folder may be read-only or locked
        % Fall back to the system temp directory instead
        out_dir = fullfile(tempdir, 'simulation_output_gwo_aco_only');
        if ~exist(out_dir, 'dir'), mkdir(out_dir); end
        warning('Output folder not writable. Using: %s', out_dir);
    else
        % Write succeeded: close and delete the probe file, it was just a test
        fclose(ftest);
        delete(probe);
    end
else
    % Folder does not exist yet: create it now
    mkdir(out_dir);
end

% STORAGE INITIALIZATION
% Master cell array collecting every logged row across all runs and node counts
% Written to the combined master CSV at the very end
all_rows = {};

% Matrix storing the final round (lifetime) of each run
% Size: nN rows x NUM_RUNS columns: one value per node count per run
LT_all = zeros(nN, NUM_RUNS);

% Cell arrays storing full time-series data for each run
% Each cell holds a vector with one value per round
% Lengths differ between runs since networks die at different rounds
RE_runs   = cell(nN, NUM_RUNS);  % Residual energy (%) per round
Dead_runs = cell(nN, NUM_RUNS);  % Dead node count per round
Sur_runs  = cell(nN, NUM_RUNS);  % Survival rate (%) per round

% OUTER LOOP: iterate over each node count
for ki = 1:nN
    % Get the current node count for this iteration
    n = node_counts(ki);

    % Print progress header to command window
    fprintf('\n=== N=%d ===\n', n);
    fprintf('  E0=%.4fJ  rmax=%d\n', E0, RMAX_FIXED);

    % Temporary cell array collecting rows just for this node count
    % Written to its own per-N CSV file after all runs for this N finish
    nc_rows = {};

    % INNER LOOP: repeat simulation NUM_RUNS times for this node count
    for run = 1:NUM_RUNS

        % Generate a random seed between 1 and 1,000,000 for this run
        % Passed to Init_Network so each run has a different node layout
        % but can be reproduced later if the seed is saved
        seed = randi(1e6);
        fprintf('  Run %d/%d (seed=%d)\n', run, NUM_RUNS, seed);

        % Initialize the network: place nodes, set parameters, define sink
        [Nodes0, sink, P] = Init_Network(n, seed);

        % Ensures all runs use exactly RMAX_FIXED rounds and E0 energy
        P.rmax = RMAX_FIXED;
        P.E0   = E0;

        % Explicitly reset each node's key fields to clean starting values
        % Necessary in case Init_Network left any non-default state
        % Redundant but necessary just in case.
        for i = 1:n
            Nodes0(i).E            = E0;   % Full battery
            Nodes0(i).alive        = 1;    % All nodes start alive
            Nodes0(i).cluster_head = -1;   % No cluster assigned yet
            Nodes0(i).relay_load   = 0;    % No relay activity yet
        end

        % Initialize the GWO wolf pack and ACO colony state
        % S is the main state struct passed through every simulation round
        S = init_gwo_aco_state_clean(Nodes0, sink, P, n);

        % Pre-allocate vectors to store metrics for every round
        % Pre-allocation is faster than dynamically growing arrays in MATLAB
        re   = zeros(1, P.rmax);  % Residual energy (%) per round
        dead = zeros(1, P.rmax);  % Dead node count per round
        sur  = zeros(1, P.rmax);  % Survival rate (%) per round

        % Assume the network survives all rounds: update if it dies early
        fin = P.rmax;

        % Tracks the last round written to CSV to prevent duplicate entries
        % when both a checkpoint and a network-death trigger at the same round
        last_chkpt_logged = 0;

        % ROUND LOOP: execute one round at a time up to rmax
        for r = 1:P.rmax

            if S.active
                % Run one full round: GWO selects CHs, ACO routes data,
                % energy is deducted, and node alive statuses are updated
                fprintf('    entering round %d...\n', r);
                S = step_gwo_aco(S, sink, P, n);
                fprintf('    finished round %d\n', r);

                % Check if all nodes have died after this round
                if S.num_alive == 0
                    S.active = false;      % Stop running further rounds
                    S.final_round = r;     % Record when the network died
                    fin = r;               % Save lifetime for this run
                end
            end

            % Count alive nodes and sum remaining energy across all nodes
            alive = sum([S.Nodes.alive]);
            e     = sum([S.Nodes.E]);

            % Calculate residual energy as % of total starting energy
            % Formula: current total energy / max possible energy x 100
            re(r) = (e / (n * P.E0)) * 100;

            % Count how many nodes have died so far
            dead(r) = n - alive;

            % Calculate what percentage of nodes are still alive
            sur(r) = (alive / n) * 100;

            % Determine if this round should be logged to CSV
            is_chkpt = (mod(r, CHKPT) == 0);  % True every 100 rounds
            is_end   = ~S.active;              % True when network just died

            if is_chkpt || is_end
                % Avoid logging the same round twice if both conditions are true
                if r ~= last_chkpt_logged
                    % Build one data row: node count, run, round, and all metrics
                    row = {n, run, r, alive, sur(r), dead(r), 100-sur(r), re(r)};
                    % Note for ok<AGROW>: MATLAB will tell us that the
                    % variable appears to be growing inside a loop.
                    % ok<AGROW> tells MATLAB that I know it is growing. I
                    % am doing it on purpose, stop warning me.
                    nc_rows{end+1}  = row; %#ok<AGROW>
                    all_rows{end+1} = row; %#ok<AGROW>
                    last_chkpt_logged = r;
                end
            end

            % Exit the round loop early if the network has died
            if ~S.active
                break;
            end
        end

        % Log a special lifetime row with round = -1 as a marker
        % This separates lifetime entries from checkpoint entries in the CSV
        lt_row = {n, run, -1, fin, NaN, NaN, NaN, NaN};
        nc_rows{end+1}  = lt_row; %#ok<AGROW>
        all_rows{end+1} = lt_row; %#ok<AGROW>

        % Trim the pre-allocated vectors down to only rounds that happened
        % Removes trailing zeros from rounds that were never executed
        re(fin+1:end)   = [];
        dead(fin+1:end) = [];
        sur(fin+1:end)  = [];

        % Store this run's lifetime and time-series data for later averaging
        LT_all(ki,run)    = fin;
        RE_runs{ki,run}   = re;
        Dead_runs{ki,run} = dead;
        Sur_runs{ki,run}  = sur;

        fprintf('    GWO-ACO LT=%d\n', fin);

        % PER-RUN FIGURE: three subplots showing metrics over rounds
        figure('Name', sprintf('GWO-ACO | N=%d Run %d/%d', n, run, NUM_RUNS), ...
            'NumberTitle','off','Position',[50 50 1100 800]);

        % Checkpoint positions along the x-axis for marker dots
        chk = CHKPT:CHKPT:numel(sur);

        % --- Subplot 1: Node Survival Rate ---
        subplot(3,1,1);
        plot(1:numel(sur), sur, 'b-', 'LineWidth', 1.5); hold on;
        if ~isempty(chk)
            % Place square markers at each checkpoint round
            plot(chk, sur(chk), 'bs', 'MarkerFaceColor', 'b', 'MarkerSize', 6);
        end
        % Vertical dotted line marking the round the network died
        plot([fin fin], [0 100], 'k:', 'LineWidth', 1.5);
        xlabel('Round'); ylabel('Alive (%)');
        title(sprintf('Node Survival Rate (LT=%d)', fin));
        grid on; ylim([0 100]);

        % --- Subplot 2: Dead Node Count ---
        subplot(3,1,2);
        plot(1:numel(dead), dead, 'r-', 'LineWidth', 1.5); hold on;
        if ~isempty(chk)
            % Place circle markers at each checkpoint round
            plot(chk, dead(chk), 'ro', 'MarkerFaceColor', 'r', 'MarkerSize', 6);
        end
        plot([fin fin], [0 n], 'k:', 'LineWidth', 1.5);
        xlabel('Round'); ylabel('Dead Nodes');
        title('Dead Node Count');
        grid on; ylim([0 n]);

        % --- Subplot 3: Residual Energy ---
        subplot(3,1,3);
        plot(1:numel(re), re, 'g-', 'LineWidth', 1.5); hold on;
        if ~isempty(chk)
            % Place diamond markers at each checkpoint round
            plot(chk, re(chk), 'gd', 'MarkerFaceColor', 'g', 'MarkerSize', 6);
        end
        plot([fin fin], [0 100], 'k:', 'LineWidth', 1.5);
        xlabel('Round'); ylabel('Residual Energy (%)');
        title('Residual Energy');
        grid on; ylim([0 100]);

        % Overall title for the figure showing run details
        sgtitle(sprintf('GWO-ACO Only | N=%d | Run %d/%d | Seed=%d', ...
            n, run, NUM_RUNS, seed), 'FontWeight', 'bold');
        drawnow;
    end

    % Write all checkpoint rows for this node count to its own CSV file
    csv_path = fullfile(out_dir, sprintf('N%d_detailed.csv', n));
    write_csv_gwo_only(nc_rows, csv_path);

    % Plot the averaged curves across all 5 runs for this node count
    plot_averaged_gwo_only(ki, n, NUM_RUNS, CHKPT, RE_runs, Dead_runs, Sur_runs, LT_all);
end

% POST-SIMULATION OUTPUTS
% Plot lifetime vs N and per-metric summary subplots for all node counts
plot_global_summary_gwo_only(node_counts, nN, NUM_RUNS, LT_all, RE_runs, Dead_runs, Sur_runs);

% Plot overlaid comparison figures — one per metric — across all N values
plot_allN_comparison_gwo_only(node_counts, nN, NUM_RUNS, RE_runs, Dead_runs, Sur_runs);

% Write every single logged row from all runs and node counts to one master CSV
write_csv_gwo_only(all_rows, fullfile(out_dir, 'ALL_NODES_master.csv'));

fprintf('\nDone. CSVs in: %s\n', out_dir);
end

% GWO-ACO STATE INITIALIZATION
% Sets up the wolf pack, ACO colony, and pheromone matrix
% before the simulation rounds begin
function S = init_gwo_aco_state_clean(Nodes, sink, P, n)

% Copy the initialized node array into the state struct
S.Nodes = Nodes;

% Network starts active: set to false when all nodes die
S.active = true;

% Assume the network survives all rounds: updated when it dies early
S.final_round = P.rmax;

% Placeholder for survival tracking: populated during simulation
S.survival = [];

% Residual energy tracker: starts at zero, updated each round
S.res_energy = 0;

% Number of wolves in the GWO pack
% More wolves = more diverse search but higher computation per round
S.num_wolves = 50;

% Number of GWO iterations per round
% Controls how many refinement steps the wolf pack takes each round
S.max_iter = 25;


% ACO PARAMETERS
% Number of ants constructing paths each round
S.ACO.num_ants = 15;

% Alpha: controls how strongly ants follow pheromone trails
% Higher alpha = ants stick more to previously good paths
S.ACO.alpha = 1;

% Beta: controls how strongly ants prefer closer / higher energy nodes
% Higher beta = ants prioritize distance and energy over pheromone
S.ACO.beta = 4;

% Rho: pheromone evaporation rate (0 to 1)
% 0.2 means 20% of pheromone evaporates each round
% Prevents the algorithm from permanently locking onto one path
S.ACO.rho = 0.2;

% Q: pheromone deposit constant
% Controls how much pheromone is added to a winning path
S.ACO.Q = 1;

% Initialize pheromone matrix to all ones
% Size (n+1) x (n+1): the extra row/column represents the sink (node n+1)
% All paths start with equal pheromone so no path is initially preferred
S.pheromone = ones(n+1, n+1);

% CH indices are computed dynamically every round based on alive nodes
% Do not set a fixed CH count here: leave empty until first round
S.ch_indices = [];

% All nodes are alive at the start
S.num_alive = n;
end


% GWO-ACO STEP — executes one full simulation round
% Phases: CH selection (GWO) -> cluster assignment ->
%         intra-cluster TX -> inter-cluster routing (ACO)
function S = step_gwo_aco(S, sink, P, n)

% Find all nodes that are currently alive
alive_idx = find([S.Nodes.alive] == 1);

% Reset relay load counter for all alive nodes at the start of each round
% Ensures load counts are fresh and not carried over from last round
for i = alive_idx
    S.Nodes(i).relay_load = 0;
end

% Update the alive node count
S.num_alive = numel(alive_idx);

% If no nodes are alive, exit immediately — nothing to do
if S.num_alive == 0
    return;
end

% EDGE CASE: Only one node left alive
% No need to cluster: just transmit directly to the sink
if S.num_alive == 1
    i = alive_idx(1);

    % Calculate distance from the last surviving node to the sink
    d = sqrt((S.Nodes(i).x - sink.x)^2 + (S.Nodes(i).y - sink.y)^2);

    % Deduct transmission energy: node sends directly to base station
    S.Nodes(i).E = max(0, S.Nodes(i).E - tx_energy(P, d));

    % Check if the node has depleted its energy and mark it dead
    if S.Nodes(i).E <= 0
        S.Nodes(i).alive = 0;
    end

    S.num_alive = sum([S.Nodes.alive]);
    S.ch_indices = [];
    return;
end

% Number of alive nodes used throughout this round
na = S.num_alive;

% DYNAMIC OPTIMAL CH COUNT (k_opt)
% Computes the ideal number of cluster heads based on alive nodes
% and their average distance to the base station

% Sum of distances from all alive nodes to the sink
d_sum = 0;
for ii = alive_idx
    d_sum = d_sum + sqrt((S.Nodes(ii).x - sink.x)^2 + (S.Nodes(ii).y - sink.y)^2);
end

% Average distance from alive nodes to the sink
dBS_alive = d_sum / na;

% Optimal CH count formula derived from energy balancing
% Balances intra-cluster and inter-cluster energy costs
k_opt = sqrt(na/(2*pi)) * sqrt(P.E_fs/P.E_mp) * (P.xm/(dBS_alive^2));

% Enforce practical bounds: at least 4% and at most 12% of alive nodes as CHs
% Prevents too few CHs (overloaded) or too many CHs (wasteful)
k_min = ceil(0.04 * na);
k_max = ceil(0.12 * na);

% Round k_opt to nearest integer and clamp within bounds
k = round(k_opt);
k = max(k_min, min(k, k_max));

% Final safety clamp: k must be at least 1 and at most na
k = max(1, min(k, na));

fprintf('alive=%d, k=%d\n', na, k);


% GWO INITIALIZATION
% Each wolf represents one candidate set of k cluster heads
% Wolves operate in alive-index space (positions 1..na)
% Initialize each wolf with k randomly chosen alive node indices
wolves = zeros(S.num_wolves, k);
for w = 1:S.num_wolves
    wolves(w,:) = randperm(na, k);
end

% Evaluate the fitness of each wolf's CH selection
fitness = zeros(1, S.num_wolves);
for w = 1:S.num_wolves
    real_ids = alive_idx(wolves(w,:));  % Convert local indices to real node IDs
    fitness(w) = gwo_fitness_aco(real_ids, S.Nodes, sink, P, n);
end

% Sort wolves by fitness: lower is better
[~, si] = sort(fitness);

% The three best wolves become alpha (best), beta (2nd), delta (3rd)
% These guide the rest of the pack toward better solutions
ap = wolves(si(1), :);
bp = wolves(si(min(2,end)), :);
dp = wolves(si(min(3,end)), :);


% GWO ITERATION LOOP
% Wolves update their positions guided by alpha, beta, delta
for iter = 1:S.max_iter

    % Linearly decrease 'a' from 2 to 0 as iterations progress
    % Controls the balance between exploration (high a) and exploitation (low a)
    a = 2 - 2*(iter / S.max_iter);

    for w = 1:S.num_wolves
        nw2 = zeros(1, k);

        for g = 1:k
            % Update position guided by alpha wolf
            r1 = rand; r2 = rand;
            X1 = ap(g) - (2*a*r1 - a) * abs(2*r2*ap(g) - wolves(w,g));

            % Update position guided by beta wolf
            r1 = rand; r2 = rand;
            X2 = bp(g) - (2*a*r1 - a) * abs(2*r2*bp(g) - wolves(w,g));

            % Update position guided by delta wolf
            r1 = rand; r2 = rand;
            X3 = dp(g) - (2*a*r1 - a) * abs(2*r2*dp(g) - wolves(w,g));

            % Average the three guided positions and round to valid node index
            raw = round((X1 + X2 + X3) / 3);

            % Clamp to valid alive-index range [1, na]
            raw = max(1, min(na, raw));
            nw2(g) = raw;
        end

        % Remove duplicate node selections in this wolf's position
        % A node cannot be its own CH twice in the same wolf
        nw2 = deduplicate_wolf(nw2, 1:na);
        wolves(w,:) = nw2;

        % Re-evaluate fitness with the updated position
        real_ids = alive_idx(nw2);
        nf = gwo_fitness_aco(real_ids, S.Nodes, sink, P, n);

        % Only update if the new position is actually better
        if nf < fitness(w)
            fitness(w) = nf;
        end
    end

    % Re-rank wolves and update alpha, beta, delta after each iteration
    [~, si] = sort(fitness);
    ap = wolves(si(1), :);
    bp = wolves(si(min(2,end)), :);
    dp = wolves(si(min(3,end)), :);
end

% Convert the best wolf's indices to real node IDs
% unique() removes any accidental duplicates
chi = unique(alive_idx(ap));

% CLUSTER ASSIGNMENT AND PRUNING
% Assign every non-CH node to its nearest CH
% Remove any CH that ends up with zero members
[chi, S.Nodes, member_map] = assign_and_prune_CHs(S.Nodes, alive_idx, chi, sink);
S.ch_indices = chi;

% If no CHs remain after pruning, update alive count and exit
if isempty(chi)
    S.num_alive = sum([S.Nodes.alive]);
    return;
end

% INTRA-CLUSTER TRANSMISSION
% Every non-CH node sends its data to its assigned cluster head
non = setdiff(alive_idx, chi);  % All alive nodes that are NOT cluster heads

for i = non
    ch = S.Nodes(i).cluster_head;

    % Only transmit if the assigned CH is valid and still alive
    if ch > 0 && ch <= n && S.Nodes(ch).alive
        d = dist2D(S.Nodes(i), S.Nodes(ch));

        % Deduct TX energy from the member node
        S.Nodes(i).E  = max(0, S.Nodes(i).E  - tx_energy(P, d));

        % Deduct RX energy from the cluster head receiving the data
        S.Nodes(ch).E = max(0, S.Nodes(ch).E - rx_energy(P));

        % Mark nodes as dead if energy is fully depleted
        if S.Nodes(i).E <= 0
            S.Nodes(i).alive = 0;
        end
        if S.Nodes(ch).E <= 0
            S.Nodes(ch).alive = 0;
        end
    end
end

% Rebuild alive index after intra-cluster deaths
alive_idx = find([S.Nodes.alive] == 1);

% Remove any CHs that died during intra-cluster transmission
chi = intersect(chi, alive_idx, 'stable');
S.ch_indices = chi;

if isempty(chi)
    S.num_alive = sum([S.Nodes.alive]);
    return;
end

% Re-run cluster assignment after deaths to keep member map accurate
[chi, S.Nodes, member_map] = assign_and_prune_CHs(S.Nodes, alive_idx, chi, sink);
S.ch_indices = chi;

if isempty(chi)
    S.num_alive = sum([S.Nodes.alive]);
    return;
end



% ACO INTER-CLUSTER ROUTING
% Each CH uses ACO to find the best path to the sink
% Path can be direct (CH -> sink) or multi-hop (CH -> CH -> sink)

% Evaporate pheromone on all edges by factor (1 - rho)
% Prevents any one path from permanently dominating
% Clamped to 1e-6 minimum to avoid zero pheromone dead zones
S.pheromone = max(S.pheromone * (1 - S.ACO.rho), 1e-6);

% Count how many member nodes each CH is aggregating data for
cs = zeros(1, n);
for cidx = 1:numel(chi)
    c = chi(cidx);
    if isKey(member_map, c)
        cs(c) = numel(member_map(c));  % Number of members in this CH's cluster
    else
        cs(c) = 0;
    end
end

% Process each cluster head
for ci = 1:numel(chi)
    src = chi(ci);

    % Skip this CH if it died during earlier processing
    if ~S.Nodes(src).alive
        continue;
    end

    % --- Baseline: direct transmission from this CH to the sink ---
    d_direct = sqrt((S.Nodes(src).x - sink.x)^2 + (S.Nodes(src).y - sink.y)^2);

    % Direct cost = TX energy + aggregation energy for all members
    direct_cost = tx_energy(P, d_direct) + P.E_agg * P.k_bits * cs(src);

    % Start with direct path as the best known option
    best_path = [src, n+1];   % n+1 represents the sink node
    best_cost = direct_cost;

    % --- ACO tries to beat the direct path with a multi-hop path ---
    for ant = 1:S.ACO.num_ants
        % Each ant constructs one candidate path from src to sink
        p = aco_build_path(src, n+1, chi, S.Nodes, sink, S.pheromone, S.ACO, P, n);

        if ~isempty(p)
            % Evaluate the cost of this ant's path
            c2 = aco_path_cost(p, S.Nodes, sink, P, n);

            % Keep this path if it is cheaper than the current best
            if c2 < best_cost
                best_cost = c2;
                best_path = p;
            end
        end
    end

    % --- Pheromone deposit on the winning path ---
    if ~isempty(best_path) && numel(best_path) >= 2

        % Calculate average relay load on the winning path
        % Used to reduce pheromone deposit on overloaded paths
        valid_nodes = best_path(best_path <= n);
        if isempty(valid_nodes)
            avg_load = 0;
        else
            avg_load = mean(arrayfun(@(x) S.Nodes(x).relay_load, valid_nodes));
        end

        % Deposit amount is inversely proportional to cost and load
        % Better (cheaper, less loaded) paths get more pheromone
        dep = S.ACO.Q / (best_cost * (1 + avg_load + eps));

        % Add pheromone to each edge in the winning path (bidirectional)
        for s = 1:numel(best_path)-1
            u = best_path(s);
            v = best_path(s+1);
            S.pheromone(u,v) = S.pheromone(u,v) + dep;
            S.pheromone(v,u) = S.pheromone(v,u) + dep;
        end

        % --- Energy deduction along the winning path ---
        for s = 1:numel(best_path)-1
            u = best_path(s);

            % Skip if u is the sink (index > n) — sink has no energy to deduct
            if u > n
                continue;
            end

            v = best_path(s+1);

            % Increment relay load on the receiving node (if it is a real node)
            if v <= n
                S.Nodes(v).relay_load = S.Nodes(v).relay_load + 1;
            end

            % Aggregation energy cost for this CH based on its cluster size
            ea = P.E_agg * P.k_bits * cs(u);

            % Calculate distance to the next hop
            if v > n
                % Next hop is the sink — use direct distance
                d = sqrt((S.Nodes(u).x - sink.x)^2 + (S.Nodes(u).y - sink.y)^2);
            else
                % Next hop is another CH node
                d = dist2D(S.Nodes(u), S.Nodes(v));
            end

            % Deduct TX energy + aggregation energy from the sending node
            S.Nodes(u).E = max(0, S.Nodes(u).E - tx_energy(P, d) - ea);

            % Deduct RX energy from the receiving node if it is a real node
            if v <= n
                S.Nodes(v).E = max(0, S.Nodes(v).E - rx_energy(P));
            end

            % Mark nodes as dead if energy fully depleted
            if S.Nodes(u).E <= 0
                S.Nodes(u).alive = 0;
            end
            if v <= n && S.Nodes(v).E <= 0
                S.Nodes(v).alive = 0;
            end
        end
    end
end

% Update total alive count after all transmissions this round
S.num_alive = sum([S.Nodes.alive]);
end


% ASSIGN AND PRUNE CLUSTER HEADS
% Assigns every non-CH alive node to its nearest CH
% Iteratively removes CHs that have no members assigned to them
function [chi, Nodes, member_map] = assign_and_prune_CHs(Nodes, alive_idx, chi, sink)

% Initialize empty member map — keys are CH node IDs, values are member lists
member_map = containers.Map('KeyType','double','ValueType','any');

% If no CHs exist, reset all nodes and return
if isempty(chi)
    for i = alive_idx
        Nodes(i).cluster_head = -1;
    end
    return;
end

% Edge case: only one node alive — it becomes its own CH with no members
if numel(alive_idx) == 1
    i = alive_idx(1);
    chi = i;
    Nodes(i).cluster_head = i;
    member_map(i) = [];
    return;
end

% Limit iterations to prevent infinite loops in edge cases
MAX_ITERS = 50;
iter_count = 0;
changed = true;

while changed && iter_count < MAX_ITERS
    changed = false;
    iter_count = iter_count + 1;

    % Reset all cluster head assignments before re-assigning
    for i = alive_idx
        Nodes(i).cluster_head = -1;
    end

    % Rebuild member map fresh each iteration
    member_map = containers.Map('KeyType','double','ValueType','any');

    % Each CH is assigned to itself
    for c = chi
        Nodes(c).cluster_head = c;
        member_map(c) = [];
    end

    % All alive nodes that are not CHs
    non = setdiff(alive_idx, chi);

    % Assign each non-CH node to its nearest cluster head
    for i = non
        bd = inf;   % Best (minimum) distance found so far
        bc = chi(1); % Best CH found so far

        for c = chi
            d = dist2D(Nodes(i), Nodes(c));
            if d < bd
                bd = d;
                bc = c;
            end
        end

        % Assign node i to the closest CH
        Nodes(i).cluster_head = bc;

        % Add node i to the CH's member list
        tmp = member_map(bc);
        tmp(end+1) = i;
        member_map(bc) = tmp;
    end

    % Check which CHs have zero members
    keep_mask = true(size(chi));
    zero_member_count = 0;

    for kk = 1:numel(chi)
        c = chi(kk);
        if ~isKey(member_map, c) || isempty(member_map(c))
            keep_mask(kk) = false;  % Mark this CH for removal
            zero_member_count = zero_member_count + 1;
        end
    end

    % All CHs have at least one member — no pruning needed, exit loop
    if zero_member_count == 0
        break;
    end

    % Safety: if ALL CHs would be pruned, keep the one closest to the sink
    % This prevents the network from having zero CHs
    if ~any(keep_mask)
        ds = inf(1, numel(chi));
        for kk = 1:numel(chi)
            c = chi(kk);
            ds(kk) = sqrt((Nodes(c).x - sink.x)^2 + (Nodes(c).y - sink.y)^2);
        end
        [~, best_idx] = min(ds);
        keep_mask(best_idx) = true;
    end

    % Remove the zero-member CHs and repeat the assignment
    chi = chi(keep_mask);
    changed = true;
end


% FINAL REBUILD — one clean assignment after pruning is complete
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



% GWO FITNESS FUNCTION
% Scores a candidate set of cluster heads
% Lower score = better CH selection
% Four weighted criteria are evaluated
function f = gwo_fitness_aco(cs, Nodes, sink, P, n)

% Fitness weights — must sum to 1.0
% w3 (energy balance) is weighted highest — energy fairness matters most
w1 = 0.25;  % Weight for average intra-cluster distance
w2 = 0.25;  % Weight for average CH-to-sink distance
w3 = 0.35;  % Weight for CH energy balance (higher = favor energy-rich CHs)
w4 = 0.15;  % Weight for minimum CH-to-sink distance (best direct hop)

% Remove duplicate CH IDs if any
cs = unique(cs);
k = numel(cs);


% CRITERION 1: Average intra-cluster distance (di)
% Measures how compact the clusters are
% Smaller di = nodes are closer to their CH = less TX energy needed
di  = 0;
cnt = 0;

for i = 1:n
    if ~Nodes(i).alive
        continue;
    end

    % Find the distance from node i to its nearest candidate CH
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

% Average the intra-cluster distances
if cnt > 0
    di = di / cnt;
end


% CRITERION 2: Average CH-to-sink distance (db)
% Measures how close the CHs are to the base station on average
% Smaller db = CHs spend less energy transmitting to the sink
db = 0;
for c = cs
    db = db + sqrt((Nodes(c).x - sink.x)^2 + (Nodes(c).y - sink.y)^2);
end
db = db / k;


% CRITERION 4 (used in w4): Minimum CH-to-sink distance
% Identifies the CH best positioned to transmit directly to sink
min_db = inf;
for c = cs
    d_s = sqrt((Nodes(c).x - sink.x)^2 + (Nodes(c).y - sink.y)^2);
    if d_s < min_db
        min_db = d_s;
    end
end


% CRITERION 3 (w3): Total energy of all candidate CHs
% Favors CH sets where nodes still have high residual energy
% Avoids selecting nearly-dead nodes as CHs
tE = sum(arrayfun(@(c) Nodes(c).E, cs));
if tE <= 0
    tE = 1e-9;  % Avoid division by zero if all CHs are nearly dead
end

% Diagonal of the field — used to normalize distances to [0, 1]
Dm = sqrt(P.xm^2 + P.ym^2);

% Combine all four criteria into one weighted fitness score
% Each term is normalized by Dm so distances are on the same scale
f = w1*(di/Dm) + w2*(db/Dm) + w3*((1/tE)*P.E0) + w4*(min_db/Dm);
end


% DEDUPLICATE WOLF POSITIONS
% Ensures no node index appears twice in a wolf's CH selection
% Duplicate CHs waste slots and reduce cluster diversity
function wolf = deduplicate_wolf(wolf, valid_pool)

% Get unique values and their first positions
[uw, ia] = unique(wolf, 'stable');

% If all values are already unique, nothing to fix
if numel(uw) == numel(wolf)
    return;
end

% Track which positions are duplicates
used = uw;
dp   = setdiff(1:numel(wolf), ia);  % Indices of duplicate positions

% Pool of unused valid node indices to replace duplicates with
ca = setdiff(valid_pool, used);

% Replace each duplicate position with an unused node index
for p = dp
    if isempty(ca)
        break;  % No replacements left — leave remaining duplicates as-is
    end
    wolf(p) = ca(1);   % Assign the next unused node
    ca(1) = [];        % Remove it from the available pool
    used(end+1) = wolf(p); %#ok<AGROW>
end
end



% ACO PATH BUILDER
% One ant constructs a path from src (a CH) to dest (the sink)
% Path moves through other CHs using pheromone + heuristic probabilities
function path = aco_build_path(src, dest, chi, Nodes, sink, ph, ACO, P, n)

path = [src];   % Path starts at the source CH
vis  = [src];   % Track visited nodes to avoid cycles
cur  = src;     % Current position of the ant

% Maximum hops allowed — prevents infinite loops in large networks
mh = numel(chi) + 2;

for hop = 1:mh

    % Candidate next nodes = unvisited CHs + the sink (always an option)
    ca = union(setdiff(chi, vis), dest);

    if isempty(ca)
        break;  % No more candidates — path is stuck
    end

    % Calculate selection probability for each candidate
    pr = zeros(1, numel(ca));

    for ci = 1:numel(ca)
        nx = ca(ci);

        if nx == dest
            % Sink is always a valid next hop
            % Desirability = 1 / distance (closer sink = more desirable)
            d   = sqrt((Nodes(cur).x - sink.x)^2 + (Nodes(cur).y - sink.y)^2);
            eta = 1 / (d + 1e-9);
        else
            d = dist2D(Nodes(cur), Nodes(nx));

            % Only consider alive relay nodes
            % Dead nodes get zero desirability — never chosen
            if Nodes(nx).alive
                % Desirability = residual energy / distance
                % High energy + short distance = preferred relay
                eta = Nodes(nx).E / (d + 1e-9);
            else
                eta = 0;
            end
        end

        % Pheromone level on this edge
        tau = ph(cur, nx);

        % Combined probability: pheromone^alpha * desirability^beta
        % alpha and beta control the balance between the two factors
        pr(ci) = (tau^ACO.alpha) * (eta^ACO.beta);
    end

    tot = sum(pr);

    if tot == 0
        % All probabilities are zero — go directly to sink as fallback
        nx = dest;
    else
        % Normalize probabilities and use roulette wheel selection
        pr = pr / tot;
        cm = cumsum(pr);
        idx = find(cm >= rand(), 1, 'first');
        if isempty(idx)
            idx = numel(ca);
        end
        nx = ca(idx);
    end

    % Move the ant to the chosen next node
    path(end+1) = nx; %#ok<AGROW>
    vis(end+1)  = nx; %#ok<AGROW>
    cur         = nx;

    % Stop if the ant has reached the sink
    if nx == dest
        break;
    end
end

% Safety: if the ant never reached the sink, force-append it
% Guarantees every path always ends at the sink
if path(end) ~= dest
    path(end+1) = dest;
end
end



% ACO PATH COST
% Calculates the total energy cost of a complete CH-to-sink path
% Includes TX energy, RX energy, and a penalty for overloaded relay nodes
function cost = aco_path_cost(path, Nodes, sink, P, n)
cost = 0;

for s = 1:numel(path)-1
    u = path(s);
    v = path(s+1);

    % Calculate hop distance
    if v > n
        % Next hop is the sink
        d = sqrt((Nodes(u).x - sink.x)^2 + (Nodes(u).y - sink.y)^2);
    else
        % Next hop is another node
        d = dist2D(Nodes(u), Nodes(v));
    end

    % Add a load penalty if the receiving node is already heavily used
    % Discourages routing through congested relay nodes
    if v <= n
        load_penalty = 0.5 * Nodes(v).relay_load;
    else
        load_penalty = 0;  % Sink has no load penalty
    end

    % Accumulate cost: TX energy + RX energy + load penalty per hop
    cost = cost + tx_energy(P, d) + rx_energy(P) + load_penalty;
end
end



% PLOT AVERAGED RESULTS PER NODE COUNT
% After all runs for one N complete, plots averaged survival,
% dead count, and residual energy with individual run traces behind
function plot_averaged_gwo_only(ki, n, NUM_RUNS, CHKPT, RE_runs, Dead_runs, Sur_runs, LT_all)

    % Inner helper: pad all run vectors to the same length then average
    % Shorter runs are extended by repeating their last value
    function avg = pad_avg(cell_row, nr)
        ml = max(cellfun(@numel, cell_row));
        mat = zeros(nr, ml);
        for rr = 1:nr
            v = cell_row{rr};
            mat(rr,:) = [v, repmat(v(end), 1, ml-numel(v))];
        end
        avg = mean(mat, 1);
    end

    % Compute averaged curves across all runs for this node count
    sa_avg = pad_avg(Sur_runs(ki,:),  NUM_RUNS);
    da_avg = pad_avg(Dead_runs(ki,:), NUM_RUNS);
    ra_avg = pad_avg(RE_runs(ki,:),   NUM_RUNS);

    % Average lifetime across all runs
    lt_avg = mean(LT_all(ki,:));

    % Checkpoint positions for marker dots on the averaged line
    chk = CHKPT:CHKPT:numel(sa_avg);

    figure('Name', sprintf('GWO-ACO | N=%d Averaged (%d runs)', n, NUM_RUNS), ...
        'NumberTitle','off','Position',[120 120 1100 820]);

    % --- Subplot 1: Node Survival Rate ---
    ax1 = subplot(3,1,1); hold(ax1,'on');

    % Plot individual run traces in light blue behind the average
    for rr = 1:NUM_RUNS
        v = Sur_runs{ki,rr};
        plot(ax1, 1:numel(v), v, 'Color', [0.75 0.75 1], 'LineWidth', 0.5, 'HandleVisibility','off');
    end

    % Bold averaged line on top
    plot(ax1, 1:numel(sa_avg), sa_avg, 'b-', 'LineWidth', 2.5, 'DisplayName','Average');
    if ~isempty(chk)
        plot(ax1, chk, sa_avg(chk), 'bs', 'MarkerFaceColor', 'b', 'MarkerSize', 7, 'HandleVisibility','off');
    end

    % Dotted vertical line at the average lifetime
    plot(ax1, [lt_avg lt_avg], [0 100], 'k:', 'LineWidth', 1.5, 'HandleVisibility','off');
    xlabel(ax1, 'Round'); ylabel(ax1, 'Alive (%)');
    title(ax1, sprintf('Node Survival Rate (avg LT = %.0f)', lt_avg));
    grid(ax1, 'on'); ylim(ax1, [0 100]);
    legend(ax1, 'Location', 'southwest');

    % --- Subplot 2: Dead Node Count ---
    ax2 = subplot(3,1,2); hold(ax2,'on');

    % Individual run traces in light red
    for rr = 1:NUM_RUNS
        v = Dead_runs{ki,rr};
        plot(ax2, 1:numel(v), v, 'Color', [1 0.8 0.8], 'LineWidth', 0.5, 'HandleVisibility','off');
    end

    plot(ax2, 1:numel(da_avg), da_avg, 'r-', 'LineWidth', 2.5, 'DisplayName','Average');
    if ~isempty(chk)
        plot(ax2, chk, da_avg(chk), 'ro', 'MarkerFaceColor', 'r', 'MarkerSize', 7, 'HandleVisibility','off');
    end
    plot(ax2, [lt_avg lt_avg], [0 n], 'k:', 'LineWidth', 1.5, 'HandleVisibility','off');
    xlabel(ax2, 'Round'); ylabel(ax2, 'Dead Nodes');
    title(ax2, 'Dead Node Count');
    grid(ax2, 'on'); ylim(ax2, [0 n]);
    legend(ax2, 'Location', 'northwest');

    % --- Subplot 3: Residual Energy ---
    ax3 = subplot(3,1,3); hold(ax3,'on');

    % Individual run traces in light green
    for rr = 1:NUM_RUNS
        v = RE_runs{ki,rr};
        plot(ax3, 1:numel(v), v, 'Color', [0.8 1 0.8], 'LineWidth', 0.5, 'HandleVisibility','off');
    end

    plot(ax3, 1:numel(ra_avg), ra_avg, 'g-', 'LineWidth', 2.5, 'DisplayName','Average');
    if ~isempty(chk)
        plot(ax3, chk, ra_avg(chk), 'gd', 'MarkerFaceColor', 'g', 'MarkerSize', 7, 'HandleVisibility','off');
    end
    plot(ax3, [lt_avg lt_avg], [0 100], 'k:', 'LineWidth', 1.5, 'HandleVisibility','off');
    xlabel(ax3, 'Round'); ylabel(ax3, 'Residual Energy (%)');
    title(ax3, 'Residual Energy');
    grid(ax3, 'on'); ylim(ax3, [0 100]);
    legend(ax3, 'Location', 'northeast');

    sgtitle(sprintf('GWO-ACO Only | N=%d | %d runs averaged', n, NUM_RUNS), ...
        'FontWeight','bold','FontSize',12);
end



% GLOBAL SUMMARY PLOTS
% After all node counts finish, plots lifetime vs N and
% per-metric subplots for every N value tested
function plot_global_summary_gwo_only(node_counts, nN, NUM_RUNS, LT_all, RE_runs, Dead_runs, Sur_runs)

    % Pad and average helper — same as above
    function avg = pad_avg(cell_row, nr)
        ml = max(cellfun(@numel, cell_row));
        mat = zeros(nr, ml);
        for rr = 1:nr
            v = cell_row{rr};
            mat(rr,:) = [v, repmat(v(end), 1, ml-numel(v))];
        end
        avg = mean(mat, 1);
    end

    % Average lifetime across runs for each node count
    LT_mean = mean(LT_all, 2)';

    % Plot network lifetime vs number of nodes
    figure('Name','GWO-ACO Only | Network Lifetime vs N','NumberTitle','off');
    plot(node_counts, LT_mean, 'bo-', 'LineWidth', 2, 'MarkerFaceColor', 'b');
    xlabel('Number of Nodes (N)');
    ylabel('Network Lifetime (Rounds)');
    title(sprintf('GWO-ACO Only | Average Lifetime vs N (%d runs)', NUM_RUNS));
    grid on;

    % Determine subplot grid layout based on how many N values were tested
    nCols = min(nN,5);
    nRows = ceil(nN/nCols);

    % Pre-create one figure per metric for the per-N subplots
    fig_re   = figure('Name','GWO-ACO Only | Residual Energy All N', ...
        'NumberTitle','off','Position',[50 50 300*nCols 300*nRows]);
    fig_dead = figure('Name','GWO-ACO Only | Dead Nodes All N', ...
        'NumberTitle','off','Position',[50 50 300*nCols 300*nRows]);
    fig_sur  = figure('Name','GWO-ACO Only | Survival Rate All N', ...
        'NumberTitle','off','Position',[50 50 300*nCols 300*nRows]);

    for ki = 1:nN
        nc = node_counts(ki);

        % Compute averaged curves for this node count
        ra = pad_avg(RE_runs(ki,:),   NUM_RUNS);
        da = pad_avg(Dead_runs(ki,:), NUM_RUNS);
        sa = pad_avg(Sur_runs(ki,:),  NUM_RUNS);

        % Add one subplot per node count in each figure
        figure(fig_re); subplot(nRows,nCols,ki);
        plot(1:numel(ra), ra, 'g-', 'LineWidth', 1.5);
        xlabel('Round'); ylabel('Energy (%)');
        title(sprintf('N=%d', nc));
        grid on; ylim([0 100]);

        figure(fig_dead); subplot(nRows,nCols,ki);
        plot(1:numel(da), da, 'r-', 'LineWidth', 1.5);
        xlabel('Round'); ylabel('Dead');
        title(sprintf('N=%d', nc));
        grid on; ylim([0 nc]);

        figure(fig_sur); subplot(nRows,nCols,ki);
        plot(1:numel(sa), sa, 'b-', 'LineWidth', 1.5);
        xlabel('Round'); ylabel('Alive (%)');
        title(sprintf('N=%d', nc));
        grid on; ylim([0 100]);
    end

    % Add overall titles to each figure
    figure(fig_re);   sgtitle(sprintf('GWO-ACO Only | Residual Energy (avg %d runs)', NUM_RUNS));
    figure(fig_dead); sgtitle(sprintf('GWO-ACO Only | Dead Node Count (avg %d runs)', NUM_RUNS));
    figure(fig_sur);  sgtitle(sprintf('GWO-ACO Only | Node Survival Rate (avg %d runs)', NUM_RUNS));
end




% CSV WRITER
% Writes all logged checkpoint and lifetime rows to a CSV file
% Falls back to the system temp folder if the target path is not writable
function write_csv_gwo_only(rows, csv_path)

% Create the output folder if it does not exist
out_folder = fileparts(csv_path);
if ~isempty(out_folder) && ~exist(out_folder, 'dir')
    mkdir(out_folder);
end

% Try to open the CSV file for writing
[fid, msg] = fopen(csv_path, 'wt');

if fid == -1
    % Could not open the file — fall back to temp directory
    [~, fname, ext] = fileparts(csv_path);
    fallback_path = fullfile(tempdir, [fname ext]);
    warning('write_csv_gwo_only:Fallback', ...
        'Could not open:\n  %s\nReason: %s\nFalling back to:\n  %s', ...
        csv_path, msg, fallback_path);

    [fid, msg2] = fopen(fallback_path, 'wt');
    if fid == -1
        % Even the fallback failed — throw an error to stop execution
        error('write_csv_gwo_only:FileOpenFailed', ...
            'Could not open fallback CSV:\n%s\nReason: %s', fallback_path, msg2);
    end
    csv_path = fallback_path;
end

% Write the column headers
fprintf(fid, ['Nodes,Run,Round,' ...
    'GWO_ACO_Alive,GWO_ACO_Alive_Pct,GWO_ACO_Dead,GWO_ACO_Dead_Pct,GWO_ACO_Energy_Pct\n']);

% Write each data row
for i = 1:numel(rows)
    r = rows{i};

    if r{3} == -1
        % Lifetime row — round = -1 is a special marker
        % Only node count, run number, and final round are meaningful here
        fprintf(fid, '%d,%d,LIFETIME,%d\n', r{1}, r{2}, r{4});
    else
        % Regular checkpoint row — write all eight columns
        fprintf(fid, '%d,%d,%d,%.2f,%.2f,%.2f,%.2f,%.2f\n', ...
            r{1}, r{2}, r{3}, r{4}, r{5}, r{6}, r{7}, r{8});
    end
end

fclose(fid);
fprintf('  CSV saved: %s\n', csv_path);
end



% CROSS-N COMPARISON PLOTS
% Overlays averaged curves from all node counts on one figure per metric
% Allows direct visual comparison of how network size affects performance
function plot_allN_comparison_gwo_only(node_counts, nN, NUM_RUNS, RE_runs, Dead_runs, Sur_runs)

    % Pad and average helper — same pattern used throughout
    function avg = pad_avg(cell_row, nr)
        ml = max(cellfun(@numel, cell_row));
        mat = zeros(nr, ml);
        for rr = 1:nr
            v = cell_row{rr};
            if isempty(v)
                continue;
            end
            mat(rr,:) = [v, repmat(v(end), 1, ml-numel(v))];
        end
        avg = mean(mat, 1);
    end

    % --- Figure 1: Dead Node Count across all N ---
    figure('Name','GWO-ACO Only | Dead Node Count Comparison Across All N', ...
        'NumberTitle','off','Position',[100 100 1000 700]);
    hold on;

    for ki = 1:nN
        da = pad_avg(Dead_runs(ki,:), NUM_RUNS);
        % Each N gets its own colored line — MATLAB auto-assigns colors
        plot(1:numel(da), da, 'LineWidth', 1.8, ...
            'DisplayName', sprintf('N=%d', node_counts(ki)));
    end

    xlabel('Round');
    ylabel('Dead Node Count');
    title('Dead Node Count Comparison for N = 50 to 500');
    grid on;
    legend('Location','eastoutside');
    hold off;

    % --- Figure 2: Node Survival Rate across all N ---
    figure('Name','GWO-ACO Only | Node Survival Rate Comparison Across All N', ...
        'NumberTitle','off','Position',[120 120 1000 700]);
    hold on;

    for ki = 1:nN
        sa = pad_avg(Sur_runs(ki,:), NUM_RUNS);
        plot(1:numel(sa), sa, 'LineWidth', 1.8, ...
            'DisplayName', sprintf('N=%d', node_counts(ki)));
    end

    xlabel('Round');
    ylabel('Node Survival Rate (%)');
    title('Node Survival Rate Comparison for N = 50 to 500');
    ylim([0 100]);
    grid on;
    legend('Location','eastoutside');
    hold off;

    % --- Figure 3: Residual Energy across all N ---
    figure('Name','GWO-ACO Only | Residual Energy Comparison Across All N', ...
        'NumberTitle','off','Position',[140 140 1000 700]);
    hold on;

    for ki = 1:nN
        ra = pad_avg(RE_runs(ki,:), NUM_RUNS);
        plot(1:numel(ra), ra, 'LineWidth', 1.8, ...
            'DisplayName', sprintf('N=%d', node_counts(ki)));
    end

    xlabel('Round');
    ylabel('Residual Energy (%)');
    title('Residual Energy Comparison for N = 50 to 500');
    ylim([0 100]);
    grid on;
    legend('Location','eastoutside');
    hold off;
end



% EUCLIDEAN DISTANCE HELPER
% Returns the straight-line distance between two nodes A and B
function d = dist2D(A, B)
d = sqrt((A.x - B.x)^2 + (A.y - B.y)^2);
end



% TRANSMISSION ENERGY
% Calculates energy cost to transmit one packet over distance d
% Uses free-space model (d^2) for short distances
% Uses multipath fading model (d^4) for long distances
% Threshold distance d0 is where the two models are equal
function e = tx_energy(P, d)

% Crossover distance — below this use free-space, above use multipath
d0 = sqrt(P.E_fs / P.E_mp);

if d < d0
    % Free-space model: energy scales with distance squared
    e = P.E_elec * P.k_bits + P.E_fs * P.k_bits * d^2;
else
    % Multipath fading model: energy scales with distance to the fourth power
    e = P.E_elec * P.k_bits + P.E_mp * P.k_bits * d^4;
end
end



% RECEPTION ENERGY
% Calculates energy cost to receive one packet
% Only the circuit energy is paid — no amplifier cost for receiving
function e = rx_energy(P)
e = P.E_elec * P.k_bits;
end