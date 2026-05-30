function Compare_All_Protocols()

clear; close all; clc;
rng('shuffle');

node_counts = 50:50:500;
nN          = numel(node_counts);
NUM_RUNS    = 5;
CHKPT       = 100;

E0          = 0.5;
RMAX_FIXED  = 5000;

out_dir = 'simulation_output_all_protocols';
if ~exist(out_dir, 'dir')
    mkdir(out_dir);
end

ALG_NAMES = {'GWO-ACO','LEACH','KEAC','RLBEEP','EEHCHR','DC-GWO'};
ALG_COLORS = [
    0.0000 0.4470 0.7410; % GWO-ACO (blue)
    0.8500 0.3250 0.0980; % LEACH (orange)
    0.9290 0.6940 0.1250; % KEAC (yellow)
    0.4940 0.1840 0.5560; % RLBEEP (purple)
    0.4660 0.6740 0.1880; % EEHCHR (green)
    0.3010 0.7450 0.9330; % DC-GWO (cyan)
];
ALG_FUNS  = {@ALG_GWOACO,@ALG_LEACH,@ALG_KEAC,@ALG_RLBEEP,@ALG_EEHCHR,@ALG_DCGWO};
nAlg      = numel(ALG_NAMES);

all_rows = {};

LT_all   = cell(1, nAlg);
RE_runs  = cell(1, nAlg);
Dead_runs= cell(1, nAlg);
Sur_runs = cell(1, nAlg);

for a = 1:nAlg
    LT_all{a}    = zeros(nN, NUM_RUNS);
    RE_runs{a}   = cell(nN, NUM_RUNS);
    Dead_runs{a} = cell(nN, NUM_RUNS);
    Sur_runs{a}  = cell(nN, NUM_RUNS);
end

for ki = 1:nN
    n = node_counts(ki);
    fprintf('\n=== N=%d ===\n', n);
    fprintf('  E0=%.4fJ  rmax=%d\n', E0, RMAX_FIXED);

    nc_rows = {};

    for run = 1:NUM_RUNS
        seed = randi(1e6);
        fprintf('  Run %d/%d (seed=%d)\n', run, NUM_RUNS, seed);

        [Nodes0, sink, P] = Init_Network(n, seed);
        P.rmax = RMAX_FIXED;
        P.E0   = E0;

        for i = 1:n
            Nodes0(i).E            = E0;
            Nodes0(i).alive        = 1;
            Nodes0(i).cluster_head = -1;
            Nodes0(i).relay_load   = 0;
        end

        outs = cell(1, nAlg);
        for a = 1:nAlg
            NodesA = Nodes0;
            outs{a} = ALG_FUNS{a}(NodesA, sink, P);

            LT_all{a}(ki,run)    = outs{a}.fin;
            RE_runs{a}{ki,run}   = outs{a}.re;
            Dead_runs{a}{ki,run} = outs{a}.dead;
            Sur_runs{a}{ki,run}  = outs{a}.sur;

            fprintf('    %-8s LT=%d\n', ALG_NAMES{a}, outs{a}.fin);
        end

        % ---------------------------------------------------------
        % CSV rows at checkpoint rounds
        % ---------------------------------------------------------
        max_fin = 0;
        for a = 1:nAlg
            max_fin = max(max_fin, outs{a}.fin);
        end

        chk_rounds = unique([CHKPT:CHKPT:max_fin, max_fin]);

        for rr = chk_rounds
            row = cell(1, 2 + 1 + 5*nAlg);
            row{1} = n;
            row{2} = run;
            row{3} = rr;

            col = 4;
            for a = 1:nAlg
                outa = outs{a};

                if rr <= numel(outa.sur)
                    alive_pct  = outa.sur(rr);
                    dead_count = outa.dead(rr);
                    dead_pct   = 100 - alive_pct;
                    energy_pct = outa.re(rr);
                    alive_cnt  = n - dead_count;
                else
                    alive_pct  = outa.sur(end);
                    dead_count = outa.dead(end);
                    dead_pct   = 100 - alive_pct;
                    energy_pct = outa.re(end);
                    alive_cnt  = n - dead_count;
                end

                row{col}   = alive_cnt;   col = col + 1;
                row{col}   = alive_pct;   col = col + 1;
                row{col}   = dead_count;  col = col + 1;
                row{col}   = dead_pct;    col = col + 1;
                row{col}   = energy_pct;  col = col + 1;
            end

            nc_rows{end+1}  = row; %#ok<AGROW>
            all_rows{end+1} = row; %#ok<AGROW>
        end

        for a = 1:nAlg
            lt_row = {n, run, ['LIFETIME_' ALG_NAMES{a}], outs{a}.fin};
            nc_rows{end+1}  = lt_row; %#ok<AGROW>
            all_rows{end+1} = lt_row; %#ok<AGROW>
        end

        figure('Name', sprintf('All Protocols | N=%d Run %d/%d', n, run, NUM_RUNS), ...
            'NumberTitle','off','Position',[50 50 1200 850]);

        subplot(3,1,1); hold on;
        for a = 1:nAlg
            sur = outs{a}.sur;
            fin = outs{a}.fin;
            chk = CHKPT:CHKPT:numel(sur);

            plot(1:numel(sur), sur, 'LineWidth', 1.5, 'Color', ALG_COLORS(a,:), 'DisplayName', ALG_NAMES{a});
            if ~isempty(chk)
                plot(chk, sur(chk), 's', 'Color', ALG_COLORS(a,:), 'MarkerSize', 6, 'HandleVisibility','off');
            end
            plot([fin fin], [0 100], ':', 'LineWidth', 1.0, 'HandleVisibility','off');
        end
        xlabel('Round'); ylabel('Alive (%)');
        title('Node Survival Rate');
        ylim([0 100]); grid on; legend('Location','eastoutside');

        subplot(3,1,2); hold on;
        for a = 1:nAlg
            dead = outs{a}.dead;
            fin = outs{a}.fin;
            chk = CHKPT:CHKPT:numel(dead);

            plot(1:numel(dead), dead, 'LineWidth', 1.5, 'Color', ALG_COLORS(a,:), 'DisplayName', ALG_NAMES{a});
            if ~isempty(chk)
                plot(chk, dead(chk), 'o', 'MarkerSize', 6, 'HandleVisibility','off');
            end
            plot([fin fin], [0 n], ':', 'LineWidth', 1.0, 'HandleVisibility','off');
        end
        xlabel('Round'); ylabel('Dead Nodes');
        title('Dead Node Count');
        ylim([0 n]); grid on; legend('Location','eastoutside');

        subplot(3,1,3); hold on;
        for a = 1:nAlg
            re = outs{a}.re;
            fin = outs{a}.fin;
            chk = CHKPT:CHKPT:numel(re);

            plot(1:numel(re), re, 'LineWidth', 1.5, 'Color', ALG_COLORS(a,:), 'DisplayName', ALG_NAMES{a});
            if ~isempty(chk)
                plot(chk, re(chk), 'd', 'MarkerSize', 6, 'HandleVisibility','off');
            end
            plot([fin fin], [0 100], ':', 'LineWidth', 1.0, 'HandleVisibility','off');
        end
        xlabel('Round'); ylabel('Residual Energy (%)');
        title('Residual Energy');
        ylim([0 100]); grid on; legend('Location','eastoutside');

        sgtitle(sprintf('Protocol Comparison | N=%d | Run %d/%d | Seed=%d', ...
            n, run, NUM_RUNS, seed), 'FontWeight','bold');
        drawnow;
    end

    csv_path = fullfile(out_dir, sprintf('N%d_detailed.csv', n));
    write_csv_all_protocols(nc_rows, csv_path, ALG_NAMES);

    plot_averaged_all_protocols(ki, n, NUM_RUNS, CHKPT, ALG_NAMES, LT_all, RE_runs, Dead_runs, Sur_runs);
end

plot_global_summary_all_protocols(node_counts, nN, NUM_RUNS, ALG_NAMES, LT_all, RE_runs, Dead_runs, Sur_runs);
plot_allN_comparison_all_protocols(node_counts, nN, NUM_RUNS, ALG_NAMES, RE_runs, Dead_runs, Sur_runs);

write_csv_all_protocols(all_rows, fullfile(out_dir, 'ALL_NODES_master.csv'), ALG_NAMES);

fprintf('\nDone. CSVs in: %s\n', out_dir);
end

function plot_averaged_all_protocols(ki, n, NUM_RUNS, CHKPT, ALG_NAMES, LT_all, RE_runs, Dead_runs, Sur_runs)
nAlg = numel(ALG_NAMES);

figure('Name', sprintf('All Protocols | N=%d Averaged (%d runs)', n, NUM_RUNS), ...
    'NumberTitle','off','Position',[120 120 1200 860]);

subplot(3,1,1); hold on;
for a = 1:nAlg
    sa_avg = pad_avg(Sur_runs{a}(ki,:), NUM_RUNS);
    lt_avg = mean(LT_all{a}(ki,:));
    chk = CHKPT:CHKPT:numel(sa_avg);

    plot(1:numel(sa_avg), sa_avg, 'LineWidth', 2.0, 'DisplayName', ALG_NAMES{a});
    if ~isempty(chk)
        plot(chk, sa_avg(chk), 's', 'MarkerSize', 6, 'HandleVisibility','off');
    end
    plot([lt_avg lt_avg], [0 100], ':', 'LineWidth', 1.0, 'HandleVisibility','off');
end
xlabel('Round'); ylabel('Alive (%)');
title('Node Survival Rate');
ylim([0 100]); grid on; legend('Location','eastoutside');

subplot(3,1,2); hold on;
for a = 1:nAlg
    da_avg = pad_avg(Dead_runs{a}(ki,:), NUM_RUNS);
    lt_avg = mean(LT_all{a}(ki,:));
    chk = CHKPT:CHKPT:numel(da_avg);

    plot(1:numel(da_avg), da_avg, 'LineWidth', 2.0, 'DisplayName', ALG_NAMES{a});
    if ~isempty(chk)
        plot(chk, da_avg(chk), 'o', 'MarkerSize', 6, 'HandleVisibility','off');
    end
    plot([lt_avg lt_avg], [0 n], ':', 'LineWidth', 1.0, 'HandleVisibility','off');
end
xlabel('Round'); ylabel('Dead Nodes');
title('Dead Node Count');
ylim([0 n]); grid on; legend('Location','eastoutside');

subplot(3,1,3); hold on;
for a = 1:nAlg
    ra_avg = pad_avg(RE_runs{a}(ki,:), NUM_RUNS);
    lt_avg = mean(LT_all{a}(ki,:));
    chk = CHKPT:CHKPT:numel(ra_avg);

    plot(1:numel(ra_avg), ra_avg, 'LineWidth', 2.0, 'DisplayName', ALG_NAMES{a});
    if ~isempty(chk)
        plot(chk, ra_avg(chk), 'd', 'MarkerSize', 6, 'HandleVisibility','off');
    end
    plot([lt_avg lt_avg], [0 100], ':', 'LineWidth', 1.0, 'HandleVisibility','off');
end
xlabel('Round'); ylabel('Residual Energy (%)');
title('Residual Energy');
ylim([0 100]); grid on; legend('Location','eastoutside');

sgtitle(sprintf('All Protocols | N=%d | %d runs averaged', n, NUM_RUNS), ...
    'FontWeight','bold','FontSize',12);
end

function plot_global_summary_all_protocols(node_counts, nN, NUM_RUNS, ALG_NAMES, LT_all, RE_runs, Dead_runs, Sur_runs)
nAlg = numel(ALG_NAMES);

figure('Name','All Protocols | Network Lifetime vs N','NumberTitle','off');
hold on;
for a = 1:nAlg
    LT_mean = mean(LT_all{a}, 2)';
    plot(node_counts, LT_mean, 'o-', 'LineWidth', 2, 'DisplayName', ALG_NAMES{a});
end
xlabel('Number of Nodes (N)');
ylabel('Network Lifetime (Rounds)');
title(sprintf('Average Lifetime vs N (%d runs)', NUM_RUNS));
grid on;
legend('Location','eastoutside');
hold off;

nCols = min(nN,5);
nRows = ceil(nN/nCols);

for a = 1:nAlg
    fig_re = figure('Name',sprintf('%s | Residual Energy All N', ALG_NAMES{a}), ...
        'NumberTitle','off','Position',[50 50 300*nCols 300*nRows]);
    fig_dead = figure('Name',sprintf('%s | Dead Nodes All N', ALG_NAMES{a}), ...
        'NumberTitle','off','Position',[50 50 300*nCols 300*nRows]);
    fig_sur = figure('Name',sprintf('%s | Survival Rate All N', ALG_NAMES{a}), ...
        'NumberTitle','off','Position',[50 50 300*nCols 300*nRows]);

    for ki = 1:nN
        nc = node_counts(ki);
        ra = pad_avg(RE_runs{a}(ki,:), NUM_RUNS);
        da = pad_avg(Dead_runs{a}(ki,:), NUM_RUNS);
        sa = pad_avg(Sur_runs{a}(ki,:), NUM_RUNS);

        figure(fig_re); subplot(nRows,nCols,ki);
        plot(1:numel(ra), ra, 'LineWidth', 1.5);
        xlabel('Round'); ylabel('Energy (%)');
        title(sprintf('N=%d', nc)); grid on; ylim([0 100]);

        figure(fig_dead); subplot(nRows,nCols,ki);
        plot(1:numel(da), da, 'LineWidth', 1.5);
        xlabel('Round'); ylabel('Dead');
        title(sprintf('N=%d', nc)); grid on; ylim([0 nc]);

        figure(fig_sur); subplot(nRows,nCols,ki);
        plot(1:numel(sa), sa, 'LineWidth', 1.5);
        xlabel('Round'); ylabel('Alive (%)');
        title(sprintf('N=%d', nc)); grid on; ylim([0 100]);
    end

    figure(fig_re);   sgtitle(sprintf('%s | Residual Energy (avg %d runs)', ALG_NAMES{a}, NUM_RUNS));
    figure(fig_dead); sgtitle(sprintf('%s | Dead Node Count (avg %d runs)', ALG_NAMES{a}, NUM_RUNS));
    figure(fig_sur);  sgtitle(sprintf('%s | Node Survival Rate (avg %d runs)', ALG_NAMES{a}, NUM_RUNS));
end
end

function plot_allN_comparison_all_protocols(node_counts, nN, NUM_RUNS, ALG_NAMES, RE_runs, Dead_runs, Sur_runs)
nAlg = numel(ALG_NAMES);

for a = 1:nAlg
    figure('Name',sprintf('%s | Dead Node Count Comparison Across All N', ALG_NAMES{a}), ...
        'NumberTitle','off','Position',[100 100 1000 700]);
    hold on;
    for ki = 1:nN
        da = pad_avg(Dead_runs{a}(ki,:), NUM_RUNS);
        plot(1:numel(da), da, 'LineWidth', 1.8, 'DisplayName', sprintf('N=%d', node_counts(ki)));
    end
    xlabel('Round'); ylabel('Dead Node Count');
    title(sprintf('%s | Dead Node Count Comparison for N = 50 to 500', ALG_NAMES{a}));
    grid on; legend('Location','eastoutside'); hold off;

    figure('Name',sprintf('%s | Node Survival Rate Comparison Across All N', ALG_NAMES{a}), ...
        'NumberTitle','off','Position',[120 120 1000 700]);
    hold on;
    for ki = 1:nN
        sa = pad_avg(Sur_runs{a}(ki,:), NUM_RUNS);
        plot(1:numel(sa), sa, 'LineWidth', 1.8, 'DisplayName', sprintf('N=%d', node_counts(ki)));
    end
    xlabel('Round'); ylabel('Node Survival Rate (%)');
    title(sprintf('%s | Node Survival Rate Comparison for N = 50 to 500', ALG_NAMES{a}));
    ylim([0 100]); grid on; legend('Location','eastoutside'); hold off;

    figure('Name',sprintf('%s | Residual Energy Comparison Across All N', ALG_NAMES{a}), ...
        'NumberTitle','off','Position',[140 140 1000 700]);
    hold on;
    for ki = 1:nN
        ra = pad_avg(RE_runs{a}(ki,:), NUM_RUNS);
        plot(1:numel(ra), ra, 'LineWidth', 1.8, 'DisplayName', sprintf('N=%d', node_counts(ki)));
    end
    xlabel('Round'); ylabel('Residual Energy (%)');
    title(sprintf('%s | Residual Energy Comparison for N = 50 to 500', ALG_NAMES{a}));
    ylim([0 100]); grid on; legend('Location','eastoutside'); hold off;
end
end

% ======================================================================
% CSV writer
% ======================================================================
function write_csv_all_protocols(rows, csv_path, ALG_NAMES)
fid = fopen(csv_path, 'w');

fprintf(fid, 'Nodes,Run,Round');
for a = 1:numel(ALG_NAMES)
    nm = matlab.lang.makeValidName(strrep(ALG_NAMES{a},'-','_'));
    fprintf(fid, ',%s_Alive,%s_Alive_Pct,%s_Dead,%s_Dead_Pct,%s_Energy_Pct', ...
        nm, nm, nm, nm, nm);
end
fprintf(fid, '\n');

for i = 1:numel(rows)
    r = rows{i};

    if ischar(r{3}) || isstring(r{3})
        fprintf(fid, '%d,%d,%s,%d\n', r{1}, r{2}, string(r{3}), r{4});
        continue;
    end

    fprintf(fid, '%d,%d,%d', r{1}, r{2}, r{3});
    for j = 4:numel(r)
        fprintf(fid, ',%.6f', r{j});
    end
    fprintf(fid, '\n');
end

fclose(fid);
fprintf('  CSV: %s\n', csv_path);
end

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