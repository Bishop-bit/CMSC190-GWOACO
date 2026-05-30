function [Nodes, sink, Params] = Init_Network(n, seed)
if nargin < 2, seed = 'shuffle'; end
rng(seed);
Params.xm = 100;
Params.ym = 100;
sink.x = Params.xm / 2;
sink.y = Params.ym / 2;
Params.E0      = 0.5;
Params.E_elec  = 50e-9;
Params.E_fs    = 10e-12;
Params.E_mp    = 0.0013e-12;
Params.k_bits  = 4000;
Params.E_agg   = 5e-9;
Params.R_sense = 15;
Params.rmax    = 3000;
Params.E_amp   = Params.E_fs;
Params.snapshot_round = 50;
Nodes = struct([]);
for i = 1:n
    Nodes(i).id           = i;
    Nodes(i).x            = rand * Params.xm;
    Nodes(i).y            = rand * Params.ym;
    Nodes(i).E            = Params.E0;
    Nodes(i).alive        = 1;
    Nodes(i).cluster_head = -1;
    Nodes(i).rho          = 0;
    Nodes(i).lastCH       = -inf;
    Nodes(i).buffer       = 0;
    Nodes(i).relay_load   = 0;
end
end
