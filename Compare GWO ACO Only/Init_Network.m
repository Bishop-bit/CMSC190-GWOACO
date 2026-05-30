function [Nodes, sink, Params] = Init_Network(n, seed)
% Init_Network the .m file that initializes the WSN field 
% (Used by Compare_GWO_ACO_Only and the Compare_All_Protocols)
% Inputs:  n    = number of sensor nodes
%          seed = random seed for reproducibility
% Outputs: Nodes  = array of sensor node structs
%          sink   = base station position
%          Params = simulation parameter struct

% If no seed is provided, use 'shuffle' (clock-based random seed)
% Ensures different node layouts each run when seed is not specified
if nargin < 2, seed = 'shuffle'; end

% Seed MATLAB's random number generator for reproducible node placement
% rng(): all randomness flows from this
rng(seed);


% SIMULATION AREA

% Define the sensing field as a 100m x 100m square area
Params.xm = 100;  % Field width in meters
Params.ym = 100;  % Field height in meters

% Place the base station (sink) at the center of the field
sink.x = Params.xm / 2;  % sink x-coordinate = 50m
sink.y = Params.ym / 2;  % sink y-coordinate = 50m

% Initial energy of each sensor node = 0.5 Joules (full battery)
Params.E0 = 0.5;           % Joules

% Energy dissipated by radio electronics to transmit or receive 1 bit
% Paid by BOTH sender and receiver regardless of distance
Params.E_elec = 50e-9;     % Joules/bit (50 nanojoules)

% Amplifier energy for free-space model (short-range transmission, d < d0)
% Energy cost scales with d^2: used when signal path is clear
Params.E_fs = 10e-12;      % Joules/bit/m^2 (10 picojoules)

% Amplifier energy for multipath fading model (long-range, d >= d0)
% Energy cost scales with d^4 — used when signal suffers reflections
Params.E_mp = 0.0013e-12;  % Joules/bit/m^4 (0.0013 picojoules)


% DATA AND AGGREGATION PARAMETERS

% Packet size in bits: each sensor transmits 4000 bits per round
Params.k_bits = 4000;      % bits per packet

% Energy cost for data aggregation at the cluster head
% CH combines data from all its members before forwarding to sink
Params.E_agg = 5e-9;       % Joules/bit (5 nanojoules)


% SENSING AND SIMULATION CONTROL

% Sensing radius of each node: how far it can detect events
Params.R_sense = 10;       % meters

% Maximum number of simulation rounds before forced termination (this is
% for space and time complexity efficiency - one of the way how to run this
% code faster)
Params.rmax = 5000;        % rounds

% Alias for E_fs: some internal functions refer to amplifier energy as E_amp
Params.E_amp = Params.E_fs;

% Round interval for saving network snapshots during simulation
% Every 50 rounds a state snapshot can be logged for analysis
Params.snapshot_round = 50; % rounds


% NODE INITIALIZATION

% Initialize Nodes as an empty struct array before the loop
Nodes = struct([]);

for i = 1:n
    % Unique identifier for each node (1 to n)
    Nodes(i).id = i;

    % Random x and y position within the 100x100m field
    % rand returns a value in [0,1], scaled to field dimensions
    Nodes(i).x = rand * Params.xm;   % x-coordinate in meters
    Nodes(i).y = rand * Params.ym;   % y-coordinate in meters

    % Current energy: starts at full battery (E0 = 0.5J)
    % Decremented every round based on transmission and reception costs
    Nodes(i).E = Params.E0;

    % Alive status: 1 = active, 0 = dead (energy depleted)
    % Node is removed from all operations once alive = 0
    Nodes(i).alive = 1;

    % Cluster head assignment: -1 means not yet assigned to any cluster
    % Updated each round during GWO cluster formation phase
    Nodes(i).cluster_head = -1;

    % Node density weight: used in some routing priority calculations
    % Initialized to 0; may be updated during simulation
    Nodes(i).rho = 0;

    % Last round this node served as a cluster head
    % Set to -inf so all nodes are fully eligible as CH in round 1
    % Using -inf instead of 0 avoids false cooldown bias at start
    Nodes(i).lastCH = -inf;

    % Data buffer: amount of data waiting to be transmitted
    % Starts empty; accumulates sensing data each round
    Nodes(i).buffer = 0;

    % Relay load counter: tracks how many times this node forwarded
    % data on behalf of other cluster heads in the current round
    % Used by ACO path cost function to penalize overloaded relays
    Nodes(i).relay_load = 0;
end
end