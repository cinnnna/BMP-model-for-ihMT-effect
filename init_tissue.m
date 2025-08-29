%%% pars = init_tissue(name)
%
% Function to generate structure containing tissue properties. Possible
% name inputs are:
%
% - 'random' random initialized value for fitting (optional)
% - 'hc' = 'hair conditioner'
% - 'egg' = simple rounded parameter values
%
% default is set to 'hc'
% The function initialises the following parameters:
%%% Feb 2019 Initialise tissue parameters
function pars = init_tissue(name)
pars = struct;
if nargin==0
    % set to default
    name = 'hc';
end
% Replace all letters with lower case
name = lower(name);
% initialise the empty structure
% Rates are all s^-1, times on seconds, magnetizations referenced to total
% M0=1
%%% Free water pool
pars.free.R1 = [];
pars.free.R2 = [];
%%% semisolid pool
pars.semi.M0 = []; % assume M0f+M0s = 1
pars.semi.R1 = []; %<--- assume same for both s1 and s2
pars.semi.R1D = [];
pars.semi.f = []; % fraction of M0s that is in s1
pars.semi.T2 = []; %<--- assume same for both s1 and s2
%%% overall exchange rate constant
pars.k = [];

switch name
    case 'random'
        rng(42); % or use rng('shuffle') for different results each time
        % Define reasonable ranges for each parameter based on literature
        ranges = struct();
        ranges.free_R1 = [0.1, 2.0];      % 1/T1 for free water (s^-1)
        ranges.free_R2 = [1, 50];         % 1/T2 for free water (s^-1) 
        ranges.semi_M0 = [0.01, 0.3];     % Semisolid fraction
        ranges.semi_R1 = [1, 20];         % Semisolid R1 (s^-1)
        ranges.semi_R1D = [10, 100];      % Dipolar R1 (s^-1)
        ranges.semi_T2 = [5e-6, 50e-6];   % Semisolid T2 (s)
        ranges.k = [1, 20];             % Exchange rate (s^-1)

        % Generate random initial values
        pars.free.R1 = ranges.free_R1(1) + diff(ranges.free_R1) * rand();
        pars.free.R2 = ranges.free_R2(1) + diff(ranges.free_R2) * rand();
        pars.semi.M0 = ranges.semi_M0(1) + diff(ranges.semi_M0) * rand();
        pars.semi.R1 = ranges.semi_R1(1) + diff(ranges.semi_R1) * rand();
        pars.semi.R1D = ranges.semi_R1D(1) + diff(ranges.semi_R1D) * rand();
        pars.semi.f = 0; % Keep this fixed %1 hc 0 eggwhite
        pars.semi.T2 = ranges.semi_T2(1) + diff(ranges.semi_T2) * rand();
        pars.k = ranges.k(1) + diff(ranges.k) * rand();
        pars.lineshape = 'Gaussian';
    
    case 'egg'
        pars.free.R1 = 0.401541;
        pars.free.R2 = 15.779814;
        pars.semi.M0 = 0.198607;
        pars.semi.R1 = 13.386881;
        pars.semi.R1D = 50.002152;
        pars.semi.f = 0;
        pars.semi.T2 = 0.000040; %15.20e-6;
        %%% overall exchange constant
        pars.k = 12.434369;
        pars.lineshape = 'Gaussian';
    
     case 'hc'
        pars.free.R1 = 1.326633; % s⁻¹ (T1f ≈ 1.68 s)  %random value
        pars.free.R2 = 10.661474; % s⁻¹ (T2f ≈ 167 ms)
        % ---- Semisolid (Zeeman + dipolar) pools -----------------------------
        pars.semi.M0 = 0.19658; % total semisolid fraction (1–M0f)
        pars.semi.R1 = 2.89518; % s⁻¹ (Zeeman pool)
        pars.semi.R1D = 88.110912; % s⁻¹ (dipolar-order pool) 300
        pars.semi.f = 1; % fraction that is dipolar-coupled
        pars.semi.T2 = 0.000009; % s (25.9 µs)
        pars.k = 19.834545; % exchange rate
        % ---- Lineshape ------------------------------------------------------
        pars.lineshape = 'SL'; % Super-Lorentzian
        
    otherwise
        error('Unknown tissue type: %s', name);
end

end