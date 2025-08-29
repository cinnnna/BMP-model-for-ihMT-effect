function [Mz_free] = ihMT_integrate(b1_band,B1_max,dt,delta,tissuepars,nband)
%%% Modified version using b1_band and delta
%%% Integrates over all entries in b1_band vector

%%% Unpack tissue parameters
M0s = tissuepars.semi.M0;   %<--- TOTAL semisolid fraction
f = tissuepars.semi.f;      %<--- f = fraction of semisolid pool that is dipolar-coupled
M0f = 1-M0s;                %<--- free fraction is 1-M0s
R1 = [tissuepars.free.R1 tissuepars.semi.R1 tissuepars.semi.R1 tissuepars.semi.R1D];
R2f = tissuepars.free.R2;   % semisolid T2, used in lineshape calculation
T2s = tissuepars.semi.T2;
k = tissuepars.k;           % Overall exchange rate for free and both semisolid pools

%% FIRST: State Equilibrium magnetization (no RF)
Mz_0 = [0 0 M0f (1-f)*M0s f*M0s 0 1]';

%% SECOND: Calculate Mz at the end of RF pulse
Mz_with_RF = calculate_steady_state(b1_band, B1_max, dt, delta, tissuepars, nband, Mz_0, M0s, M0f, f, R1, R2f, T2s, k);

%% THIRD: Normalize like experimental data (Mz/M0)
Mz_free = Mz_with_RF(3)/Mz_0(3);

end

function Mz_out = calculate_steady_state(b1_band, B1_max, dt, delta, tissuepars, nband, Mz_0, M0s, M0f, f, R1, R2f, T2s, k);
% gamma for RF calculation
gam = 267.5221; %< rad /s /uT

%%% lineshape - using scalar delta
switch tissuepars.lineshape
    case 'SL'
        [G,w_loc] = SuperLorentzian_lineshape(T2s,delta,'interpzero');% seconds
    case 'Gaussian'
        [G,w_loc] = gauss_lineshape(T2s,delta);% seconds
end

%% Lambda matrix and C are time invariant
% Free pool transverse components - no off resonance here
La = [-R2f 0;0 -R2f];
% the rest
Lb = [-k*M0s-R1(1) k*M0f k*M0f 0;k*(1-f)*M0s -k*M0f-R1(2) 0 0;...
    k*f*M0s 0 -k*M0f-R1(3) 0;0 0 0 -R1(4)];
Lambda = blkdiag(La,Lb);
C = [0 0 R1(1)*M0f R1(2)*(1-f)*M0s R1(3)*f*M0s 0]';

%% Get number of time points from b1_band
nt = size(b1_band, 1);

% Initialize evolution matrix
Xtilde_rf = eye(7);

%% Loop over all entries in b1_band vector
for tt = 1:nt
    % Get current b1 value
    b1_value = b1_band(tt);
    
    % Take the complex pulse components for free pool
    b1x = real(b1_value); % Real part (x-component) of complex pulse
    b1y = imag(b1_value); % Imaginary part (y-component) of complex pulse
    
    OmegaFree = gam*[0 0 -b1y; 0 0 b1x; b1y -b1x 0];
    
    %Semisolid pool
    switch nband
       case '1band'
           w1 = gam*abs(b1_value); 
       case '2band'
           w1 = gam*abs(b1_value)/sqrt(2); 
    end

    %W = pi*gam^2*abs(b1_value)^2*G;
    W = pi*w1^2*G;
    
    % Calculate frequency offset effect
    if abs(delta) == 0
        D = 0;   % W*D = 0 when on-resonance
    else
        D = 2*pi*delta/w_loc;
    end

    %OmegaSemi calculation depends on nband
    switch nband
        case '1band'
            OmegaSemi = [[-W 0 0];[0 -W W*D];[0 W*D -W*D^2]]; 
        case '2band'
            OmegaSemi = [[-2*W 0 0];[0 -2*W 0];[0 0 -2*W*D^2]];
    end
    
    Omega = blkdiag(OmegaFree, OmegaSemi);
    
    % Make overall evolution matrix
    Atilde = cat(1, [(Lambda+Omega) C], zeros(1,7));
    
    Xtilde_rf = expm(Atilde*dt) * Xtilde_rf;
end

% Calculate final state
Mz_out = Xtilde_rf * Mz_0;

end