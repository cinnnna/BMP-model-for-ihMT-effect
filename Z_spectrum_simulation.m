% ========================================================================
% Z-spectrum simulation - MODIFIED FOR 2-BAND SYSTEM
% –6 kHz … +6 kHz • 63 sample points
% ========================================================================
%% pulse parameters
pulse_duration = 5; % 5 s saturation pulse
npoints = 100000; % samples in the pulse shape
B1_max = 6.9; % µT (peak B1)
nband = '1band';
shape = 'square'; 
dt = pulse_duration/npoints; %dwell time

switch nband
    case '1band'
        B1_max = B1_max;        
    case '2band'
        B1_max = sqrt(2)*B1_max;
end

%% frequency offsets to sweep
offset_vec = linspace(0, 6e3, 31); % Hz
Mz_vec = zeros(size(offset_vec)); % pre-allocate result

%% tissue parameters
tissuepars = init_tissue('hc');
tissuepars.lineshape = 'SL';

%% sweep over offsets
for k = 1:length(offset_vec)
    delta = abs(offset_vec(k)); % positive magnitude

    % build the unit-max shape for THIS offset
    pulse_shape = gen_pulse(pulse_duration, npoints, delta, nband, shape);
    
    b1_band = B1_max * pulse_shape(:); % N×1, µT

    % propagate once through Bloch-McConnell
    Mz_vec(k) = BMP_integrate(b1_band, B1_max, dt, delta, tissuepars, nband);
end


%% plot the Z-spectrum
% figure('DockControls','on');
plot(offset_vec/1e3, Mz_vec, 'LineWidth', 1.6);
hold on;
% plot(offset_vec/1e3, Mz_vec, 'o'); hold on;
xlabel('\Delta (kHz)');
ylabel('$|M_{z,\mathrm{free}}|$ after RF pulse', 'Interpreter','latex');
title('Simulated Z-spectrum (9.4 T)');
grid on;