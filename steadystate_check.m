%% Steady-State Test - Single-band HC at Multiple Powers
clear; clc; close all;

%% Common Parameters
pulse_duration = 5; % seconds
npoints = 100000;
dt = pulse_duration/npoints;
shape = 'square';
delta = 10000; % Hz offset
tissue_type = 'hc'; % tissue type
n_samples = 100; % Number of points to plot
nband = '2band';
sample_times = linspace(0, pulse_duration, n_samples);

%% Tissue parameters
tissuepars = init_tissue(tissue_type);
tissuepars.lineshape = 'SL';

%% Test cases 3 power matched set of single- and dual- band data
%B1_powers = [13.8, 6.9, 3.45]; % μT
B1_powers = [19.62, 9.76, 4.88];
colors = {'r', 'b', 'g'};

%% Initialize figure
figure('Position', [100, 100, 800, 600]);
hold on;

fprintf('Single-band HC steady-state tests at %d Hz:\n', delta);
fprintf('=========================================\n');

for idx = 1:length(B1_powers)
    B1_max = B1_powers(idx);
    
    % Generate pulse
    pulse_shape = te_gen_MB_pulse(pulse_duration, npoints, delta, nband, shape);
    b1_band = B1_max * pulse_shape(:);
    
    % Track evolution
    Mz_evolution = zeros(n_samples, 1);
    
    for i = 1:n_samples
        if i == 1
            Mz_evolution(i) = 1; % Initial value (normalized)
        else
            n_current = round(sample_times(i)/dt);
            if n_current > 0
                b1_truncated = b1_band(1:n_current);
                Mz_evolution(i) = te_new_Dualcase_ssSPGR_ihMT_integrate(b1_truncated, B1_max, dt, delta, tissuepars, nband);
            end
        end
    end
    
    % Plot evolution
    plot(sample_times, Mz_evolution, [colors{idx} '-'], 'LineWidth', 2, ...
        'DisplayName', sprintf('%.1f μT', B1_max));
    
    % Plot final value line
    final_value = Mz_evolution(end);
    plot([0, pulse_duration], [final_value, final_value], [colors{idx} '--'], ...
        'LineWidth', 1, 'HandleVisibility', 'off');
    
    % Check steady state
    last_10_percent = round(0.9*n_samples):n_samples;
    variation = std(Mz_evolution(last_10_percent));
    fprintf('B1 = %4.1f μT: Final Mz = %.6f, Variation = %.2e %s\n', ...
        B1_max, final_value, variation, ...
        iif(variation < 1e-5, '✓', '⚠'));
end

%% Format plot
xlabel('Time (s)', 'FontSize', 12);
ylabel('Mz_{free}/M_0', 'FontSize', 12);
title('Single-band HC Steady-State Test at 10 kHz', 'FontSize', 14);
legend('Location', 'best', 'FontSize', 11);
grid on;
ylim([0, 1]);

% Helper function
function result = iif(condition, true_val, false_val)
    if condition
        result = true_val;
    else
        result = false_val;
    end
end