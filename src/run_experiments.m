
%% Expermiment 1
clear;
close all;
clc;
parameters;

I_TONIC = 1.0;
wei = 1;
wie = 1;
use_stochastic_GABA = false;

[Ve1, Vi1, T1] = run_network(I_TONIC, wei, wie, use_stochastic_GABA);

figure;
plot(T1, Ve1, 'b', 'LineWidth', 1.5); hold on;
plot(T1, Vi1, 'r', 'LineWidth', 1.5);
title('Experiment 1: Baseline Oscillatory Activity', 'FontSize', 14);
xlabel('Time (ms)', 'FontSize', 12);
ylabel('Membrane Potential (mV)', 'FontSize', 12);
legend({'Excitatory (Ve)', 'Inhibitory (Vi)'}, 'Location', 'Best', 'FontSize', 10);
grid on;
xlim([T1(1), T1(end)]);
set(gca, 'FontSize', 12, 'Box', 'off');

%% Experiment 2: Parametric Sweep over ITonic, wEI, wIE
%% Experiment 2a: Sweep over I_TONIC
clear;
close all;
clc;
parameters;

I_TONIC_vals = linspace(0.02, 4.0, 30);
wei = 1;
wie = 1;
n_vals = length(I_TONIC_vals);
results_I = zeros(n_vals, 6);
results_Ii = zeros(n_vals, 6);
all_Ve = cell(n_vals, 1);
all_Vi = cell(n_vals, 1);

for i = 1:n_vals
    I_TONIC = I_TONIC_vals(i);
    use_stochastic_GABA = false;
    [Ve, Vi, T] = run_network(I_TONIC, wei, wie, use_stochastic_GABA);
    all_Ve{i} = Ve;
    all_Vi{i} = Vi;

    spikes_E = find(Ve(1:end-1) < 0 & Ve(2:end) >= 0);
    freq_E = 0;
    if length(spikes_E) > 1
        freq_E = 1000 / mean(diff(T(spikes_E)));
    end
    if freq_E == 0
        rest_E = mean(Ve(T >= dt_BEFORE & T < dt_BEFORE + dt_TONIC));
    else
        rest_E = NaN;
    end
    results_I(i, :) = [I_TONIC, rest_E, max(Ve), min(Ve), max(Ve)-min(Ve), freq_E];

    spikes_I = find(Vi(1:end-1) < 0 & Vi(2:end) >= 0);
    freq_I = 0;
    if length(spikes_I) > 1
        freq_I = 1000 / mean(diff(T(spikes_I)));
    end
    if freq_I == 0
        rest_I = mean(Vi(T >= dt_BEFORE & T < dt_BEFORE + dt_TONIC));
    else
        rest_I = NaN;
    end
    results_Ii(i, :) = [I_TONIC, rest_I, max(Vi), min(Vi), max(Vi) - min(Vi), freq_I];
end

I_vals = results_I(:, 1);
resting = results_I(:, 2);
Vmax = results_I(:, 3);
Vmin = results_I(:, 4);
amplitude = results_I(:, 5);
frequency = results_I(:, 6);

Vi_rest = results_Ii(:, 2);
Vi_max = results_Ii(:, 3);
Vi_min = results_Ii(:, 4);
Vi_amp = results_Ii(:, 5);
Vi_freq = results_Ii(:, 6);

% Plot 2a.1: Excitatory neuron summary
figure;
plot(I_vals, frequency, 'o-', 'LineWidth', 1.5);
xlabel('I_{TONIC}'); ylabel('Firing Frequency (Hz)');
title('Firing Frequency vs I_{TONIC}'); grid on;

figure;
subplot(2,2,1);
plot(I_vals, frequency, 'o-', 'LineWidth', 1.5);
xlabel('I_{TONIC}'); ylabel('Firing Frequency (Hz)');
title('Firing Frequency vs I_{TONIC}'); grid on;

subplot(2,2,2);
plot(I_vals, amplitude, 'o-', 'LineWidth', 1.5);
xlabel('I_{TONIC}'); ylabel('Spike Amplitude (mV)');
title('Amplitude vs I_{TONIC}'); grid on;

subplot(2,2,3);
plot(I_vals, resting, 'o-', 'LineWidth', 1.5);
xlabel('I_{TONIC}'); ylabel('Resting Potential (mV)');
title('Resting V_m vs I_{TONIC}'); grid on;

subplot(2,2,4);
plot(I_vals, Vmax, 'o-', I_vals, Vmin, 's-', 'LineWidth', 1.5);
xlabel('I_{TONIC}'); ylabel('Membrane Potential (mV)');
legend('Max V_m', 'Min V_m'); title('Max & Min V_m vs I_{TONIC}'); grid on;

sgtitle('Sweep over I_{TONIC} (Excitatory Neuron)');

% Plot 2a.2: Inhibitory neuron summary
figure;
subplot(2,2,1);
plot(I_vals, Vi_freq, 'o-', 'LineWidth', 1.5);
xlabel('I_{TONIC}'); ylabel('Firing Frequency (Hz)');
title('Firing Frequency vs I_{TONIC} (Vi)'); grid on;

subplot(2,2,2);
plot(I_vals, Vi_amp, 'o-', 'LineWidth', 1.5);
xlabel('I_{TONIC}'); ylabel('Spike Amplitude (mV)');
title('Amplitude vs I_{TONIC} (Vi)'); grid on;

subplot(2,2,3);
plot(I_vals, Vi_rest, 'o-', 'LineWidth', 1.5);
xlabel('I_{TONIC}'); ylabel('Resting V_m (mV)');
title('Resting V_m vs I_{TONIC} (Vi)'); grid on;

subplot(2,2,4);
plot(I_vals, Vi_max, 'o-', I_vals, Vi_min, 's-', 'LineWidth', 1.5);
xlabel('I_{TONIC}'); ylabel('Membrane Potential (mV)');
legend('Max V_m', 'Min V_m'); title('Max & Min V_m (Vi)'); grid on;

sgtitle('Sweep over I_{TONIC} (Inhibitory Neuron)');

% Plot 2a.3: Example traces for Ve and Vi at 4 I_TONIC values
example_indices = [3, 8, 12, 18];
figure;
for k = 1:4
    idx = example_indices(k);
    Ve = all_Ve{idx};
    Vi = all_Vi{idx};

    subplot(2, 2, k);
    plot(T, Ve, 'b', 'LineWidth', 1.2); hold on;
    plot(T, Vi, 'r', 'LineWidth', 1.2);
    title(sprintf('I_{TONIC} = %.2f', I_TONIC_vals(idx)));
    xlabel('Time (ms)'); ylabel('Membrane Potential (mV)');
    legend('Ve', 'Vi'); grid on;
    xlim([T(1), T(end)]);
end
sgtitle('Membrane Potentials at Different I_{TONIC} Values');

example_indices = [3, 8, 12, 18];
colors = {'b', 'g', 'r', 'm'};
line_width = 2;
font_size = 12;

figure;
for k = 1:4
    idx = example_indices(k);
    Ve = all_Ve{idx};
    dVe_dt = diff(Ve) / dt;
    Ve_mid = Ve(1:end-1);

    subplot(2, 2, k);
    plot(Ve_mid, dVe_dt, 'Color', colors{k}, 'LineWidth', line_width);
    title(sprintf('Phase Plane: I_{TONIC} = %.2f', I_TONIC_vals(idx)), 'FontSize', font_size);
    xlabel('V_e (mV)', 'FontSize', font_size);
    ylabel('dV_e/dt (mV/ms)', 'FontSize', font_size);
    grid on;
    set(gca, 'FontSize', font_size, 'Box', 'off');
end

sgtitle('Excitatory Neuron Phase Plane Across I_{TONIC}', 'FontSize', 14);


%% Experiment 2b: Sweep over wEI (E to I strength)
clear;
close all;
clc;
parameters;
wEI_vals = linspace(0.2, 2.0, 20);
I_TONIC = 1.0;
wie = 1;

n_vals = length(wEI_vals);
results_wEI_E = zeros(n_vals, 6);
results_wEI_I = zeros(n_vals, 6);
all_Ve_wEI = cell(n_vals, 1);
all_Vi_wEI = cell(n_vals, 1);

for i = 1:n_vals
    wei = wEI_vals(i);
    use_stochastic_GABA = false;
    [Ve, Vi, T] = run_network(I_TONIC, wei, wie, use_stochastic_GABA);
    all_Ve_wEI{i} = Ve;
    all_Vi_wEI{i} = Vi;

    spikes_E = find(Ve(1:end-1) < 0 & Ve(2:end) >= 0);
    freq_E = 0;
    if length(spikes_E) > 1
        freq_E = 1000 / mean(diff(T(spikes_E)));
    end
    if freq_E == 0
        rest_E = mean(Ve(T >= dt_BEFORE & T < dt_BEFORE + dt_TONIC));
    else
        rest_E = NaN;
    end
    results_wEI_E(i,:) = [wei, rest_E, max(Ve), min(Ve), max(Ve) - min(Ve), freq_E];

    spikes_I = find(Vi(1:end-1) < 0 & Vi(2:end) >= 0);
    freq_I = 0;
    if length(spikes_I) > 1
        freq_I = 1000 / mean(diff(T(spikes_I)));
    end
    if freq_I == 0
        rest_I = mean(Vi(T >= dt_BEFORE & T < dt_BEFORE + dt_TONIC));
    else
        rest_I = NaN;
    end
    results_wEI_I(i,:) = [wei, rest_I, max(Vi), min(Vi), max(Vi) - min(Vi), freq_I];
end

wei_vals = results_wEI_E(:, 1);
E_rest = results_wEI_E(:, 2);
E_max = results_wEI_E(:, 3);
E_min = results_wEI_E(:, 4);
E_amp = results_wEI_E(:, 5);
E_freq = results_wEI_E(:, 6);

I_rest = results_wEI_I(:, 2);
I_max = results_wEI_I(:, 3);
I_min = results_wEI_I(:, 4);
I_amp = results_wEI_I(:, 5);
I_freq = results_wEI_I(:, 6);

% Plots 2b.1 for E neuron (wEI sweep)
figure;
subplot(2,2,1);
plot(wei_vals, E_freq, 'o-');
xlabel('w_{EI}');
ylabel('Freq (Hz)');
title('E: Freq vs w_{EI}');
grid on;

subplot(2,2,2);
plot(wei_vals, E_amp, 'o-');
xlabel('w_{EI}');
ylabel('Amplitude (mV)');
title('E: Amp vs w_{EI}');
grid on;

subplot(2,2,3);
plot(wei_vals, E_rest, 'o-');
xlabel('w_{EI}');
ylabel('Rest V_m');
title('E: Rest vs w_{EI}');
grid on;

subplot(2,2,4);
plot(wei_vals, E_max, 'o-', wei_vals, E_min, 's-');
xlabel('w_{EI}');
ylabel('Membrane Potential');
legend('Max', 'Min');
title('E: Max/Min V_m vs w_{EI}');
grid on;

sgtitle('Sweep over w_{EI} (Excitatory Neuron)');

figure;
subplot(2,2,1);
plot(wei_vals, I_freq, 'o-');
xlabel('w_{EI}');
ylabel('Freq (Hz)');
title('I: Freq vs w_{EI}');
grid on;

subplot(2,2,2); plot(wei_vals, I_amp, 'o-');
xlabel('w_{EI}');
ylabel('Amplitude (mV)');
title('I: Amp vs w_{EI}');
grid on;

subplot(2,2,3);
plot(wei_vals, I_rest, 'o-');
xlabel('w_{EI}');
ylabel('Rest V_m');
title('I: Rest vs w_{EI}');
grid on;

subplot(2,2,4);
plot(wei_vals, I_max, 'o-', wei_vals, I_min, 's-');
xlabel('w_{EI}');
ylabel('Membrane Potential');
legend('Max', 'Min');
title('I: Max/Min V_m vs w_{EI}');
grid on;

sgtitle('Sweep over w_{EI} (Inhibitory Neuron)');

% Plots 2b.3 traces for 4 values of wEI
figure;
example_indices = [3, 8, 12, 18];
for k = 1:4
    idx = example_indices(k);
    subplot(2,2,k);
    plot(T, all_Ve_wEI{idx}, 'b', 'LineWidth', 1.2); hold on;
    plot(T, all_Vi_wEI{idx}, 'r', 'LineWidth', 1.2);
    xlabel('Time (ms)'); ylabel('Membrane Potential (mV)');
    title(sprintf('w_{EI} = %.2f', wEI_vals(idx)));
    legend('Ve','Vi'); grid on;
end
sgtitle('Membrane Potentials for Selected w_{EI} Values');

%% Experiment 2c: Sweep over wIE (I to E strength)
clear;
close all;
clc;
parameters;

wIE_vals = linspace(0.2, 2.0, 20);
I_TONIC = 1.0;
wei = 1;

n_vals = length(wIE_vals);
results_wIE_E = zeros(n_vals, 6);
results_wIE_I = zeros(n_vals, 6);
all_Ve_wIE = cell(n_vals, 1);
all_Vi_wIE = cell(n_vals, 1);

for i = 1:n_vals
    wie = wIE_vals(i);
    use_stochastic_GABA = false;
    [Ve, Vi, T] = run_network(I_TONIC, wei, wie, use_stochastic_GABA);
    all_Ve_wIE{i} = Ve;
    all_Vi_wIE{i} = Vi;

    spikes_E = find(Ve(1:end-1) < 0 & Ve(2:end) >= 0);
    freq_E = 0;
    if length(spikes_E) > 1
        freq_E = 1000 / mean(diff(T(spikes_E)));
    end
    if freq_E == 0
        rest_E = mean(Ve(T >= dt_BEFORE & T < dt_BEFORE + dt_TONIC));
    else
        rest_E = NaN;
    end
    results_wIE_E(i,:) = [wie, rest_E, max(Ve), min(Ve), max(Ve) - min(Ve), freq_E];

    spikes_I = find(Vi(1:end-1) < 0 & Vi(2:end) >= 0);
    freq_I = 0;
    if length(spikes_I) > 1
        freq_I = 1000 / mean(diff(T(spikes_I)));
    end
    if freq_I == 0
        rest_I = mean(Vi(T >= dt_BEFORE & T < dt_BEFORE + dt_TONIC));
    else
        rest_I = NaN;
    end
    results_wIE_I(i,:) = [wie, rest_I, max(Vi), min(Vi), max(Vi) - min(Vi), freq_I];
end

wie_vals = results_wIE_E(:, 1);
E_rest = results_wIE_E(:, 2);
E_max = results_wIE_E(:, 3);
E_min = results_wIE_E(:, 4);
E_amp = results_wIE_E(:, 5);
E_freq = results_wIE_E(:, 6);

I_rest = results_wIE_I(:, 2);
I_max = results_wIE_I(:, 3);
I_min = results_wIE_I(:, 4);
I_amp = results_wIE_I(:, 5);
I_freq = results_wIE_I(:, 6);

% Plots 2c.1 for E neuron (wIE sweep)
figure;
plot(wie_vals, E_freq, 'o-');
xlabel('w_{IE}');
ylabel('Freq (Hz)');
title('Sweep over wIE E-neuron: Freq vs w_{IE}');
grid on;

figure;
plot(wie_vals, I_freq, 'o-');
xlabel('w_{IE}');
ylabel('Freq (Hz)');
title('Sweep over wIE I: Freq vs w_{IE}');
grid on;

% Plots 2c.3 traces for 4 values of wIE
figure;
example_indices = [3, 8, 12, 18];
for k = 1:4
    idx = example_indices(k);
    subplot(2,2,k);
    plot(T, all_Ve_wIE{idx}, 'b', 'LineWidth', 1.2); hold on;
    plot(T, all_Vi_wIE{idx}, 'r', 'LineWidth', 1.2);
    xlabel('Time (ms)'); ylabel('Membrane Potential (mV)');
    title(sprintf('w_{IE} = %.2f', wIE_vals(idx)));
    legend('Ve','Vi'); grid on;
end
sgtitle('Membrane Potentials for Selected w_{IE} Values');

%% Experiment 3: Stochastic GABA Condition
clear;
close all;
clc;
parameters;
I_TONIC = 1.0;
wei = 1;
wie = 1;
use_stochastic_GABA = true;

[Ve3, Vi3, T3] = run_network(I_TONIC, wei, wie, use_stochastic_GABA);

figure;
plot(T3, Ve3, 'b', 'LineWidth', 1.5); hold on;
plot(T3, Vi3, 'r', 'LineWidth', 1.5);
title('Experiment 3: Oscillations with Stochastic GABAa Input', 'FontSize', 14);
xlabel('Time (ms)', 'FontSize', 12);
ylabel('Membrane Potential (mV)', 'FontSize', 12);
legend({'Excitatory (Ve)', 'Inhibitory (Vi)'}, 'Location', 'Best', 'FontSize', 10);
grid on;
xlim([T1(1), T1(end)]);
set(gca, 'FontSize', 12, 'Box', 'off');
%% Testing Martyna's code
clear;
close all;
clc;

% Define common model parameters
I_TONIC = 1.0;       % Lower tonic input to reduce excessive inhibition
wEI = 5;          % Stronger E → I
wIE = 3.0;          % Strong I → E
a = 5;              % Very steep sigmoid (strong nonlinearity)
theta = 0.5;         % Threshold for sigmoid
tau_E = 20;         % Faster excitatory dynamics
tau_I = 10;         % Slower inhibitory dynamics
tau_sE = 10;        % Shorter excitatory synaptic filter
tau_sI = 30; 

[rE, rI, T] = run_firing_rates(I_TONIC, wEI, wIE, a, theta, tau_E, tau_I, tau_sE, tau_sI);

% Plot firing rate dynamics
figure;
plot(T, rE, 'b', 'LineWidth', 1.5); hold on;
plot(T, rI, 'r', 'LineWidth', 1.5);
xlabel('Time (ms)', 'FontSize', 12);
ylabel('Firing Rate (a.u.)', 'FontSize', 12);
legend('r_E (Excitatory)', 'r_I (Inhibitory)', 'Location', 'Best', 'FontSize', 10);
title('Experiment 4: Firing Rate Model (Sigmoid Wilson-Cowan)', 'FontSize', 14);
grid on;
set(gca, 'FontSize', 12, 'Box', 'off');
xlim([T(1), T(end)]);