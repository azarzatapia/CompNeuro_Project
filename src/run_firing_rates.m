function [rE, rI, T] = run_firing_rates(I_TONIC, wEI, wIE, a, theta, tau_E, tau_I, tau_sE, tau_sI)
parameters;
t_max = dt_TONIC;
T = 0:dt:t_max;
n_t = length(T);

sigmoid = @(x) 1 ./ (1 + exp(-a * (x - theta)));

% Initialize variables
rE = zeros(1, n_t);
rI = zeros(1, n_t);
sE = zeros(1, n_t); % Filtered excitation
sI = zeros(1, n_t); % Filtered inhibition

% Initial conditions
rE(1) = 0.25;
rI(1) = 0.25;
sE(1) = rE(1);
sI(1) = rI(1);

for k = 2:n_t
    % Synaptic filtering
    dsE = (-sE(k-1) + rE(k-1)) / tau_sE;
    dsI = (-sI(k-1) + rI(k-1)) / tau_sI;
    sE(k) = sE(k-1) + dt * dsE;
    sI(k) = sI(k-1) + dt * dsI;

    % Inputs with synaptic filtering
    E_input = I_TONIC - wIE * sI(k-1);
    I_input = wEI * sE(k-1);

    % Firing rate updates
    drE = (-rE(k-1) + sigmoid(E_input)) / tau_E;
    drI = (-rI(k-1) + sigmoid(I_input)) / tau_I;
    rE(k) = rE(k-1) + dt * drE;
    rI(k) = rI(k-1) + dt * drI;
end
end

