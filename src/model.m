clear all;
close all;

% Load parameters
parameters;

% Time vector setup
t_max = dt_BEFORE + dt_TONIC + dt_AFTER;
T = 0:dt:t_max;
n_t = length(T);

% Injected Current
I_inj = zeros(1, n_t);
i1 = round(dt_BEFORE / dt) + 1;
i2 = i1 + round(dt_TONIC / dt);
I_inj(1:i1-1) = I_BEFORE;
I_inj(i1:i2-1) = I_TONIC;
I_inj(i2:end) = I_AFTER;

% Initialize variables for both neurons (E and I)
Ve = zeros(1, n_t); Vi = zeros(1, n_t);
I_Le = zeros(1, n_t); I_Li = zeros(1, n_t);
I_Nae = zeros(1, n_t); I_Nai = zeros(1, n_t);
I_Ke = zeros(1, n_t); I_Ki = zeros(1, n_t);
m_Nae = zeros(1, n_t); m_Nai = zeros(1, n_t);
h_Nae = zeros(1, n_t); h_Nai = zeros(1, n_t);
n_Ke = zeros(1, n_t); n_Ki = zeros(1, n_t);
p_AMPA_EI = zeros(1, n_t); I_AMPA_EI = zeros(1, n_t);
p_NMDA_EI = zeros(1, n_t); I_NMDA_EI = zeros(1, n_t);
p_GABAa_IE = zeros(1, n_t); I_GABAa_IE = zeros(1, n_t);
p_GABAb_IE = zeros(1, n_t); I_GABAb_IE = zeros(1, n_t);

E_spikes = zeros(1, n_t);
I_spikes = zeros(1, n_t);

% Initial conditions
Ve(1) = V_L; Vi(1) = V_L;
m_Nae(1) = 1 / (1 + exp(-(Ve(1) + 30)/9.5));
h_Nae(1) = 1/(1 + exp((Ve(1)+ 53)/7));
n_Ke(1) = 1/(1 + exp(-(Ve(1) + 30)/10));
m_Nai(1) = 1/(1 + exp(-(Vi(1) + 30)/9.5));
h_Nai(1) = 1/(1 + exp((Vi(1)+ 53)/7));
n_Ki(1) = 1/(1 + exp(-(Vi(1) + 30)/10));

% Main simulation loop
for k = 2:n_t
    I_Le(k) = g_L * (Ve(k-1) - V_L);

    % Sodium current
    m_Nae(k) = 1 / (1 + exp(-(Ve(k-1) + 30)/9.5)); % Instantaneous activation
    h_Nae_inf = 1 / (1 + exp((Ve(k-1)+ 53)/7));
    tau_h = 0.37 + 2.78 * (1/(1 + exp((Ve(k-1) + 40.5)/6)));
    h_Nae(k) = h_Nae(k-1) + dt * (h_Nae_inf - h_Nae(k-1)) / tau_h;
    I_Nae(k) = g_Na * m_Nae(k-1)^3 * h_Nae(k-1) * (Ve(k-1) - V_Na);

    % Potassium current
    n_Ke_inf = 1/(1 + exp(-(Ve(k-1) + 30)/10));
    tau_n = 0.37 + 1.85 * (1/(1 + exp((Ve(k-1) + 27)/15)));
    n_Ke(k) = n_Ke(k-1) + dt * (n_Ke_inf - n_Ke(k-1)) / tau_n;
    I_Ke(k) = g_K * n_Ke(k-1)^4 * (Ve(k-1) - V_K);

    % Synaptic currents from I neuron
    I_GABAa_IE(k) = g_GABAa * wie * p_GABAa_IE(k-1) * (Ve(k-1) - V_GABAa);
    I_GABAb_IE(k) = g_GABAb * wie * p_GABAb_IE(k-1) * (Ve(k-1) - V_GABAb);

    % Total current for E neuron
    I_total_E = I_inj(k) - I_Le(k) - I_Nae(k) - I_Ke(k) - I_GABAa_IE(k) - I_GABAb_IE(k);
    Ve(k) = Ve(k-1) + (dt/C) * I_total_E;

    % Leak current
    I_Li(k) = g_L * (Vi(k-1) - V_L);

    % Sodium current
    m_Nai(k) = 1 / (1 + exp(-(Vi(k-1) + 30)/9.5)); % Instantaneous activation
    h_Nai_inf = 1 / (1 + exp((Vi(k-1)+ 53)/7));
    tau_h = 0.37 + 2.78 * (1/(1 + exp((Vi(k-1) + 40.5)/6)));
    h_Nai(k) = h_Nai(k-1) + dt * (h_Nai_inf - h_Nai(k-1)) / tau_h;
    I_Nai(k) = g_Na * m_Nai(k-1)^3 * h_Nai(k-1) * (Vi(k-1) - V_Na);

    % Potassium current
    n_Ki_inf = 1/(1 + exp(-(Vi(k-1) + 30)/10));
    tau_n = 0.37 + 1.85 * (1/(1 + exp((Vi(k-1) + 27)/15)));
    n_Ki(k) = n_Ki(k-1) + dt * (n_Ki_inf - n_Ki(k-1)) / tau_n;
    I_Ki(k) = g_K * n_Ki(k-1)^4 * (Vi(k-1) - V_K);

    % Synaptic currents from E neuron
    x_NMDA = 1/(1 + 0.42*exp(-0.062*Vi(k-1)));
    I_AMPA_EI(k) = g_AMPA * wei * p_AMPA_EI(k-1) * (Vi(k-1) - V_AMPA);
    I_NMDA_EI(k) = g_NMDA * wei * p_NMDA_EI(k-1) * x_NMDA * (Vi(k-1) - V_NMDA);

    % Total current for I neuron
    I_total_I = - I_Li(k) - I_Nai(k) - I_Ki(k) - I_AMPA_EI(k) - I_NMDA_EI(k);
    Vi(k) = Vi(k-1) + (dt/C) * I_total_I;

    E_spikes(k) = (Ve(k-1) < 0 && Ve(k) >= 0);
    I_spikes(k) = (Vi(k-1) < 0 && Vi(k) >= 0);

    if E_spikes(k)
        p_AMPA_EI(k) = p_AMPA_EI(k-1) + alpha_AMPA*(1-p_AMPA_EI(k-1));
        p_NMDA_EI(k) = p_NMDA_EI(k-1) + alpha_NMDA*(1-p_NMDA_EI(k-1));
    else
        p_AMPA_EI(k) = p_AMPA_EI(k-1)*exp(-dt*beta_AMPA);
        p_NMDA_EI(k) = p_NMDA_EI(k-1)*exp(-dt*beta_NMDA);
    end

    if I_spikes(k)
        p_GABAa_IE(k) = p_GABAa_IE(k-1) + alpha_GABAa*(1-p_GABAa_IE(k-1));
        p_GABAb_IE(k) = p_GABAb_IE(k-1) + alpha_GABAb*(1-p_GABAb_IE(k-1));
    else
        p_GABAa_IE(k) = p_GABAa_IE(k-1)*exp(-dt*beta_GABAa);
        p_GABAb_IE(k) = p_GABAb_IE(k-1)*exp(-dt*beta_GABAb);
    end
end

% Visualization
figure;
plot(T, Ve, 'b', 'LineWidth', 1.5);
hold on;
plot(T, Vi, 'r', 'LineWidth', 1.5);
xlabel('Time (ms)');
ylabel('Membrane Potential (mV)');
title('E-I Network Dynamics');
legend('Excitatory (E)', 'Inhibitory (I)');
grid on;

% figure;
% plot(T, I_NMDA_EI, 'g'); hold on;  % Excitatory (E→I, negative = depolarizing)
% plot(T, I_AMPA_EI, 'r');          % Inhibitory (I→E, positive = hyperpolarizing)
% legend('NMDA (E->I)','AMPA (E->I)');

% figure;
% plot(T, I_GABAb_IE, 'g'); hold on;  % Excitatory (E→I, negative = depolarizing)
% plot(T, I_GABAa_IE, 'r');          % Inhibitory (I→E, positive = hyperpolarizing)
% legend('GABAb (E->I)','GABAa (E->I)');