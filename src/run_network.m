function [Ve, Vi, T] = run_network(I_TONIC, wei, wie, use_stochastic_GABA)
% Load parameters
parameters;

% Time vector setup
t_max = dt_BEFORE + dt_TONIC + dt_AFTER;
T = 0:dt:t_max;
n_t = length(T);
% Neurotransmitters pulse duration
pulse_dur = 1; % ms
dk_NT = round(pulse_dur/dt); 

% Injected Current
I_inj = zeros(1, n_t);
i1 = round(dt_BEFORE / dt) + 1;
i2 = i1 + round(dt_TONIC / dt);
I_inj(1:i1-1) = I_BEFORE;
I_inj(i1:i2-1) = I_TONIC;
I_inj(i2:end) = I_AFTER;

% Function for alpha and beta of h,n and m of HH neurons
m_Na_alpha = @(Ve, Vt) (-0.32*(Ve - Vt - 13)) / (exp((-(Ve - Vt - 13)) / 4) - 1);
m_Na_beta = @(Ve, Vt) (0.28*(Ve - Vt - 40)) / (exp((Ve - Vt - 40) / 5) - 1);

h_Na_alpha = @(Ve, Vt) 0.128 * exp(-(Ve - Vt -17) / 18);
h_Na_beta = @(Ve, Vt) 4 / (1 + exp(-(Ve - Vt - 40) / 5));

n_K_alpha = @(Ve, Vt) (-0.032 * (Ve - Vt -15)) / (exp(-(Ve - Vt - 15)/5) - 1);
n_K_beta = @(Ve, Vt) 0.5 * exp(-(Ve - Vt -10) / 40);

% Magnesium Blocking for NMDA (Destexhe Paper)
Mg_conc = 1.0; % Magnesium Blocking Concentration
x_NMDA = @(V) 1 ./ (1 + Mg_conc * exp(-0.062 * V) / 3.57);

% Initialize variables for both neurons (E and I) and Neurotransmitters
Ve = zeros(1, n_t); Vi = zeros(1, n_t); Ve(1) = Ve_L; Vi(1) = Vi_L;
I_Le = zeros(1, n_t); I_Li = zeros(1, n_t);
I_Nae = zeros(1, n_t); I_Nai = zeros(1, n_t);
I_Ke = zeros(1, n_t); I_Ki = zeros(1, n_t);
m_Nae = zeros(1, n_t); m_Nai = zeros(1, n_t);
h_Nae = zeros(1, n_t); h_Nai = zeros(1, n_t);
n_Ke = zeros(1, n_t); n_Ki = zeros(1, n_t);
r_AMPA = zeros(1, n_t); I_AMPA = zeros(1, n_t); T_AMPA = zeros(1, n_t);
r_NMDA = zeros(1, n_t); I_NMDA = zeros(1, n_t); T_NMDA = zeros(1, n_t);
r_GABAa = zeros(1, n_t); I_GABAa = zeros(1, n_t); T_GABAa = zeros(1, n_t);
r_GABAb = zeros(1, n_t); I_GABAb = zeros(1, n_t); T_GABAb = zeros(1, n_t);


% Initial conditions
Ve(1) = Ve_L; Vi(1) = Vi_L;

m_Nae_alpha_init = m_Na_alpha(Ve(1), Vte);
m_Nae_beta_init = m_Na_beta(Ve(1), Vte);
m_Nae_steady_state_init = m_Nae_alpha_init / (m_Nae_alpha_init + m_Nae_beta_init);
m_Nae(1) = m_Nae_steady_state_init;

h_Nae_alpha_init = h_Na_alpha(Ve(1), Vte);
h_Nae_beta_init = h_Na_beta(Ve(1), Vte);
h_Nae_steady_state_init = h_Nae_alpha_init / (h_Nae_alpha_init + h_Nae_beta_init);
h_Nae(1) = h_Nae_steady_state_init;

n_Ke_alpha_init = n_K_alpha(Ve(1), Vte);
n_Ke_beta_init = n_K_beta(Ve(1), Vte);
n_Ke_steady_state_init = n_Ke_alpha_init / (n_Ke_alpha_init + n_Ke_beta_init);
n_Ke(1) = n_Ke_steady_state_init;

m_Nai_alpha_init = m_Na_alpha(Vi(1), Vti);
m_Nai_beta_init = m_Na_beta(Vi(1), Vti);
m_Nai_steady_state_init = m_Nai_alpha_init / (m_Nai_alpha_init + m_Nai_beta_init);
m_Nai(1) = m_Nai_steady_state_init;

h_Nai_alpha_init = h_Na_alpha(Vi(1), Vti);
h_Nai_beta_init = h_Na_beta(Vi(1), Vti);
h_Nai_steady_state_init = h_Nai_alpha_init / (h_Nai_alpha_init + h_Nai_beta_init);
h_Nai(1) = h_Nai_steady_state_init;

n_Ki_alpha_init = n_K_alpha(Vi(1), Vti);
n_Ki_beta_init = n_K_beta(Vi(1), Vti);
n_Ki_steady_state_init = n_Ki_alpha_init / (n_Ki_alpha_init + n_Ki_beta_init);
n_Ki(1) = n_Ki_steady_state_init;

% Main simulation loop
for k = 2:n_t
    % ----------------------- Excitatory Neuron --------------------------%
    
    % ---- Leakage Current ---- %

    I_Le(k) = ge_L * (Ve(k-1) - Ve_L);

    % ------- Sodium E Current ------- %

    m_Nae_alpha = m_Na_alpha(Ve(k-1), Vte);
    m_Nae_beta = m_Na_beta(Ve(k-1), Vte);
    m_Nae_steady_state = m_Nae_alpha / (m_Nae_alpha + m_Nae_beta);
    m_Nae_tau = 1 / (m_Nae_alpha + m_Nae_beta);
    m_Nae(k) = m_Nae(k-1) + dt * (m_Nae_steady_state - m_Nae(k-1)) / m_Nae_tau;

    h_Nae_alpha = h_Na_alpha(Ve(k-1), Vte);
    h_Nae_beta = h_Na_beta(Ve(k-1), Vte);
    h_Nae_steady_state =  h_Nae_alpha / (h_Nae_alpha + h_Nae_beta);
    h_Nae_tau = 1 / (h_Nae_alpha + h_Nae_beta);
    h_Nae(k) = h_Nae(k-1) + dt * (h_Nae_steady_state - h_Nae(k-1)) / h_Nae_tau;

    I_Nae(k) = ge_Na * m_Nae(k-1)^3 * h_Nae(k-1) * (Ve(k-1) - Ve_Na);

    % ------- Potassium E Current ------- %

    n_Ke_alpha = n_K_alpha(Ve(k-1), Vte);
    n_Ke_beta = n_K_beta(Ve(k-1), Vte);
    n_Ke_steady_state = n_Ke_alpha / (n_Ke_alpha + n_Ke_beta);
    n_Ke_tau = 1 / (n_Ke_alpha + n_Ke_beta);
    n_Ke(k) = n_Ke(k-1) + dt * (n_Ke_steady_state - n_Ke(k-1)) / n_Ke_tau;

    I_Ke(k) = ge_K * n_Ke(k-1)^4 * (Ve(k-1) - Ve_K);

    % ---- Neurotransmitter Currents ---- %

    dr_GABAa = GABAa_alpha * T_GABAa(k-1) * (1 - r_GABAa(k-1)) - GABAa_beta * r_GABAa(k-1);
    r_GABAa(k) = r_GABAa(k-1) + dt * dr_GABAa;
    I_GABAa(k) = GABAa_g * r_GABAa(k) * (Ve(k-1) - GABAa_Erev);

    dr_GABAb = GABAb_alpha * T_GABAb(k-1) * (1 - r_GABAb(k-1)) - GABAb_beta * r_GABAb(k-1);
    r_GABAb(k) = r_GABAb(k-1) + dt * dr_GABAb;
    I_GABAb(k) = GABAb_g * r_GABAb(k) * (Ve(k-1) - GABAb_Erev);

    % -------------------------------- %
    % ------- TOTAL E CURRENT ------- %
    % -------------------------------- %
    I_total_E = I_inj(k) - I_Le(k) - I_Nae(k) - I_Ke(k) - wie * (I_GABAa(k) + I_GABAb(k));
    Ve(k) = Ve(k-1) + (dt/C) * I_total_E;

    % ----------------------- Inhibitory Neuron --------------------------%

    I_Li(k) = gi_L * (Vi(k-1) - Vi_L);

    % ------- Sodium I Current ------- %

    m_Nai_alpha = m_Na_alpha(Vi(k-1), Vti);
    m_Nai_beta = m_Na_beta(Vi(k-1), Vti);
    m_Nai_steady_state = m_Nai_alpha / (m_Nai_alpha + m_Nai_beta); 
    m_Nai_tau = 1 / (m_Nai_alpha + m_Nai_beta);
    m_Nai(k) = m_Nai(k-1) + dt * (m_Nai_steady_state - m_Nai(k-1)) / m_Nai_tau;


    h_Nai_alpha = h_Na_alpha(Vi(k-1), Vti);
    h_Nai_beta = h_Na_beta(Vi(k-1), Vti);
    h_Nai_steady_state =  h_Nai_alpha / (h_Nai_alpha + h_Nai_beta);
    h_Nai_tau = 1 / (h_Nai_alpha + h_Nai_beta);
    h_Nai(k) = h_Nai(k-1) + dt * (h_Nai_steady_state - h_Nai(k-1)) / h_Nai_tau;

    I_Nai(k) = gi_Na * m_Nai(k-1)^3 * h_Nai(k-1) * (Vi(k-1) - Vi_Na);

    % ------- Potassium I Current ------- %

    n_Ki_alpha = n_K_alpha(Vi(k-1), Vti);
    n_Ki_beta = n_K_beta(Vi(k-1), Vti);
    n_Ki_steady_state = n_Ki_alpha / (n_Ki_alpha + n_Ki_beta);
    n_Ki_tau = 1 / (n_Ki_alpha + n_Ki_beta);
    n_Ki(k) = n_Ki(k-1) + dt * (n_Ki_steady_state - n_Ki(k-1)) / n_Ki_tau;

    I_Ki(k) = gi_K * n_Ki(k-1)^4 * (Vi(k-1) - Vi_K);

    % ---- Neurotransmitter Currents ---- %

    dr_AMPA = AMPA_alpha * T_AMPA(k-1) * (1 - r_AMPA(k-1)) - AMPA_beta * r_AMPA(k-1);
    r_AMPA(k) = r_AMPA(k-1) + dt * dr_AMPA;
    I_AMPA(k) = AMPA_g * r_AMPA(k) * (Ve(k-1) - AMPA_Erev);

    dr_NMDA = NMDA_alpha * T_NMDA(k-1) * (1 - r_NMDA(k-1)) - NMDA_beta * r_NMDA(k-1);
    r_NMDA(k) = r_NMDA(k-1) + dt * dr_NMDA;
    I_NMDA(k) = NMDA_g * r_NMDA(k) * x_NMDA(Ve(k-1)) * (Ve(k-1) - NMDA_Erev);

    % -------------------------------- %
    % ------- TOTAL I CURRENT ------- %
    % -------------------------------- %
    I_total_I = - I_Li(k) - I_Nai(k) - I_Ki(k) - wei * (I_NMDA(k) + I_AMPA(k)) ;
    Vi(k) = Vi(k-1) + (dt/C) * I_total_I;


    % Square transmitter pulse upon E spike 
    if Ve(k-1) < 0 && Ve(k) >= 0
        for n = 0:dk_NT-1
            if (k+n) <= n_t
                T_AMPA(k+n) = T_AMPA(k+n) + A_release;
                T_NMDA(k+n) = T_NMDA(k+n) + A_release;
            end
        end
    end

    % Square transmitter pulse upon I spike 
    if Vi(k-1) < 0 && Vi(k) >= 0
        if use_stochastic_GABA
            % Add stochasticity for n_GABA_spikes
            for s = 1:n_GABA_spikes
                delay_ms = rand() * GABA_delay_range;
                delay_steps = round(delay_ms / dt);
                spike_time = k + delay_steps;

                for n = 0:(dk_NT - 1)
                    idx = spike_time + n; 
                    if spike_time <= n_t
                        T_GABAa(idx) = T_GABAa(idx) + A_release;
                        T_GABAb(idx) = T_GABAb(idx) + A_release;
                    end
                end
            end
        else
            % Deterministic GABAa and GABAb response
            for n = 0:(dk_NT - 1)
                idx = k + n;
                if idx <= n_t
                    T_GABAa(idx) = T_GABAa(idx) + A_release;
                    T_GABAb(idx) = T_GABAb(idx) + A_release;
                end
            end
        end
    end
end

