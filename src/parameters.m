% Time step

dt = 0.0025; % ms

% Current Injection Protocol [ms]

dt_BEFORE = 100;
dt_TONIC = 1000;
dt_AFTER = 300;

I_AFTER = 0;
I_BEFORE = 0;

% Membrane Capacitance

C = 1;

% ---------------- EXCITATORY NEURON PARAMETERS (Pyramidal Neuron) ---------------- %

% E Threshold Potential
Vte = -56.2;

% E Leakage
ge_L = 0.05; % mS.cm-2
Ve_L = -70; % mV

% E Na --------

ge_Na = 24; % mS.cm-2
Ve_Na = 50; % mV

% E K --------

ge_K = 3; % mS.cm-2
Ve_K = -90; % mV

% ---------------- INHIBITORY NEURON PARAMETERS (Interneuron) ---------------- %

% IThreshold Potential
Vti = -67.9;

gi_L = 0.05; % mS.cm-2
Vi_L = -70; % mV

% I Na --------

gi_Na = 24; % mS.cm-2
Vi_Na = 50; % mV

% I K --------

gi_K = 3; % mS.cm-2
Vi_K = -90; % mV

% ---------------- NEUROTRANSMITTERS ---------------- %
% times faster than closing
rise_decay_factor = 0.1;

% Transmitter release (if you want to customize later)
tau_release = 0.3;   % ms
A_release = 1.0; % Concentration of Neurotransmitter

% AMPA
AMPA_tau    = 2.5;
AMPA_alpha  = rise_decay_factor / AMPA_tau;
AMPA_beta   = 1 / AMPA_tau;
AMPA_Erev   = 0;
AMPA_g      = 5;   % Total for ~1000 synapses

% NMDA
NMDA_tau    = 50;
NMDA_alpha  = rise_decay_factor / NMDA_tau;
NMDA_beta   = 1 / NMDA_tau;
NMDA_Erev   = 0;
NMDA_g      = 5;

% GABAa
GABAa_tau   = 10;
GABAa_alpha = rise_decay_factor / GABAa_tau;
GABAa_beta  = 1 / GABAa_tau;
GABAa_Erev  = -60;
GABAa_g     = 30;

% GABAb 
GABAb_tau   = 150;
GABAb_alpha = rise_decay_factor / GABAb_tau;
GABAb_beta  = 1 / GABAb_tau;
GABAb_Erev  = -90;
GABAb_g     = 5;

% --- Stochastic GABA Parameters ---
n_GABA_spikes = 5; 
GABA_delay_range = 200;% ms
