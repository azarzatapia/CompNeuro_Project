% Time step

dt = 0.025; % ms

% Current Injection Protocol [ms]

dt_BEFORE = 200;
dt_TONIC = 1000;
dt_AFTER = 200;

I_BEFORE = 0;
I_AFTER = 0;

% Parameters to Play with

I_TONIC = 1.1;
wei = 2.2;
wie = 3;

% Membrane Capacitance

C = 1;

% I Leak

g_L = 0.05;
V_L = -70;

% I Na

g_Na = 24;
V_Na = 50;

% I K

g_K = 3;
V_K = -90;

% Factors important to specific receptors

rise_decay_factor = 0.1;
V_Star = V_L + I_TONIC/g_L;

% I AMPA

V_AMPA = 0;
tau_AMPA = 2.5;
alpha_AMPA = rise_decay_factor/tau_AMPA;
beta_AMPA = 1/tau_AMPA;

% I NMDA

V_NMDA = 0;
tau_NMDA = 10;
alpha_NMDA = rise_decay_factor/tau_NMDA;
beta_NMDA = 1/tau_NMDA;

% I GABAa

V_GABAa = -60;
tau_GABAa = 50;
alpha_GABAa = rise_decay_factor/tau_GABAa;
beta_GABAa = 1/tau_GABAa;

% I GABAb

V_GABAb = -90;
tau_GABAb = 150;
alpha_GABAb = rise_decay_factor/tau_GABAb;
beta_GABAb = 1/tau_GABAb;

% Conductances for receptors

g_Scale = 50;
g_AMPA = g_Scale / abs(V_Star - V_AMPA);
x_NMDA_rest = 1/(1 + 0.42*exp(-0.062*V_Star));
g_NMDA = g_Scale / (x_NMDA_rest * abs(V_Star - V_NMDA));
g_GABAa = g_Scale / abs(V_Star - V_GABAa);
g_GABAb = g_Scale / abs(V_Star - V_GABAb);
