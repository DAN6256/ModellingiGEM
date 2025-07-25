%% PARAMETERS AND INITIAL CONDITIONS

% General simulation parameters
tspan = [0 50];

% Biological parameters
u_max = 1.7;             % Max bacterial growth rate (1/h)
delta = 0.1;             % Bacterial death rate (1/h)
K_B0 = 1e10;             % Initial carrying capacity (cells/mL)
B_0 = 1e7;               % Initial bacterial concentration (cells/mL)

% Urease dynamics
initial_urease = 0;      % Initial urease amount
alpha = 0.05;            % Urease production rate (1/h)
beta = 0.0000000005;     % Urease degradation rate (1/h)
u_alpha = 0.000000001;   % Unused (can be removed)



% Urea hydrolysis
urea_0 = 10000;          % Initial urea concentration (mM)
urea_broken_0 = 0;       % Initial hydrolyzed urea (mM)
K_m = 0.3;               % Michaelis constant (mM)
kcat = 1;               % Turnover rate (1/h)

% Carbonate formation
kpptn = 0.01;            % Precipitation rate constant
Ca = 1.0;                % Calcium concentration (mM)
yCO3 = 1.0;              % Yield of carbonate per urea hydrolyzed


p_max = 10;              % Maximum CaCO3 accumulation

carbonate_0 = 0;
Ca2_0 = Ca ;
CaCO3_0 = 0;


%% SOLVERS AND PLOTTING

U0 = [B_0; initial_urease; urea_0; carbonate_0; Ca2_0; CaCO3_0];
[t, Y] = ode45(@(t, U) odeCarbonateFormation(t, U, K_B0, delta, u_max, u_alpha, beta, K_m, kcat, yCO3, kpptn), tspan, U0);
figure
yyaxis right
plot(t,Y(:,5),LineWidth=2); hold on;
plot(t,Y(:,6), LineWidth=2);hold on;
xlabel("Time(hours");
ylabel("Amount of Ca^{2+} and CaCO_3 (mM?)");
yyaxis left
plot(t,Y(:,1),LineWidth=2);
ylabel("Number of cells(cells/mL)");
title("Plot of Bacteria, Calcium Ions and Calcium Carbonate Dynamics");
legend("Ca^{2+}","CaCO_3","Bacteria");





%% MODEL FUNCTIONS


function dbdt = odefcn_dynamic_KB(t, B, u_max, K_B0, delta, alpha, p_max, t_max)
    P_CaCO3 = p_max * (1.2 * t / t_max);
    K_B = K_B0 * (1 - alpha * (P_CaCO3 / p_max));
    dbdt = u_max * B * (1 - B / K_B) - delta * B;
end


function dUdt = odeCarbonateFormation(t, U, K_B, delta, u_max, alpha, beta, K_m, kcat, yCO3, kpptn)
    B = U(1); E = U(2); S = U(3); C = U(4); Ca = U(5); CaCO3 = U(6) ;

    dUdt = zeros(6,1);
    dUdt(1) = u_max * B * (1 - B / K_B) - delta * B;           % Bacteria
    dUdt(2) = (alpha * B) - (beta * E);                            % Urease
    v = (kcat * E * S) / (K_m + S);                            % Michaelis-Menten rate
    dUdt(3) = -v;                                              % Urea consumed
    dUdt(4) = yCO3 * abs(v) - kpptn * C * Ca;                  % Carbonate formation
    dUdt(5) = - kpptn * C * Ca;                                % Ca2+ formation
    dUdt(6) = kpptn * C * Ca;                                  %CaCO3

end




