% urea hydrolysis
% ureaseDynamics
% Parameters
u_max = 1.7;
K_B0 = 1e10; 
B_0 = 1e7;   

alpha = 0.05;    
p_max = 10;      
delta = 0.1;
tspan = [0 50];
u_alpha = 0.000000001;
beta = 0.0000000005;
initial_urease = 0;
U_0 = [B_0; initial_urease];

% Original model (constant K_B)
[t, B] = ode45(@(t, B) odefcn_constant_KB(t, B, u_max, K_B0, delta), tspan, B_0);

% Dynamic model: K_B changes with time as CaCO3 builds up
[t1, B1] = ode45(@(t, B) odefcn_dynamic_KB(t, B, u_max, K_B0, delta, alpha, p_max, tspan(end)), tspan, B_0);
[t2, B2] = ode45(@(t, U) odeUrease(t, U, K_B0, delta, u_max, u_alpha, beta), tspan, U_0);

% Plotting
figure
plot(t, B, 'r', 'LineWidth', 2); hold on;
plot(t1, B1, 'b--', 'LineWidth', 2);
xlabel("Time (h)")
ylabel("Number of bacteria (cells/mL)")
legend("Constant K_B", "Dynamic K_B")
title("Bacteria Population with Static vs Dynamic Carrying Capacity")
grid on;

figure
plot(t, B, 'r', 'LineWidth', 2)
ylabel("Number of bacteria (cells/mL)")
xlabel("Time (h)")
title("Bacteria Population with Static Carrying Capacity")
grid on;

figure
plot(t1, B1, 'r', 'LineWidth', 2)
ylabel("Number of bacteria (cells/mL)")
title("Bacteria Population with Dynamic Carrying Capacity")
xlabel("Time (h)")
grid on;

figure
yyaxis right
plot(t2, B2(:,2), 'r', 'LineWidth', 2); 
ylabel("Amount of urease")

yyaxis left
plot(t2, B2(:,1), 'b--', 'LineWidth', 2); 
ylabel("Number of bacteria (cells/mL)")

xlabel("Time (h)")
legend("Bacteria","Urease Dynamics")
title("Bacteria vs Urease Dynamics")
grid on;





% Next step: Simulate with dynamic K_B?

% Urea Hydrolysis Simulation
% Parameters
U_max = 0.5;      % Maximum rate (Vmax), will be replaced by kcat * urease
K_m = 0.3;        % Michaelis constant for urease (mM)
urea_0 = 10000;     % Initial urea concentration (mM)

% Combine initial conditions for [Bacteria], [Urease], [Urea]
urea_broken_0 = 0;  % Initially no urea has been broken down
U3_0 = [B_0; initial_urease; urea_0; urea_broken_0];


[t3, Y] = ode45(@(t, U) odeUreaseHydrolysis(t, U, K_B0, delta, u_max, u_alpha, beta, K_m), tspan, U3_0);

% Plot urea degradation
figure
plot(t3, Y(:,3), 'g', 'LineWidth', 2)
xlabel('Time (h)')
ylabel('Urea Concentration (mM)')
title('Urea Hydrolysis via Michaelis-Menten Kinetics; Available vrs Consumed')
grid on
hold on
plot(t3, Y(:,4), 'm', 'LineWidth', 2)
legend('Urea Remaining', 'Urea Broken Down')

figure
yyaxis left
plot(t3, Y(:,1), 'g', 'LineWidth', 2)
xlabel('Time (h)')
ylabel('Amount of bacteria (cell/mL)')
hold on;
yyaxis right
plot(t3, Y(:,3), 'g', 'LineWidth', 2)
hold on;
plot(t3, Y(:,4), 'm', 'LineWidth', 2)
ylabel("Amount of Urea (mM)")
legend('Bacteria Population', "Urea Remaining",'Urea Broken Down')





% Static K_B function
function dbdt = odefcn_constant_KB(t, B, u_max, K_B, delta)
    dbdt = u_max * B * (1 - B / K_B) - delta * B;
end

% Dynamic K_B function based on CaCO3 build-up
function dbdt = odefcn_dynamic_KB(t, B, u_max, K_B0, delta, alpha, p_max, t_max)
    P_CaCO3 = p_max * (1.2*t / t_max);         
    K_B = K_B0 * (1 - alpha * (P_CaCO3 / p_max));  
    dbdt = u_max * B * (1 - B / K_B) - delta * B;
end

function dudt = odeUrease(t, U, K_B, delta, u_max, alpha, beta)
  B = U(1);
  urease = U(2); 

  dudt = zeros(2,1);
  dudt(1) = u_max * B * (1 - B / K_B) - delta * B;
  dudt(2) = (alpha * B) - (beta * urease);
end

function dUdt = odeUreaseHydrolysis(t, U, K_B, delta, u_max, alpha, beta, K_m)
    B = U(1);         % Bacteria
    urease = U(2);    % Urease
    urea = U(3);      % Urea
    urea_broken = U(4);  % Urea broken down

    kcat = 1.0;  % Turnover rate (s⁻¹) – adjust to your system

    dUdt = zeros(4,1);
    
    % Bacteria growth
    dUdt(1) = u_max * B * (1 - B / K_B) - delta * B;

    % Urease production and decay
    dUdt(2) = (alpha * B) - (beta * urease);

    % Urea hydrolysis (Michaelis-Menten)
    v_urea = (kcat * urease * urea) / (K_m + urea);

    dUdt(3) = -v_urea;              % Urea consumed
    dUdt(4) = v_urea;               % Urea broken down

end