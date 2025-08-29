% MICP pH dynamics model simulation (closed system only)

% ---- Parameters ----
P.Km       = 5e-3;        % M (urease Km)
P.Vmax     = 0.05;        % M/h per OD (urea hydrolysis max rate)
P.Ka_NH4   = 10^(-9.25);  % NH4+ <-> NH3 + H+
P.Ka1      = 10^(-6.3);   % H2CO3 <-> HCO3- + H+
P.Ka2      = 10^(-10.33); % HCO3- <-> CO3^2- + H+
P.Ksp      = 2.8e-9;      % CaCO3 solubility (mol^2/L^2)
P.k_precip = 0.5;         % L/(mol·h) (precipitation kinetics)
P.mu       = 1.2;         % 1/h (max specific growth rate)[May change this]
P.B_max    = 5.0;         % OD (carrying capacity for biomass)

% ---- Initial conditions ----
U0    = 0.30;         % M (initial urea)
Ca0   = 0.30;         % M (initial Ca2+)
P.Cl0 = 2*Ca0;        % M (Cl- from CaCl2, charge balance)
N0    = 0.0;          % M (initial NH3 + NH4+)
C0    = 1e-8;         % M (initial total inorganic carbon, small)
B0    = 0.005;        % OD (initial biomass concentration)

y0    = [U0 N0 C0 Ca0 B0];
tspan = [0 50];       % hours
options = odeset('NonNegative',1:5);

% ---- Solve closed-system scenarios ----
sol1 = ode45(@(t,y) odefun(t,y,false,P), tspan, y0, options); % closed, constant B
sol2 = ode45(@(t,y) odefun(t,y,true, P), tspan, y0, options);  % closed, growing B

% ---- Evaluate solution and compute pH ----
t  = linspace(tspan(1), tspan(2), 601);
Y1 = deval(sol1, t);  pH1 = arrayfun(@(i) calc_pH(Y1(2,i), Y1(3,i), Y1(4,i), P), 1:numel(t));
Y2 = deval(sol2, t);  pH2 = arrayfun(@(i) calc_pH(Y2(2,i), Y2(3,i), Y2(4,i), P), 1:numel(t));


% ------ Plot pH for the 2 scenarios --------
figure;
plot(t, pH1, 'b-', 'LineWidth', 1.5); hold on;
plot(t, pH2, 'r--', 'LineWidth', 1.5); 
ylabel("pH"); 
ylim([6.5 9.5]);
xlabel("Time (h)"); 
title("pH Dynamics for Growing and Constant Biomass"); 
legend("pH-Constant Biomass", "pH-Growing Biomass")
grid on;
xlim([0 2])

%----pH vs Bacteria Dynamics----
figure;
subplot(1,2,1)
yyaxis right
plot(t, pH1, 'b-', 'LineWidth', 1.5); hold on;
ylabel("pH"); 
xlabel("Time (h)"); 
yyaxis left
plot(t, Y1(5,:).*10^9, 'r--', 'LineWidth', 1.5)
ylabel("Biomass (cells)")
title("pH Dynamics Vs Constant Biomass"); 
legend("pH Dynamics", "Constant Biomass")
grid on;

subplot(1,2,2)
yyaxis right
plot(t, pH2, 'b-', 'LineWidth', 1.5); hold on;
ylabel("pH"); 
xlabel("Time (h)"); 
yyaxis left
plot(t, Y2(5,:).*10^9, 'r--', 'LineWidth', 1.5)
ylabel("Biomass(cells)")
title("pH Dynamics Vs Growing Biomass"); 
legend("pH Dynamics", "Growing Biomass")
grid on;
hold off;



% ---- Urea, Ca2+, NH4+, and pH Dynamics ----
figure;
subplot(1,2,1); % Constant B
yyaxis left
plot(t, Y1(1,:), 'g-',  t, Y1(4,:), 'c--', t, Y1(2,:), 'm-.', 'LineWidth', 2);
ylabel("Concentration (M)");
yyaxis right
plot(t, pH1, 'b:', 'LineWidth', 2); ylabel("pH");
xlabel("Time (h)"); title("Constant Biomass"); grid on;
legend("Urea","Ca^{2+}","NH_{4}^{+}","pH", 'Location',"best");

subplot(1,2,2); % Growing B
yyaxis left
plot(t, Y2(1,:), 'g-',  t, Y2(4,:), 'c--', t, Y2(2,:), 'm-.', 'LineWidth', 2);
ylabel("Concentration (M)");
yyaxis right
plot(t, pH2, 'b:', 'LineWidth', 2); ylabel("pH");
xlabel("Time (h)"); title("Growing Biomass"); grid on;
legend("Urea","Ca^{2+}","NH_{4}^{+}","pH", 'Location',"best");

sgtitle("Urea, Ca^{2+}, NH_{4}^{+}, and pH");
 
% Local functions
function dy = odefun(~, y, isGrowing, P)
    % ODE system for closed system MICP
    U = y(1);   Ntot = y(2);   Ctot = y(3);   Ca = y(4);   B = y(5);
    % Urea hydrolysis (Michaelis-Menten)
    rate_urea = P.Vmax * B * U / (P.Km + U);
    % Solve charge balance to find [H+]
    H = calc_H(Ntot, Ctot, Ca, P);
    % Carbonate speciation and CaCO3 precipitation
    CO3 = Ctot * (P.Ka1*P.Ka2) / (H^2 + P.Ka1*H + P.Ka1*P.Ka2);
    IAP = Ca * CO3;
    rate_precip = P.k_precip * max(0, IAP - P.Ksp);   % precipitation if supersaturated
    % No CO2 outgassing in closed system (no rate_gas term)
    % ODEs
    dU  = - rate_urea;
    dN  = + 2*rate_urea;
    dC  = + rate_urea - rate_precip;
    dCa = - rate_precip;
    dB  = 0;
    if isGrowing
        dB = P.mu * B * (1 - B/P.B_max);  % logistic growth
    end
    dy = [dU; dN; dC; dCa; dB];
end

function H = calc_H(Ntot, Ctot, Ca, P)
    % Solve electroneutrality: find [H+] that balances charges
    f = @(H) 2*Ca + (Ntot/(1 + P.Ka_NH4/H)) + H ...
             - ( P.Cl0 ...
               + (Ctot * (P.Ka1*H)   / (H^2 + P.Ka1*H + P.Ka1*P.Ka2)) ...  
               + 2*(Ctot * (P.Ka1*P.Ka2)/ (H^2 + P.Ka1*H + P.Ka1*P.Ka2)) ... 
               + (1e-14 / H) );
    % Find root H in [1e-12, 1e-3] (pH 12 to 3)
    pHgrid = linspace(3, 12, 400);
    Hgrid = 10.^(-pHgrid);
    fvals = arrayfun(f, Hgrid);
    idx = find(fvals(1:end-1).*fvals(2:end) < 0, 1);
    if ~isempty(idx)
        H = fzero(f, [Hgrid(idx+1), Hgrid(idx)]);
    else
        [~, i0] = min(abs(fvals));
        H0 = Hgrid(i0);
        try
            H = fzero(f, H0);
        catch
            H = max(min(H0, 1e-3), 1e-12);
        end
    end
end

function pH = calc_pH(Ntot, Ctot, Ca, P)
    H = calc_H(Ntot, Ctot, Ca, P);
    pH = -log10(H);
end
