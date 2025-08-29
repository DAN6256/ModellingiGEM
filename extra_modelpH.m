% MICP pH dynamics model simulation (robust pH solver)

% ---- Parameters ----
P.Km        = 5e-3;        % M (urease Km) 5e-3
P.Vmax      = 0.05;        % M/h per OD (urea hydrolysis max rate)
P.Ka_NH4    = 10^(-9.25);  % NH4+ <-> NH3 + H+
P.Ka1       = 10^(-6.3);   % H2CO3 <-> HCO3- + H+
P.Ka2       = 10^(-10.33); % HCO3- <-> CO3^2- + H+
P.Ksp       = 2.8e-9;      % CaCO3 solubility (mol^2/L^2)
P.k_precip  = 0.5;         % L/(mol·h) (precipitation kinetics)
P.mu        = 1.2;%4;         % 1/h (max specific growth rate) --CHANGING FROM 0.8 TO 4
P.B_max     = 5; %1.0;         % OD units (carrying capacity)

% ---Notes
%B_max = 10^10 for 1 ml volume meaning 5*10^10 for 5 ml
%B_0 = 1.2*10^9 cells

% ---- Initial conditions ----
U0    = 0.30;         % M (urea)
Ca0   = 0.30;         % M (Ca2+) Will look more here
P.Cl0 = 2*Ca0;        % M (from CaCl2 charge balance)
N0    = 0.0;          % M (NH3 + NH4+)
C0    = 1e-8;         % M (small starting DIC) --I changed a value here was -6
B0    = 0.005;%0.12;%0.1;          % OD (biomass)

y0    = [U0 N0 C0 Ca0 B0];
tspan = [0 50];        % hours
options = odeset('NonNegative',1:5);

% ---- Solve the four scenarios ----
sol1 = ode45(@(t,y) odefun(t,y,false,false,P), tspan, y0, options); % closed, constant B
sol2 = ode45(@(t,y) odefun(t,y,false,true, P), tspan, y0, options); % closed, growth
sol3 = ode45(@(t,y) odefun(t,y,true, false,P), tspan, y0, options); % open,   constant B
sol4 = ode45(@(t,y) odefun(t,y,true, true, P), tspan, y0, options); % open,   growth

% ---- Evaluate and compute pH ----
t  = linspace(tspan(1), tspan(2), 601);
Y1 = deval(sol1, t);  pH1 = arrayfun(@(i) calc_pH(Y1(2,i), Y1(3,i), Y1(4,i), P), 1:numel(t));
Y2 = deval(sol2, t);  pH2 = arrayfun(@(i) calc_pH(Y2(2,i), Y2(3,i), Y2(4,i), P), 1:numel(t));
Y3 = deval(sol3, t);  pH3 = arrayfun(@(i) calc_pH(Y3(2,i), Y3(3,i), Y3(4,i), P), 1:numel(t));
Y4 = deval(sol4, t);  pH4 = arrayfun(@(i) calc_pH(Y4(2,i), Y4(3,i), Y4(4,i), P), 1:numel(t));

% ---- Plot ----
%figure;
%yyaxis left
%plot(t, pH2, LineWidth=2);
%hold on;
%ylabel("pH dynamics")
%plot(t,Y2(5,:), LineWidth=2);
%ylabel("Bacteria")
%xlabel("Time(h)")
%legend("pH", "Bacteria")
%hold off;


figure; hold on;
plot(t.*60, pH1, 'b-',  'LineWidth',1.5); % closed, constant
plot(t.*60, pH2, 'b--', 'LineWidth',1.5); % closed, growth
%plot(t, pH3, 'r-',  'LineWidth',1.5); % open,   constant
%plot(t, pH4, 'r--', 'LineWidth',1.5); % open,   growth
xlabel('Time (h)'); ylabel('pH');
%legend('Closed, constant B','Closed, growing B','Open, constant B','Open, growing B','Location','best');
legend('Closed, constant B','Closed, growing B','Location','best');
xlim([0 70])

title('MICP pH dynamics simulation'); grid on; hold off;

% ================== Local functions ==================

function dy = odefun(~, y, isOpen, isGrowing, P)
    % y = [U, Ntot, Ctot, Ca, B]
    U    = y(1);
    Ntot = y(2);
    Ctot = y(3);
    Ca   = y(4);
    B    = y(5);

    % Ureolysis rate (Michaelis-Menten)
    rate_urea = P.Vmax * B * U/(P.Km + U);   % M/h

    % Equilibrium-based carbonate speciation to get [CO3^2-]
    H   = calc_H(Ntot, Ctot, Ca, P);         % [H+]
    CO3 = Ctot * (P.Ka1*P.Ka2) ./ (H.^2 + P.Ka1.*H + P.Ka1*P.Ka2);

    % Precipitation rate if supersaturated
    IAP = Ca * CO3;
    rate_precip = P.k_precip * max(0, IAP - P.Ksp);  % M/h

    % CO2 gas loss (open case: remove CO2 as produced by ureolysis)
    rate_gas = 0;
    if isOpen
        rate_gas = rate_urea;
    end

    % ODEs
    dU   = - rate_urea;
    dN   = + 2*rate_urea;
    dC   = + rate_urea - rate_precip - rate_gas;
    dCa  = - rate_precip;
    dB   = 0; %I have to do some investigation
    if isGrowing
        dB = P.mu * B * (1 - B/P.B_max);
    end

    dy = [dU; dN; dC; dCa; dB];
end

function pH = calc_pH(Ntot, Ctot, Ca, P)
    H  = calc_H(Ntot, Ctot, Ca, P);
    %display(H);
    pH = -log10(H);
end

function H = calc_H(Ntot, Ctot, Ca, P)
    % Electroneg: 2[Ca2+] + [NH4+] + [H+] = [Cl-] + [HCO3-] + 2[CO3^2-] + [OH-]
    f = @(H) 2*Ca + (Ntot./(1 + P.Ka_NH4./H)) + H ...
             - ( P.Cl0 ...
               +  (Ctot .* (P.Ka1.*H)   ./ (H.^2 + P.Ka1.*H + P.Ka1*P.Ka2)) ...  % [HCO3-]
               + 2*(Ctot .* (P.Ka1*P.Ka2)./ (H.^2 + P.Ka1.*H + P.Ka1*P.Ka2)) ... % 2[CO3^2-]
               + (1e-14 ./ H) );                                                % [OH-]

    % Robust bracket search over pH = 3..12
    pHmin = 3; pHmax = 12; nGrid = 400;
    Hgrid = 10.^(-linspace(pHmin, pHmax, nGrid));
    fvals = arrayfun(f, Hgrid);

    % Find sign change
    idx = find(fvals(1:end-1).*fvals(2:end) < 0, 1, 'first');

    try
        if ~isempty(idx)
            % bracket around sign change
            H = fzero(f, [Hgrid(idx+1), Hgrid(idx)]);
        else
            % no sign change: use min |f| as starting guess
            [~,i0] = min(abs(fvals));
            H0 = Hgrid(i0);
            % try unbracketed fzero; if it fails, return H0
            H = fzero(f, H0);
        end
    catch
        % fallback: minimum |f| point
        [~,i0] = min(abs(fvals));
        H = max(min(Hgrid(i0), 1e-3), 1e-12);  % clamp inside [1e-12,1e-3]
    end
end
