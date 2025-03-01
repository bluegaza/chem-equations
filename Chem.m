% MATLAB Script 

% Given Data
F_A0 = 20; % mol/s, inlet molar flow rate of A
T0 = 450; % K, inlet temperature
P0 = 10; % atm, initial pressure
V_max = 10; % dm^3, max reactor volume for plotting
Cp_A = 40; % J/mol-K
Cp_B = 25; % J/mol-K
Cp_C = 15; % J/mol-K
Delta_H_A = -70e3; % J/mol
Delta_H_B = -50e3; % J/mol
Delta_H_C = -40e3; % J/mol
R = 8.314; % J/mol-K
k_450 = 0.133; % Rate constant at 450 K
delta_H_rxn = Delta_H_B + Delta_H_C - Delta_H_A; % Reaction enthalpy
E = 31.4e3; % J/mol, activation energy

% Define rate equation
rate = @(X, T) -k_450 * exp((E/R) * (1/450 - 1/T)) * (1 - X);

% Define ODEs for PFR
pfr_odes = @(V, Y) [ -rate(Y(1), Y(2)) / F_A0;
                     (-delta_H_rxn * rate(Y(1), Y(2))) / (F_A0 * (Cp_A + Cp_B + Cp_C)) ];                 
% Solve ODEs
V_span = linspace(0, V_max, 100);
X0 = [0; T0];
[V, sol] = ode45(pfr_odes, V_span, X0);
X_sol = sol(:,1);
T_sol = sol(:,2);

% Plot results
figure;
subplot(1,2,1);
plot(V, X_sol, 'b');
xlabel('Reactor Volume (dm^3)'); ylabel('Conversion X');
title('Conversion vs Reactor Volume');

subplot(1,2,2);
plot(V, T_sol, 'r');
xlabel('Reactor Volume (dm^3)'); ylabel('Temperature (K)');
title('Temperature vs Reactor Volume');

% CSTR Analysis
X_CSTR = 0.9;
V_CSTR = F_A0 * (1 - X_CSTR) / -rate(X_CSTR, T0);
disp(['Required CSTR Volume: ', num2str(V_CSTR), ' dm^3']);

% Batch Reactor Analysis
batch_odes = @(t, Y) [ rate(1 - Y(1)/10, Y(2)) * 10;
                       (-delta_H_rxn * rate(1 - Y(1)/10, Y(2))) / (Cp_A + Cp_B + Cp_C) ];
                   
time_span = linspace(0, 10, 100);
Y0_batch = [10; T0];
[t, batch_sol] = ode45(batch_odes, time_span, Y0_batch);
N_A_sol = batch_sol(:,1);
T_batch_sol = batch_sol(:,2);

% Plot Batch Reactor Results
figure;
subplot(1,2,1);
plot(t, N_A_sol, 'b');
xlabel('Time (min)'); ylabel('Moles of A');
title('Moles of A vs Time');

subplot(1,2,2);
plot(t, T_batch_sol, 'r');
xlabel('Time (min)'); ylabel('Temperature (K)');
title('Temperature vs Time');

% Pressure Drop Analysis in PBR
rho_b = 1; % kg/dm^3
alpha = 0.019;

pbr_odes = @(V, Y) [ -rate(Y(1), Y(2)) / F_A0;
                      (-delta_H_rxn * rate(Y(1), Y(2))) / (F_A0 * (Cp_A + Cp_B + Cp_C));
                      -alpha * (Y(3) / P0) ];

Y0_pbr = [0; T0; P0];
[V, pbr_sol] = ode45(pbr_odes, V_span, Y0_pbr);
X_pbr_sol = pbr_sol(:,1);
T_pbr_sol = pbr_sol(:,2);
P_pbr_sol = pbr_sol(:,3);

% Plot PBR results
figure;
subplot(1,3,1);
plot(V, X_pbr_sol, 'b');
xlabel('Reactor Volume (dm^3)'); ylabel('Conversion X');
title('PBR Conversion vs Volume');

subplot(1,3,2);
plot(V, T_pbr_sol, 'r');
xlabel('Reactor Volume (dm^3)'); ylabel('Temperature (K)');
title('PBR Temperature vs Volume');

subplot(1,3,3);
plot(V, P_pbr_sol, 'g');
xlabel('Reactor Volume (dm^3)'); ylabel('Pressure (atm)');
title('PBR Pressure vs Volume');