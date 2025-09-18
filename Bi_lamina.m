
% Copyright (c) 2025 Ahmad Reza Ravangard. All rights reserved.
% Citation required: "Bilamina Analytical Simulation — Ahmad Reza Ravangard. Code available at: https://github.com/arava001/Bilamina_Analytical_Simulation"
% No use, copying, modification or distribution without written permission.

clc; clear; close all 

filename = 'temp_cycle.csv';
[t,T] = exp_temp(filename, 1);
t = t*60;
global alpha_init 
alpha_init = 0.8;
step = 30;
% time_ = [0 11 145.333 290] * 60;                 
% Temp_ = [75 175 175 34.8];                            

DOC = computeDOC(t, T);
%% Neat resin
% [DOC,t,T] = Kamal_Sourour_model_cure_kinetic (Temp_,time_); % t[sec] T[C]
[Tg] = glass_transition_temp (DOC); %Tg[C]



figure Name 'tt'
hold on 
plot(t/60,T)

% Plots
figure Name 'results1';
subplot(2,3,1)
yyaxis right
plot(t / 60, DOC,  'Color', 'r', 'LineWidth', 2);
xlabel('Time [min]');
ylabel('Degree of Cure (\alpha) [-]');
grid on;
ax = gca;
ax.YAxis(2).Color = 'r';

yyaxis left 
plot(t / 60, T, 'b', 'LineWidth', 2); % Temperature
hold on;
plot(t / 60, Tg, 'b--', 'LineWidth', 2); % Tg
ylabel('T  and T_g [°C]' );
legend( 'Temperature (T)', 'Glass Transition Temperature (Tg)','', 'Location', 'best');
set(gca, 'FontSize',14)

subplot(2,3,2)
plot(DOC,Tg,'LineWidth',2)
xlabel('Degree of Cure (\alpha) [-]')
ylabel ('T_g [C]')
grid on 
set(gca, 'FontSize',14)


beta = 0.03; % Matrix chemical shrinkage 
z = T-Tg; % does not matter[C or K]
ind_b = find (z<0);  %below Tg is chemical shrinkage 
ind_a = find (z>=0); % abot Tg both themal expansion and chemical shrinkage happen


CTE_b = NaN+0.*DOC;
CTE_a = NaN+0.*DOC;

% CTE_b (ind_b) = (-451.*DOC(ind_b) + 619.54).*10^(-6) ; %[1/C]
CTE_b (ind_b) = (-410.*DOC(ind_b) + 619.54).*10^(-6) ; %[1/C]

CTE_a (ind_a) = 84.6e-6;


CTE = NaN+0.*DOC;
CTE(ind_b) = CTE_b (ind_b);
CTE(ind_a) = CTE_a (ind_a);

subplot(2,3,3)
yyaxis right 
plot(t/60, T,'r','LineWidth',2);
ylabel('Temperature [C]');
grid on;

yyaxis left 
plot(t/60, z,'b','LineWidth',2);
xlabel('Time [min]')
ylabel('T-Tg [C]');
set(gca, 'FontSize',14)
ax = gca;
ax.YAxis(2).Color = 'r';

subplot(2,3,4)
yyaxis left 
plot(DOC(ind_a),T(ind_a),'b',LineWidth=2)
hold on 
plot(DOC(ind_a),Tg(ind_a),'b--',LineWidth=2)
xlabel ('Degree of Cure (\alpha)')
ylabel ('T & Tg (C)')
grid on 
hold off

yyaxis right
% plot(DOC, CTE_b,'*','Linewidth', 3)
% hold on
plot(DOC(ind_a), CTE_a(ind_a),'r','Linewidth', 1)
ylim([mean(CTE_a(ind_a))-0.00001;mean(CTE_a(ind_a))+0.00001])
grid on 
legend ('T', 'Tg','above Tg','','',Location='best')
xlabel('Degree of Cure (\alpha) [-]')
ylabel('CTE [1/C]')
ax = gca;
ax.YAxis(2).Color = 'r';
set(gca, 'FontSize',14)

subplot(2,3,5)
yyaxis left 
plot(DOC(ind_b),T(ind_b),'b',LineWidth=2)
hold on 
plot(DOC(ind_b),Tg(ind_b),'b--',LineWidth=2)
xlabel ('Degree of Cure (\alpha)')
ylabel ('T & Tg (C)')
grid on 
hold off

yyaxis right
% plot(DOC, CTE_b,'*','Linewidth', 3)
% hold on
plot(DOC(ind_b), CTE_b(ind_b),'r','Linewidth', 1)
ylim([mean(CTE_b(ind_b))-0.00001;mean(CTE_b(ind_b))+0.00001])
grid on 
legend ('T', 'Tg','below Tg','','',Location='best')
xlabel('Degree of Cure (\alpha) [-]')
ylabel('CTE [1/C]')
ax = gca;
ax.YAxis(2).Color = 'r';
set(gca, 'FontSize',14)



eps_ch = beta*DOC; %[-]

subplot(2,3,6)
plot(DOC,eps_ch,'k',LineWidth=2)
xlabel ('Degree of Cure (\alpha)')
ylabel ('$\varepsilon^{ch}$',Interpreter='latex',FontSize=27)
grid on 
hold off
set(gca, 'FontSize',14)
set(gcf, 'Units', 'normalized', 'OuterPosition', [0 0 1 1]); % Set the figure to full screen
%===============================================================
E0_1  = 211.7;       % [MPa]
E0_2  = 329.14;      % [MPa]
delta_T_1 = 35.29;   % [C]
delta_T_2 = -17.95;  % [C]
tau_1 = 3.97;        % [C]
tau_2 = 11.46;       % [C]

mem_E_1 = E0_1./(1+exp(((T-Tg)-delta_T_1)./(tau_1))); % [MPa]
mem_E_2 = E0_2./(1+exp(((T-Tg)-delta_T_2)./(tau_2))); % [MPa]
EH_prime = mem_E_1+mem_E_2; % [MPa]

EC_prime = -6.58*T+3701; % [MPa]

tc = 300*60; % time that cooling begins [Sec]
c = 0.02*60; % for smooth transition of E' from isothermal to cooling

E_prime = (EC_prime-EH_prime)./(1+exp((tc-t)/c)); % [MPa]

E_prime = E_prime+EH_prime; % [MPa]

%% Lamina 
% Given material constants from the table in the paper
E1f = 275e+3;       % Longitudinal Young's modulus of fiber [MPa]
E2f = 153e+3;        % Transverse Young's modulus of fiber [MPa]
G12f = 26e+3;       % Shear modulus of fiber [MPa]

v12f = 0.26;     % Major Poisson's ratio of fiber
v23f = 0.26;     % Poisson's ratio in the transverse plane of fiber
alpha1f = -1e-6; % CTE in longitudinal direction of fiber [1/°C]
alpha2f = 15e-6; % CTE in transverse direction of fiber [1/°C]

% Additional constants inferred from the document text
vm = 0.36;       % Poisson's ratio of matrix (resin)
ff = 0.6;        % Fiber volume fraction
fm = 1 - ff;     % Matrix volume fraction
Em = E_prime;    % Assumed value for the matrix modulus [MPa]
alpha_m = CTE;   % CTE of matrix [1/°C]
beta_m = beta;   % Chemical shrinkage coefficient of matrix

% calculations 

Gm = Em ./ (2 .* (1 + vm));  % [MPa]

G23f = E2f / (2 * (1 + v23f));   % [MPa]    

km = Em ./ (2 .* (1 - vm) - 4 .* vm.^2);  % [MPa]

kf = (E1f * 15) ./ (2 * (1 - 0.26) * E1f - 4 * 0.26^2 * 15); % [MPa]

k2 = ((kf + Gm) .* km + (kf - km) .* Gm * ff) ./ ((kf + Gm) - (kf - km) .* ff);  % [MPa]

E1 = E1f * ff + Em .* fm + (4 * (vm - v12f).^2 .* km .* kf .* Gm ./ ((kf + Gm) .* km + (kf - km) .* Gm .* ff)) * ff * fm;% [MPa]

G12 = Gm .* ((G12f + Gm) + (G12f - Gm) .* ff) ./ ((G12f + Gm) - (G12f - Gm) .* ff);  % [MPa]
G13 = G12; % [MPa]

G23 = Gm .* ((G23f + Gm) .* km + 2 .* G23f .* Gm + (G23f - Gm) .* km * ff) ./ ((G23f + Gm) .* km + 2 .* G23f .* Gm - (G23f - Gm) .* (km + 2 * Gm) * ff);% [MPa]

v12 = v12f * ff + vm .* fm + ((vm - v12f) .* (km - kf) .* Gm ./ ((kf + Gm) .* km + (kf - km) .* Gm * ff)) .* ff .* fm;
v13 = v12; 

E2 = 1 ./ (1./(4 .* k2) + 1/(4 * G23f) + (v12f^2) ./ E1f);  % [MPa]
E3 = E2;% [MPa]

v23 = 1 - (E2 ./ (2 .* k2)) - 2 * v12f^2 .* (E2 / E1f);

alpha2 = (alpha2f + v12f * alpha1f) * ff + (1 + vm) .* alpha_m .* fm - (v12f * ff + vm .* fm) * alpha1f;
alpha3 = alpha2;
alpha1 = (alpha1f * E1f * ff + alpha_m .* Em .* fm) ./ (E1f * ff + Em .* fm);  % Element-wise

beta2 = (1 + vm) .* beta_m .* fm - v12f * beta_m;
beta3 = beta2;
beta1 = (beta_m .* Em .* fm) ./ (E1f * ff + Em .* fm);  % Element-wise

% ================================================================================================
% % Calculation of Delta (Δ')
% v21 = v12 .* (E2 ./ E1);
% v32 = v23 .* (E3 ./ E2);
% v31 = v13 .* (E3 ./ E1);
% 
% Delta = (1 - v12 * v21 - v23 * v32 - v31 .* v13 - 2 * v21 .* v32 .* v13) ./ (E1 * E2 * E3);
% 
% % Homogenized stiffness matrix elements
% C11 = (1 - v23^2) / (E2 * E3 * Delta);
% C22 = (1 - v13 * v31) / (E1 * E3 * Delta);
% C33 = (1 - v12 * v21) / (E1 * E2 * Delta);
% 
% C12 = (v21 + v31 * v23) / (E2 * E3 * Delta);
% C13 = (v31 + v21 * v32) / (E2 * E3 * Delta);
% C23 = (v32 + v12 * v31) / (E1 * E2 * Delta);
% 
% C44 = G23;
% C55 = G13;
% C66 = G12;
% ================================================================================================
%% Laminate 
% calculating the thickness h1(t) and curvature
%input 
h1_0 = 1*10^(-3); % [m], initial thickness of h1
h1 = calculate_thickness(h1_0, beta, DOC, t, T, CTE);
% =======================
L =  0.1016;         % [m], strip length
h2 = 0.125*10^(-3);  % [m], thickness of substrate
w1 = 0.0254;         % [m], sample width
w2 = 0.0254;         % [m], substrate width
E2 = 135e+3;         % [MPa] substrate   

w = w1;              % m, sample width
I2 = w*h2^3/12;      % [m^4]    
A2 = w*h2;           % [m^2] 
I1 = w/12.*(h1).^3;
A1 = w.*h1; 
h = h1+h2;
b = beta*ones(length(DOC));


dens = 1550;

q = w*h*dens*9.80665;
C1 = -q*L^3/6;
C2 = 3/24*q*L^4;
delta_W = (1 ./ (1e+6*E1 .* I1 + 1e+6*E2 * I2)) .* (q * L^4 / 8);
x= 0:L/20:L;
w_weight = -1./(1e+6*E1.*I1+1e+6*E2*I2).*(q.*x.^4/24 + C1.*x + C2);

% curvature 
KC = 0.*DOC;
KT = 0.*DOC;
for i = 1:length(t)-1
    t1 = t(i); 
    t2 = t(i+1);
    a1 = DOC(i);
    a2 = DOC(i+1);
    T1 = T(i);
    T2 = T(i+1);
    dt = t2-t1;
    da = a2-a1;
    dT = T2-T1;
    da_dt = da/dt;
    dT_dt = dT/dt;

    tt1 = b(i);
    tt2 = b(i+1);
    mm1 = -CTE(i);
    mm2 = -CTE(i+1);


    denom1 = 2/h(i);
    denom1 = denom1*(E1(i)*I1(i)+E2*I2); 
    denom1 = denom1*((1/(E1(i)*A1(i)))+(1/(E2*A2)));
    denom1 = denom1+h(i)/2;
    z1 = tt1/denom1; 
    l1 = mm1/denom1; 

    denom2 = 2/h(i+1);
    denom2 = denom2*(E1(i+1)*I1(i+1)+E2*I2); 
    denom2 = denom2*((1/(E1(i+1)*A1(i+1)))+(1/(E2*A2)));
    denom2 = denom2+h(i+1)/2;
    z2 = tt2/denom2; 
    l2 = mm2/denom2;
    nnn =10;
    ds1 = linspace(a1,a2,nnn); 
    ds2 = linspace(T1,T2,nnn);  
    dz = linspace(z1,z2,nnn);
    dl = linspace(l1,l2,nnn);

    dkc = trapz(ds1,dz);
    KC(i+1) = KC(i) + dkc;
    dkt = trapz(ds2,dl);
    KT(i+1) = KT(i)+dkt;    
end
K_total = KC+KT; 
delta_CP = K_total*L^2./2;
delta = delta_CP-delta_W;

%========================================================================
figure Name 'results3'
subplot(2,3,1)
plot(DOC,E_prime,LineWidth=2)
xlabel('DOC')
ylabel("E' [MPa]")
grid on 
set(gca, 'FontSize',14)

subplot(2,3,2)
hold on
plot(DOC,E1/1000,LineWidth=2)
xlabel('DOC')
ylabel("E [GPa]")
legend('E1')
grid on 
set(gca, 'FontSize',14)

subplot(2,3,3)
hold on
plot(DOC,E2/1000,LineWidth=2)
xlabel('DOC')
ylabel("E [GPa]")
legend('E2')
grid on 
set(gca, 'FontSize',14)

subplot(2,3,4)
hold on
plot(DOC,Gm,LineWidth=2)
xlabel('DOC')
ylabel("Gm [MPa]")
grid on 
set(gca, 'FontSize',14)

subplot(2,3,5)
hold on
plot(DOC,G12/1000,LineWidth=2)
xlabel('DOC')
ylabel("G [GPa]")
legend('G12')
grid on 
set(gca, 'FontSize',14)

subplot(2,3,6)
hold on
plot(DOC,G23/1000,LineWidth=2)
xlabel('DOC')
ylabel("G [GPa]")
legend('G23')
grid on 
set(gca, 'FontSize',14)
set(gcf, 'Units', 'normalized', 'OuterPosition', [0 0 1 1]); % Set the figure to full screen
%========================================================================
figure Name 'results4'

subplot(2,3,1)
plot(DOC,alpha1,LineWidth=2)
xlabel('DOC')
ylabel("\alpha_1")
grid on 
set(gca, 'FontSize',14)

subplot(2,3,2)
plot(DOC,alpha2,LineWidth=2)
xlabel('DOC')
ylabel("\alpha_2")
grid on 
set(gca, 'FontSize',14)

subplot(2,3,3)
plot(DOC,beta1,LineWidth=2)
xlabel('DOC')
ylabel("\beta_1")
grid on 
set(gca, 'FontSize',14)
set(gcf, 'Units', 'normalized', 'OuterPosition', [0 0 1 1]); % Set the figure to full screen

%========================================================================
figure Name 'results5'

subplot(2,3,1)
yyaxis left 
plot (t/60,T)
ylabel('Temperature [C]')

yyaxis right 
plot(t/60,h*1000,LineWidth=2)
xlabel('Time [min]')
ylabel("h [mm]")
grid on 
set(gca, 'FontSize',14)

subplot(2,3,2)
yyaxis left 
plot (t/60,T)
ylabel('Temperature [C]')

yyaxis right 
plot(t/60,KC,LineWidth=2)
xlabel('Time [min]')
ylabel("K_{CH} [1/m]")
grid on 
set(gca, 'FontSize',14)



subplot(2,3,3)
yyaxis left 
plot (t/60,T)
ylabel('Temperature [C]')

yyaxis right 
plot(t/60,KT,LineWidth=2)
xlabel('Time [min]')
ylabel("K_{CTE} [1/m]")
grid on 
set(gca, 'FontSize',14)

subplot(2,3,4)
yyaxis left 
plot (t/60,T)
ylabel('Temperature [C]')

yyaxis right 
plot(t/60,KT,LineWidth=2)
xlabel('Time [min]')
ylabel("K_{total} [mm]")
grid on 
set(gca, 'FontSize',14)
set(gcf, 'Units', 'normalized', 'OuterPosition', [0 0 1 1]); % Set the figure to full screen



%% Experimental 

table = readtable ("Time vs deflection_test5.xlsx");


figure Name 'result6'
subplot(2,3,1)
yyaxis left 
plot (t/60,T,LineWidth=2)
ylabel('Temperature [C]')

yyaxis right
plot(t/60,delta_CP*1000,LineWidth=2)
xlabel('Time [min]')
ylabel("\delta_{CP} [mm]")
grid on 
set(gca, 'FontSize',14)

subplot(2,3,2)
yyaxis left 
plot (t/60,T,LineWidth=2)
ylabel('Temperature [C]')

yyaxis right
plot(t/60,delta_W*1000,LineWidth=2)
xlabel('Time [min]')
ylabel("\delta_{W} [mm]")
grid on 
set(gca, 'FontSize',14)

subplot(2,3,3)
yyaxis left 
plot (t/60,T,LineWidth=2)
ylabel('Temperature [C]')

yyaxis right
plot(t/60,delta*1000,LineWidth=2)
hold on 


plot(table.Time_min_(1:step:end), table.Deflection_mm_(1:step:end), '--*');

xlabel('Time [min]')
ylabel("\delta [mm]")
legend ('','model','experiment')
grid on 
set(gca, 'FontSize',14)

subplot(2,3,5)
plot(t / 60, T, 'b', 'LineWidth', 2); % Temperature
hold on;
plot(t / 60, Tg, 'b--', 'LineWidth', 2); % Tg
ylabel('T  and T_g [°C]' );
legend( 'Temperature (T)', 'Glass Transition Temperature (Tg)','', 'Location', 'best');
set(gca, 'FontSize',14)
set(gcf, 'Units', 'normalized', 'OuterPosition', [0 0 1 1]); % Set the figure to full screen


subplot(2,3,6)
plot(t / 60, DOC, 'b', 'LineWidth', 2); % Temperature

set(gca, 'FontSize',14)
set(gcf, 'Units', 'normalized', 'OuterPosition', [0 0 1 1]); % Set the figure to full screen

%% Function definitions


function h1 = calculate_thickness(h1_0, beta, DOC, t, T, CTE)
    % Calculates the thickness h1 over time based on chemical shrinkage
    % and thermal expansion/contraction.
    %
    % Inputs:
    %   h1_0 - Initial thickness of h1 (m)
    %   beta - Chemical shrinkage coefficient of the lamina
    %   DOC  - Degree of cure vector over time
    %   t    - Time vector
    %   T    - Temperature vector over time
    %   CTE  - Coefficient of thermal expansion vector over time
    %
    % Output:
    %   h1   - Thickness vector over time

    % Initialize thickness and change in thickness arrays
    H = zeros(length(DOC), 1);
    dH = zeros(length(DOC), 1);
    
    % Set initial thickness
    H(1) = h1_0;
    
    % Loop over each time step to calculate thickness changes
    for i = 2:length(DOC)
        dt = t(i) - t(i-1);
        
        if dt == 0
            disp(['dt=0 at index ', num2str(i)]);
            doc_dt = 0;
            dT_dt = 0;
        else
            doc_dt = (DOC(i) - DOC(i-1)) / dt;
            dT_dt = (T(i) - T(i-1)) / dt;
        end
        
        % Compute the integral term
        Int = CTE(i) * dT_dt - beta * doc_dt;
        
        % Change in thickness over this time step
        dH(i) = h1_0 * Int * dt;
        
        % Check for NaN values
        if isnan(dH(i))
            disp(['dH is NaN at index ', num2str(i)]);
        end
        
        % Update thickness
        H(i) = H(i-1) + dH(i);
    end
    
    % Output the calculated thickness
    h1 = H;
end





function [Tg_values] = glass_transition_temp (DOC)
    % Parameters for calculating Tg
    Tg0 = -27.1;       % Initial Tg of uncured resin
    Tg_inf = 125;   % Ultimate Tg of fully cured resin
    lambda =  0.392;   % Fitting parameter for Tg calculation
    
    % Calculate Tg as a function of degree of cure
    Tg_values = Tg0 + ((Tg_inf - Tg0) * lambda .* DOC) ./ (1 - (1 - lambda) .* DOC);
end

function [alpha_total,t_total,temperature] = Kamal_Sourour_model_cure_kinetic (Temp_,time_)
% Kamal-Sourour model parameters
A1 = 0;
E1 = 0;
A2 = 6879;
E2 = 53878;
m = 0.32;
n = 1.66;
R = 8.314;
alpha0 = 0.7;


% A1 = 0;
% E1 = 0;
% A2 = 24.58;
% E2 = 34.076e+3;
% m = 0.32;
% n = 1.66;
% R = 8.314;
% alpha0 = 0.01;

time_spans = {[time_(1), time_(2)], [time_(2), time_(3)], [time_(3), time_(4)]};
temperatures = [Temp_(1), Temp_(2), Temp_(3)];

t_total = [];
alpha_total = [];
temperature =[];
for i = 1:length(time_spans)
    tspan = time_spans{i};
    TEMP = temperatures(i);

    t_ini = Temp_(i);
    t_fin =Temp_(i+1);

    cure_kinetics = @(t, alpha) (A1 * exp(-E1 / (R * (TEMP + 273.15))) + ...
                                 A2 * exp(-E2 / (R * (TEMP + 273.15))) * alpha^m) * (1 - alpha)^n;
    [t_segment, alpha_segment] = ode45(cure_kinetics, tspan, alpha0);
    t_total = [t_total; t_segment];
    alpha_total = [alpha_total; alpha_segment];
    alpha0 = alpha_segment(end);
    t1 = temp_make(tspan, [t_ini t_fin],length(alpha_segment));
    temperature = [temperature t1];
end
temperature = temperature';
end


function [time,temperature] = exp_temp(filename, flag)
    if flag
        tab = readtable(filename);
        timeArray = datetime(tab.Time, 'InputFormat', 'hh:mm:ss a');
        intervalsInMinutes = minutes(timeArray - timeArray(1));
        figure('Name', 'tt');
        subplot (1, 2, 1)
        plot(intervalsInMinutes, tab.Temp_C_, 'LineWidth', 2);
        grid on;
        xlabel('Time [min]');
        ylabel('T [°C]');
        set(gca, 'FontSize', 18);
        ind2 = (find(intervalsInMinutes <= 2.833333333333333));
        intervalsInMinutes(ind2) = [];
        intervalsInMinutes = intervalsInMinutes - intervalsInMinutes(1);
        T = tab.Temp_C_; 
        T(ind2) = [];
        temperature = T;
        time = intervalsInMinutes;
        subplot (1, 2, 2)
        plot(intervalsInMinutes, T, 'LineWidth', 2);
        grid on;
        xlabel('Time [min]');
        ylabel('T [°C]');
        set(gca, 'FontSize', 18);
        rate = diff(T); 
        intervalsForRate = intervalsInMinutes(2:end);
        ind = find(abs(rate) > 1);
        significantTimes = intervalsForRate(ind);
        significantTemps = T(ind + 1);
        significantRates = rate(ind);
        disp(table(significantTimes, significantTemps, significantRates, ...
                   'VariableNames', {'Time_min', 'Temperature_C', 'Rate_C_per_min'}));
    end
end

function [T] = temp_make(time_span, Temp_span,m)
    time = linspace(time_span(1), time_span(2), m); 
    rate = Temp_span(2) - Temp_span(1);
    rate = rate / (time_span(2) - time_span(1)); 
    T = rate * (time - time_span(1)) + Temp_span(1);   
end

function alpha_resized = computeDOC(t, T)
    % Compute Degree of Cure (DOC) based on given time and temperature arrays.
    % Inputs:
    %   t - Time array (seconds)
    %   T - Temperature array (Celsius, same size as t)
    % Output:
    %   alpha_resized - Degree of Cure (DOC) at discrete time points in t
    global alpha_init
    % Parameters
    A1 = 0; % Pre-exponential factor 1
    E1 = 0; % Activation energy 1
    A2 = 6879; % Pre-exponential factor 2
    E2 = 53878; % Activation energy 2
    m = 0.32; % Reaction order parameter
    n = 1.66; % Reaction order parameter
    R = 8.314; % Universal gas constant
    alpha0 = alpha_init; % Initial condition for alpha

    dt = 0.01; 
    t_fine = t(1):dt:t(end); 

    alpha_fine = zeros(size(t_fine)); 
    alpha_fine(1) = alpha0; 

    for i = 1:length(t_fine)-1

        idx = find(t <= t_fine(i), 1, 'last'); 
        TEMP_current = T(idx);

        d_alpha_dt = (A1 * exp(-E1 / (R * (TEMP_current + 273.15))) + ...
                      A2 * exp(-E2 / (R * (TEMP_current + 273.15))) * alpha_fine(i)^m) * (1 - alpha_fine(i))^n;

        alpha_fine(i+1) = alpha_fine(i) + dt * d_alpha_dt;
    end

    alpha_resized = interp1(t_fine, alpha_fine, t, 'pchip'); % Use 'pchip' interpolation
end
