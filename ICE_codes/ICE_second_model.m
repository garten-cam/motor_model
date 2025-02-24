clear variables
% High contrast colors
% First triad
teal   = '#45DEC7';
pink   = '#DE459E';
yellow = '#DED645';
% Second triad
orange = '#DE7045';
green  = '#7FDE45';
purple = '#8545DE';
figpath = "/home/cgarcia/Documents/ICE_Modelling/ICE_doc/figures/";

% -------------------------------------------------------------------------
% Modeling the internal combustion engine from dynamo tests traditional
% -------------------------------------------------------------------------

% 1. import the complete dataset 
ice_data_table = ICE_import( ...
	"ICE_data_2013_11_18@15_22_02.csv", ...
	[2 44808]);%3402 % 20s 4002
% 1.1 change the torque to Nm
% ice_data_table.Torque = ice_data_table.Torque./10;
%%
resample_period = 10;
resample_index = 1:resample_period:size(ice_data_table.BoardTime,1);%
edmd_data = struct(...
	'y', [ice_data_table.EngineRPM_A(resample_index).*2*pi/60, ice_data_table.Torque(resample_index)], ...
	'u', ice_data_table.AbsorberRPM_C(resample_index).*0.10472, ...
	't', ice_data_table.BoardTime(resample_index));
%
% Define the model
t = edmd_data.t;                     % Time vector
omega_d_measured = edmd_data.u / max(edmd_data.u); % Measured dynamometer RPM (rad/s)
T_d_measured = edmd_data.y(:,2) / max(edmd_data.y(:,2));         % Measured absorber torque (Nm)
omega_e_measured = edmd_data.y(:,1) / max(edmd_data.y(:,1)); % Measured engine RPM (rad/s)


% Initial parameter guesses
tau1_init = 0.5;  % Engine response time constant
tau2_init = 5;  % Torque response time constant
k1_init = 1.0;    % Engine speed gain
k2_init = 10.0;   % Torque gain

params_init = [tau1_init, tau2_init, k1_init, k2_init];

% Set constraints (time constants > 0 for stability)
lb = [1e-3, 1e-3, 0, 0];  
ub = [Inf, Inf, Inf, Inf];  

% Optimize using fmincon
options = optimoptions('fmincon', 'Algorithm', 'sqp', 'MaxFunctionEvaluations', 10000);
opt_params = fmincon(@(p) cost_function(p, t, omega_d_measured, omega_e_measured, T_d_measured), ...
                     params_init, [], [], [], [], lb, ub, [], options);

% Extract optimized values
tau1_opt = opt_params(1);
tau2_opt = opt_params(2);
k1_opt = opt_params(3);
k2_opt = opt_params(4);

% fprintf('Optimized Parameters:\ntau1 = %.4f, tau2 = %.4f, k1 = %.4f, k2 = %.4f\n', tau1_opt, tau2_opt, k1_opt, k2_opt);

% Simulate with optimized parameters
omega_d_interp = @(tq) interp1(t, omega_d_measured, tq, 'linear', 'extrap');
ode_fun_opt = @(t, u) engine_model(t, u, opt_params, omega_d_interp);
[t_sim, u_sim] = ode45(ode_fun_opt, t, [omega_e_measured(1); T_d_measured(1)]);

omega_e_sim = u_sim(:,1);
T_d_sim = u_sim(:,2);
%%
% Plot results
trad_appx_fig = figure(1);
clf
colororder({teal, orange})
tiledlayout("vertical","TileSpacing","tight","Padding","tight")
nexttile
hold on
yyaxis left
plot(t, omega_e_measured, "Color",teal, 'LineWidth', 1.5);
plot(t, omega_e_sim, '-.',"Color",pink, 'LineWidth', 1.5);
ylabel('Engine [rad/s] (normalized)', "Interpreter","latex");
xlabel("time [s]")
yyaxis right
plot(t, T_d_measured, "Color",orange, 'LineWidth', 1.5); hold on;
plot(t, T_d_sim, '-.',"Color",green, 'LineWidth', 1.5);
xlabel('time [s]',"Interpreter","latex");
ylabel('Torque [Nm] (normalized)',"Interpreter","latex");
legend('Speed', 'Speed approximation','Torque', 'Torque approximation',"Location","northoutside");
nexttile
plot(t, omega_d_measured,"k",LineWidth=1.5)
xlabel("time [s]","Interpreter","latex")
ylabel("Absorber [rad/s] (normalized)","Interpreter","latex")
set(gcf,"PaperPosition",[0,0,10,13])
saveas(trad_appx_fig,strcat(figpath, "trad_appx_fig.fig"))
saveas(trad_appx_fig,strcat(figpath, "trad_appx_fig.eps"), "epsc")
%%
function dudt = engine_model(t, u, params, omega_d_interp)
    % Extract parameters
    tau1 = params(1);
    tau2 = params(2);
    k1 = params(3);
    k2 = params(4);
    
    % Get interpolated absorber RPM (input)
    omega_d = omega_d_interp(t);
    
    % State variables
    omega_e = u(1);
    T_d = u(2);
    
    % ODE equations with separate time constants
    d_omega_e = (k1 * omega_d - omega_e) / tau1;
    d_T_d = (k2 / omega_d - T_d) / tau2;
    
    % Output derivative
    dudt = [d_omega_e; d_T_d];
end

function error = cost_function(params, t, omega_d_measured, omega_e_measured, T_d_measured)
    % Interpolate absorber RPM to match solver time steps
    omega_d_interp = @(tq) interp1(t, omega_d_measured, tq, 'linear', 'extrap');
    
    % Initial conditions (first measured values)
    u0 = [omega_e_measured(1); T_d_measured(1)];
    
    % Solve ODE
    ode_fun = @(t, u) engine_model(t, u, params, omega_d_interp);
    [t_sim, u_sim] = ode45(ode_fun, t, u0);
    
    % Extract simulated values
    omega_e_sim = u_sim(:,1);
    T_d_sim = u_sim(:,2);
    
    % Compute cost (sum of squared errors)
    error = sum((omega_e_sim - omega_e_measured).^2) + sum((T_d_sim - T_d_measured).^2);
end