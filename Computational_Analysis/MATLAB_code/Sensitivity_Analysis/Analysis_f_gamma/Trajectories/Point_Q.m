%%% === SCRIPT SETUP FOR TRAJECTORY SENSITIVITY FROM POINT Q (CONFLICTUAL CASE) === %%%
% --- Define System Parameters (from Table "params_extended") ---
% These are the fixed parameters for the "regime shift" scenario.
gamma = 0.5;
omega = 0.01;
phi_FS = 0.01;
phi_NPM = 0.25;
m = 0.1;
n = -0.4;
V = 100;
W = 120;
C = 6;
K_m = 5;

% --- Simulation Configuration ---
% Define the range of the parameter 'f_gamma' to be analyzed.
% The theoretical range for a valid interior saddle point in this regime is (41, 48).
f_gamma_values = 42:0.5:47; 
% Define the single initial condition for all trajectories
initial_condition = [0.5, 0.5]; % Starting point Q

% --- Plotting Setup ---
options = odeset('RelTol', 1e-6, 'AbsTol', 1e-9);
figure('Units', 'inches', 'Position', [0 0 10 10]); 
hold on;
colors = cool(length(f_gamma_values));
legend_handles = gobjects(length(f_gamma_values), 1);

% --- Plot Trajectories for Each 'f_gamma' ---
for i = 1:length(f_gamma_values)
    f_gamma = f_gamma_values(i);
    
    % Create a function handle for the ODE with the current 'f_gamma'
    ode_function = @(t, y) replicator_dynamics_Q_conflict(t, y, f_gamma);
    
    % Solve the differential equation from the single starting point 'Q'
    [~, Y] = ode45(ode_function, [0 2000], initial_condition, options);
    
    % Create the legend label using LaTeX interpreter
    legend_label = sprintf('$f_{\\gamma} = %.1f$', f_gamma);
    
    % Plot the phase trajectory and store its handle for the legend
    h = plot(Y(:,1), Y(:,2), '-', 'Color', colors(i,:), 'LineWidth', 2, 'DisplayName', legend_label);
    legend_handles(i) = h;
end

% --- Plot Annotations ---
% Plot the initial starting point 'Q'
plot(initial_condition(1), initial_condition(2), 'ks', 'MarkerSize', 8, 'MarkerFaceColor', 'y');
text(initial_condition(1), initial_condition(2) - 0.04, '$Q$', 'FontSize', 14, 'HorizontalAlignment', 'center','Interpreter','latex');

% Plot the pure-strategy equilibria with correct stability colors for the conflictual case
% E2 and E3 are stable (green), E1 and E4 are unstable (red)
plot(0, 0, 'ko', 'MarkerSize', 8, 'MarkerFaceColor', 'r', 'MarkerEdgeColor', 'r');
text(0-0.02, 0-0.02, '$E_1$', 'FontSize', 14, 'VerticalAlignment', 'top', 'HorizontalAlignment', 'right', 'Interpreter', 'latex');
plot(0, 1, 'ko', 'MarkerSize', 8, 'MarkerFaceColor', 'g', 'MarkerEdgeColor', 'g');
text(0-0.02, 1+0.02, '$E_2$', 'FontSize', 14, 'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'right', 'Interpreter', 'latex');
plot(1, 0, 'ko', 'MarkerSize', 8, 'MarkerFaceColor', 'g', 'MarkerEdgeColor', 'g');
text(1+0.02, 0-0.02, '$E_3$', 'FontSize', 14, 'VerticalAlignment', 'top', 'HorizontalAlignment', 'left', 'Interpreter', 'latex');
plot(1, 1, 'ko', 'MarkerSize', 8, 'MarkerFaceColor', 'r', 'MarkerEdgeColor', 'r');
text(1+0.02, 1+0.02, '$E_4$', 'FontSize', 14, 'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'left', 'Interpreter', 'latex');

% --- Final Figure Formatting ---
axis([0 1 0 1]); axis square;
xlabel('Founder''s Strategy, $p$ (Prob. Concentration)','Interpreter','latex');
ylabel('Successor''s Strategy, $q$ (Prob. Manager)','Interpreter','latex');
title('Trajectory Sensitivity from Point $Q$ to Varying $f_{\gamma}$ (Conflictual Case)','Interpreter','latex');
grid on; box on; xticks(0:0.1:1); yticks(0:0.1:1);
set(gca, 'FontSize', 12);
hold off;

% --- Create, Format, and Add Padding to Legend ---
lgd = legend(legend_handles, 'Interpreter', 'latex', 'Location', 'northeastoutside');
lgd.Title.String = '\hspace{.5em}Trajectory for varying $f_{\gamma}$\hspace{.5em}';
lgd.Title.Interpreter = 'latex';

% --- Save Figure ---
ax = gca;
ax.LineWidth = 1.2;
outputFileName = 'Sensitivity_Q.pdf';
exportgraphics(ax, outputFileName, 'ContentType', 'vector');
disp(['Figure saved to: ' fullfile(pwd, outputFileName)]);

%%% --- LOCAL FUNCTIONS --- %%%
function dydt = replicator_dynamics_Q_conflict(~, y, f_gamma)
    % System Parameters for Conflictual Case (from Table "params_extended")
    gamma = 0.5;
    omega = 0.01;
    phi_FS = 0.01;
    phi_NPM = 0.25;
    m = 0.1;
    n = -0.4;
    V = 100;
    W = 120;
    C = 6;
    K_m = 5;
    
    p = y(1);
    q = y(2);
    dydt = zeros(2,1);
    
    % Replicator dynamics equation for p
    dydt(1) = p*(1 - p)*( q*(phi_FS - n)*V - (1 - q)*n*W - f_gamma );
    
    % Replicator dynamics equation for q
    X_FSM = ( p + (1 - p)*gamma*(1 - omega - phi_FS + n) )*V - C;
    X_FSS = ( p*(1 - omega - phi_NPM + m) + (1 - p)*gamma*(1 - omega - phi_NPM + m + n) )*W - K_m;
    dydt(2) = q*(1 - q)*( X_FSM - X_FSS );
end
