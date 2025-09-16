%%% === SCRIPT SETUP FOR TRAJECTORY SENSITIVITY ANALYSIS FROM POINT D === %%%
% --- Define System Parameters ---
gamma = 0.5;
omega = 0.01;
phi_FS = 0.01;
phi_NPM = 0.25;
m = 0.1;
V = 100;
W = 120;
C = 6;
K_m = 5;
f_gamma = 45;

% --- Simulation Configuration ---
% UPDATED to the newly calculated valid range for n for the Regime Shift scenario
n_values = -0.43:0.005:-0.38; 
% Define the single initial condition for all trajectories
initial_condition = [0.5, 0.4];

% --- Plotting Setup ---
options = odeset('RelTol', 1e-6, 'AbsTol', 1e-9);
figure('Units', 'inches', 'Position', [0 0 10 10]); 
hold on;
colors = cool(length(n_values));
legend_handles = gobjects(length(n_values), 1);

% --- Plot Trajectories for Each 'n' ---
for i = 1:length(n_values)
    n = n_values(i);
    
    % Create a function handle for the ODE with the current 'n'
    ode_function = @(t, y) replicator_dynamics(t, y, n);
    
    % Solve the differential equation from the single starting point 'E'
    [~, Y] = ode45(ode_function, [0 2000], initial_condition, options);
    
    % Create the legend label using LaTeX interpreter
    legend_label = sprintf('$n = %.2f$', n);
    
    % Plot the phase trajectory and store its handle for the legend
    h = plot(Y(:,1), Y(:,2), '-', 'Color', colors(i,:), 'LineWidth', 2, 'DisplayName', legend_label);
    legend_handles(i) = h;
end

% --- Plot Annotations ---
% Plot the initial starting point 'E'
plot(initial_condition(1), initial_condition(2), 'ks', 'MarkerSize', 8, 'MarkerFaceColor', 'y');
text(initial_condition(1), initial_condition(2) - 0.04, '$E$', 'FontSize', 14, 'HorizontalAlignment', 'center','Interpreter','latex');

% Plot the pure-strategy equilibria with stability colors
plot(0, 0, 'ko', 'MarkerSize', 8, 'MarkerFaceColor', 'r');
text(0-0.02, 0-0.02, '$E_1$', 'FontSize', 14, 'VerticalAlignment', 'top', 'HorizontalAlignment', 'right', 'Interpreter', 'latex');
plot(0, 1, 'ko', 'MarkerSize', 8, 'MarkerFaceColor', 'g', 'MarkerEdgeColor', 'r');
text(0-0.02, 1+0.02, '$E_2$', 'FontSize', 14, 'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'right', 'Interpreter', 'latex');
plot(1, 0, 'ko', 'MarkerSize', 8, 'MarkerFaceColor', 'g', 'MarkerEdgeColor', 'r');
text(1+0.02, 0-0.02, '$E_3$', 'FontSize', 14, 'VerticalAlignment', 'top', 'HorizontalAlignment', 'left', 'Interpreter', 'latex');
plot(1, 1, 'ko', 'MarkerSize', 8, 'MarkerFaceColor', 'r');
text(1+0.02, 1+0.02, '$E_4$', 'FontSize', 14, 'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'left', 'Interpreter', 'latex');

% --- Final Figure Formatting ---
axis([0 1 0 1]); axis square;
xlabel('Founder''s Strategy, $p$ (Prob. Concentration)','Interpreter','latex');
ylabel('Successor''s Strategy, $q$ (Prob. Manager)','Interpreter','latex');
title('Trajectory Sensitivity to Conflict Intensity $n$ from Point $E$','Interpreter','latex');
grid on; box on; xticks(0:0.1:1); yticks(0:0.1:1);
set(gca, 'FontSize', 12);
hold off;

% --- Create, Format, and Add Padding to Legend ---
lgd = legend(legend_handles, 'Interpreter', 'latex', 'Location', 'northeastoutside');
lgd.Title.String = '\hspace{.5em}Trajectory for varying $n$\hspace{.5em}';
lgd.Title.Interpreter = 'latex';

% --- Save Figure ---
ax = gca;
ax.LineWidth = 1.2;
outputFileName = 'Sensitivity_E.pdf';
exportgraphics(ax, outputFileName, 'ContentType', 'vector');
disp(['Figure saved to: ' fullfile(pwd, outputFileName)]);


%%% --- LOCAL FUNCTIONS --- %%%
function dydt = replicator_dynamics(~, y, n)
    % UPDATED parameters for Regime Shift simulation
    gamma=0.5; omega=0.01; phi_FS=0.01; phi_NPM=0.25; m=0.1;
    V=100; W=120; C=6; K_m=5; f_gamma=45;
    
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