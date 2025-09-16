%%% === SCRIPT SETUP FOR TRAJECTORY SENSITIVITY FROM POINT L === %%%
% --- Define System Parameters (from Table "params_extended") ---
% These are the fixed parameters for the "regime shift" scenario.
gamma = 0.5;
omega = 0.01;
phi_FS = 0.01;
m = 0.1;
n = -0.4;
V = 100;
W = 120;
C = 6;
K_m = 5;
f_gamma = 45;

% --- Simulation Configuration ---
% Define the range of the parameter 'phi_NPM' to be analyzed, based on the manuscript.
phi_NPM_values = 0.23:0.002:0.26; 
% Define the single initial condition for all trajectories
initial_condition = [0.9, 0.4]; % Starting point L

% --- Plotting Setup ---
options = odeset('RelTol', 1e-6, 'AbsTol', 1e-9);
figure('Units', 'inches', 'Position', [0 0 10 10]); 
hold on;
colors = cool(length(phi_NPM_values));
legend_handles = gobjects(length(phi_NPM_values), 1);

% --- Plot Trajectories for Each 'phi_NPM' ---
for i = 1:length(phi_NPM_values)
    phi_NPM = phi_NPM_values(i);
    
    % Create a function handle for the ODE with the current 'phi_NPM'
    ode_function = @(t, y) replicator_dynamics_L(t, y, phi_NPM);
    
    % Solve the differential equation from the single starting point 'L'
    [~, Y] = ode45(ode_function, [0 2000], initial_condition, options);
    
    % Create the legend label using LaTeX interpreter
    legend_label = sprintf('$\\phi_{NPM} = %.3f$', phi_NPM);
    
    % Plot the phase trajectory and store its handle for the legend
    h = plot(Y(:,1), Y(:,2), '-', 'Color', colors(i,:), 'LineWidth', 2, 'DisplayName', legend_label);
    legend_handles(i) = h;
end

% --- Plot Annotations ---
% Plot the initial starting point 'L'
plot(initial_condition(1), initial_condition(2), 'ks', 'MarkerSize', 8, 'MarkerFaceColor', 'y');
text(initial_condition(1) - 0.04, initial_condition(2) - 0.04, '$L$', 'FontSize', 14, 'HorizontalAlignment', 'center','Interpreter','latex');

% Plot the pure-strategy equilibria with CORRECTED stability colors for the conflictual case
% E2 and E3 are stable (green), E1 and E4 are unstable (red)
plot(0, 0, 'ko', 'MarkerSize', 8, 'MarkerFaceColor', 'r', 'MarkerEdgeColor', 'r'); % Unstable
text(0-0.02, 0-0.02, '$E_1$', 'FontSize', 14, 'VerticalAlignment', 'top', 'HorizontalAlignment', 'right', 'Interpreter', 'latex');
plot(0, 1, 'ko', 'MarkerSize', 8, 'MarkerFaceColor', 'g', 'MarkerEdgeColor', 'g'); % Stable
text(0-0.02, 1+0.02, '$E_2$', 'FontSize', 14, 'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'right', 'Interpreter', 'latex');
plot(1, 0, 'ko', 'MarkerSize', 8, 'MarkerFaceColor', 'g', 'MarkerEdgeColor', 'g'); % Stable
text(1+0.02, 0-0.02, '$E_3$', 'FontSize', 14, 'VerticalAlignment', 'top', 'HorizontalAlignment', 'left', 'Interpreter', 'latex');
plot(1, 1, 'ko', 'MarkerSize', 8, 'MarkerFaceColor', 'r', 'MarkerEdgeColor', 'r'); % Unstable
text(1+0.02, 1+0.02, '$E_4$', 'FontSize', 14, 'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'left', 'Interpreter', 'latex');

% --- Final Figure Formatting ---
axis([0 1 0 1]); axis square;
xlabel('Founder''s Strategy, $p$ (Prob. Concentration)','Interpreter','latex');
ylabel('Successor''s Strategy, $q$ (Prob. Manager)','Interpreter','latex');
title('Trajectory Sensitivity to Agency Cost $\phi_{NPM}$ from Point $L$','Interpreter','latex');
grid on; box on; xticks(0:0.1:1); yticks(0:0.1:1);
set(gca, 'FontSize', 12);
hold off;

% --- Create, Format, and Add Padding to Legend ---
lgd = legend(legend_handles, 'Interpreter', 'latex', 'Location', 'northeastoutside');
lgd.Title.String = '\hspace{.5em}Trajectory for varying $\phi_{NPM}$\hspace{.5em}';
lgd.Title.Interpreter = 'latex';

% --- Save Figure ---
ax = gca;
ax.LineWidth = 1.2;
outputFileName = 'Sensitivity_L.pdf';
exportgraphics(ax, outputFileName, 'ContentType', 'vector');
disp(['Figure saved to: ' fullfile(pwd, outputFileName)]);

%%% --- LOCAL FUNCTIONS --- %%%
function dydt = replicator_dynamics_L(~, y, phi_NPM)
    % UPDATED fixed parameters for the "regime shift" scenario
    gamma = 0.5;
    omega = 0.01;
    phi_FS = 0.01;
    m = 0.1;
    n = -0.4;
    V = 100;
    W = 120;
    C = 6;
    K_m = 5;
    f_gamma = 45;
    
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