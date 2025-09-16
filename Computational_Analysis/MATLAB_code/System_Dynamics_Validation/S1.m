%%% === SCRIPT SETUP === %%%
% --- Define System Parameters ---
% These are the fixed parameters for the 'S1' simulation.
gamma = 0.5;        
omega = 0.28;
phi_FS = 0.2;
phi_NPM = 0.25;
m = 0.1;
n = 0.1;
V = 100;
W = 90;
C = 45;
K_m = 20;
f_gamma = 1;

% --- Plotting Setup ---
options = odeset('RelTol', 1e-6, 'AbsTol', 1e-9);
% Initialize figure
figure('Units', 'inches', 'Position', [0 0 10 10]); 
hold on;

%%% --- Plot Trajectories from Various Initial Conditions --- %%%
% Define the ODE function handle once
ode_function = @(t, y) replicator_dynamics(t, y);

% Select four representative colors from the 'cool' colormap
% cool(4) generates a 4x3 matrix, where each row is an RGB color vector
% The sequence progresses from cyan to magenta
cool_colors = cool(4); 

% Maintain the original grouping of initial conditions, but remove hardcoded colors
initial_conditions = {
    % Group 1
    {[0.05, 0.95], [0.1 0.95], [0.15 0.92], [0.2 0.9], [0.3 0.8]};
    % Group 2
    {[0.2,0.95], [0.3 0.85], [0.3 0.95], [0.4 0.7], [0.4 0.95], [0.6 0.95]};
    % Group 3
    {[0.7 0.15], [0.75 0.15], [0.75 0.02], [0.8 0.05], [0.9 0.02]};
    % Group 4
    {[0.7 0.23], [0.8 0.15], [0.9 0.05], [0.9 0.2], [0.95 0.1]}
};

% Loop through the four groups and plot trajectories using the new cool colors
for i = 1:length(initial_conditions)
    conditions = initial_conditions{i};
    % Assign a specific color from the colormap to the current group
    current_color = cool_colors(i, :); 
    
    for j = 1:length(conditions)
        ic = conditions{j};
        [~, Y] = ode45(ode_function, [0 1000], ic, options);
        % Explicitly specify the color in the plot command while keeping the '-o' style
        plot(Y(:,1), Y(:,2), '-o', 'Color', current_color, 'LineWidth', 2);
    end
end

%%% --- Plot Annotations --- %%%
% Plot the four pure-strategy equilibria
plot(0, 0, 'ko', 'MarkerSize', 8, 'MarkerFaceColor', 'g');
text(0-0.02, 0-0.02, '$E_1$', 'FontSize', 14, 'VerticalAlignment', 'top', 'HorizontalAlignment', 'right', 'Interpreter', 'latex');
plot(0, 1, 'ko', 'MarkerSize', 8, 'MarkerFaceColor', 'r', 'MarkerEdgeColor', 'r');
text(0-0.02, 1+0.02, '$E_2$', 'FontSize', 14, 'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'right', 'Interpreter', 'latex');
plot(1, 0, 'ko', 'MarkerSize', 8, 'MarkerFaceColor', 'r', 'MarkerEdgeColor', 'r');
text(1+0.02, 0-0.02, '$E_3$', 'FontSize', 14, 'VerticalAlignment', 'top', 'HorizontalAlignment', 'left', 'Interpreter', 'latex');
plot(1, 1, 'ko', 'MarkerSize', 8, 'MarkerFaceColor', 'g');
text(1+0.02, 1+0.02, '$E_4$', 'FontSize', 14, 'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'left', 'Interpreter', 'latex');

% Calculate the saddle point E5 using the formulas
numerator_p = C - K_m - gamma*(1-omega-phi_FS+n)*V + gamma*(1-omega-phi_NPM+m+n)*W;
denominator_p = (1 - gamma*(1-omega-phi_FS+n))*V - ((1-omega-phi_NPM+m) - gamma*(1-omega-phi_NPM+m+n))*W;
p_star = numerator_p / denominator_p;
numerator_q = n * W + f_gamma;
denominator_q = (phi_FS - n) * V + n * W;
q_star = numerator_q / denominator_q;
% Plot the saddle point E5
plot(p_star, q_star, 'kx', 'MarkerSize', 10, 'LineWidth', 2); 
text(p_star-0.01, q_star-0.02, '$E_5$', 'FontSize', 14, 'VerticalAlignment', 'top', 'HorizontalAlignment', 'center', 'Interpreter', 'latex');

%%% --- Final Figure Formatting --- %%%
axis([0 1 0 1]);
axis square;
xlabel('Founder''s Strategy, $p$ (Prob. Concentration)','Interpreter','latex');
ylabel('Successor''s Strategy, $q$ (Prob. Manager)','Interpreter','latex');
title('Phase Portrait of Governance Dynamics (Baseline Case)','Interpreter','latex');
grid on;
box on;
xticks(0:0.1:1);
yticks(0:0.1:1);
set(gca, 'FontSize', 12);
hold off;

%%% --- Save Figure --- %%%
ax = gca;
ax.LineWidth = 1; 

% Define the output filename
outputFileName = 'S1.pdf';

% Save the figure using exportgraphics for a high-quality vector image
exportgraphics(ax, outputFileName, 'ContentType', 'vector');
disp(['Figure saved to: ' fullfile(pwd, outputFileName)]);


%%% --- LOCAL FUNCTIONS --- %%%
function dydt = replicator_dynamics(~, y)
    % Define system parameters
    gamma = 0.5;        
    omega = 0.28;
    phi_FS = 0.2;
    phi_NPM = 0.25;
    m = 0.1;
    n = 0.1;
    V = 100;          
    W = 90;     
    C = 45;
    K_m = 20;
    f_gamma = 1;

    % Use p and q for better readability
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