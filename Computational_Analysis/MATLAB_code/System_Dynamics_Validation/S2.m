%%% === SCRIPT SETUP FOR S2 (Regime Shift) === %%%
% --- Define System Parameters ---
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
f_gamma = 45;

% --- Plotting Setup ---
options = odeset('RelTol', 1e-6, 'AbsTol', 1e-9);
figure('Units', 'inches', 'Position', [0 0 10 10]); 
hold on;

%%% --- Plot Trajectories from Various Initial Conditions --- %%%
ode_function = @(t, y) replicator_dynamics(t, y);

% Select four representative colors from the 'cool' colormap
% The sequence progresses from cyan to magenta.
cool_colors = cool(4); 

% Maintain the original grouping of initial conditions, but remove hardcoded colors
initial_conditions = {
    % Group 1, converges to E3
    {[0.02, 0.15], [0.05, 0.12], [0.05, 0.15], [0.1 0.2], [0.15 0.2]};
    % Group 2, converges to E3
    {[0.05, 0.05], [0.05, 0.1], [0.1, 0.1], [0.15 0.05], [0.2, 0.2]};
    % Group 3, converges to E2
    {[0.95, 0.95], [0.9, 0.85], [0.9, 0.82], [0.95, 0.85], [0.98, 0.9]};
    % Group 4, converges to E2
    {[0.98, 0.85], [0.95, 0.8], [0.9, 0.75], [0.85, 0.7], [0.8, 0.6]}
};

% Loop through the four groups and plot trajectories using the new cool colors
for i = 1:length(initial_conditions)
    conditions = initial_conditions{i};
    % Assign a specific color from the colormap to the current group
    current_color = cool_colors(i, :); 
    
    for j = 1:length(conditions)
        ic = conditions{j};
        [~, Y] = ode45(ode_function, [0 1500], ic, options);
        % Explicitly specify the color in the plot command while keeping the '-o' style
        plot(Y(:,1), Y(:,2), '-o', 'Color', current_color, 'LineWidth', 1.5);
    end
end

%%% --- Plot Annotations (E2 and E3 are now stable) --- %%%
plot(0, 0, 'ko', 'MarkerSize', 8, 'MarkerFaceColor', 'r', 'MarkerEdgeColor', 'r'); % Unstable
text(0-0.02, 0-0.02, '$E_1$', 'FontSize', 14, 'VerticalAlignment', 'top', 'HorizontalAlignment', 'right', 'Interpreter', 'latex');
plot(0, 1, 'ko', 'MarkerSize', 8, 'MarkerFaceColor', 'g'); % Stable
text(0-0.02, 1+0.02, '$E_2$', 'FontSize', 14, 'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'right', 'Interpreter', 'latex');
plot(1, 0, 'ko', 'MarkerSize', 8, 'MarkerFaceColor', 'g'); % Stable
text(1+0.02, 0-0.02, '$E_3$', 'FontSize', 14, 'VerticalAlignment', 'top', 'HorizontalAlignment', 'left', 'Interpreter', 'latex');
plot(1, 1, 'ko', 'MarkerSize', 8, 'MarkerFaceColor', 'r', 'MarkerEdgeColor', 'r'); % Unstable
text(1+0.02, 1+0.02, '$E_4$', 'FontSize', 14, 'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'left', 'Interpreter', 'latex');

% Calculate and plot the saddle point E5 (optional but good for completeness)
numerator_p = C - K_m - gamma*(1-omega-phi_FS+n)*V + gamma*(1-omega-phi_NPM+m+n)*W;
denominator_p = (1 - gamma*(1-omega-phi_FS+n))*V - ((1-omega-phi_NPM+m) - gamma*(1-omega-phi_NPM+m+n))*W;
p_star = numerator_p / denominator_p;
numerator_q = n * W + f_gamma;
denominator_q = (phi_FS - n) * V + n * W;
q_star = numerator_q / denominator_q;
plot(p_star, q_star, 'x', 'MarkerSize', 10, 'MarkerEdgeColor','k', 'LineWidth', 2); 
text(p_star+0.02, q_star+0.01, '$E_5$', 'FontSize', 14, 'VerticalAlignment', 'top', 'HorizontalAlignment', 'center', 'Interpreter', 'latex');

%%% --- Final Figure Formatting --- %%%
axis([0 1 0 1]); axis square;
xlabel('Founder''s Strategy, $p$ (Prob. Concentration)','Interpreter','latex');
ylabel('Successor''s Strategy, $q$ (Prob. Manager)','Interpreter','latex');
title('Phase Portrait under Conflict: Regime Shift to $E_2$ and $E_3$','Interpreter','latex');
grid on; box on; xticks(0:0.1:1); yticks(0:0.1:1);
set(gca, 'FontSize', 12);
hold off;

%%% --- Save Figure --- %%%
ax = gca; ax.LineWidth = 1.2;
outputFileName = 'S2.pdf';
exportgraphics(ax, outputFileName, 'ContentType', 'vector');
disp(['Figure saved to: ' fullfile(pwd, outputFileName)]);

%%% --- LOCAL FUNCTIONS --- %%%
function dydt = replicator_dynamics(~, y)
    % Parameters for Regime Shift simulation
    gamma=0.5; omega=0.01; phi_FS=0.01; phi_NPM=0.25; m=0.1;
    n=-0.4; V=100; W=120; C=6; K_m=5; f_gamma=45;
   
    p = y(1); q = y(2);
    dydt = zeros(2,1);
    
    dydt(1) = p*(1 - p)*( q*(phi_FS - n)*V - (1 - q)*n*W - f_gamma );
    
    % Dynamic equation for q
    X_FSM = ( p + (1 - p)*gamma*(1 - omega - phi_FS + n) )*V - C;
    X_FSS = ( p*(1 - omega - phi_NPM + m) + (1 - p)*gamma*(1 - omega - phi_NPM + m + n) )*W - K_m;
    dydt(2) = q*(1 - q)*( X_FSM - X_FSS );
end