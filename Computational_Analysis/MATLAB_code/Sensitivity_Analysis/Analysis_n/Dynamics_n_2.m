%%% === SCRIPT SETUP FOR CONFLICTUAL CASE SENSITIVITY ANALYSIS === %%%
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
n_values = -0.43:0.005:-0.38; 
p_star = zeros(size(n_values));
q_star = zeros(size(n_values));

% --- Calculate Saddle Point Coordinates for each 'n' ---
for i = 1:length(n_values)
    n = n_values(i);
    
    % The formulas for p* and q* are the corrected ones and remain valid
    numerator_q = n*W + f_gamma;
    denominator_q = (phi_FS - n)*V + n*W;
    q_star(i) = numerator_q / denominator_q;
    
    % p* formula
    numerator_p = C - K_m - gamma*(1-omega-phi_FS+n)*V + gamma*(1-omega-phi_NPM+m+n)*W;
    denominator_p = (1 - gamma*(1-omega-phi_FS+n))*V - ((1-omega-phi_NPM+m) - gamma*(1-omega-phi_NPM+m+n))*W;
    p_star(i) = numerator_p / denominator_p;
end

% --- Plotting Setup ---
eps = 1e-5;
options = odeset('RelTol', 1e-6, 'AbsTol', 1e-9);
figure('Units', 'inches', 'Position', [0 0 10 10]); 
hold on;
colors = cool(length(n_values)); 
legend_handles = gobjects(length(n_values), 1);

% --- Calculate and Plot Manifolds for Each 'n' ---
for i = 1:length(n_values)
    n = n_values(i);
    p_val = p_star(i);
    q_val = q_star(i);
    
    J_E = calculate_jacobian(p_val, q_val, n);
    [V_eig, D_eig] = eig(J_E);
    lambda = diag(D_eig);
    [~, idx_pos] = max(real(lambda));
    [~, idx_neg] = min(real(lambda));
    v_unstable = V_eig(:, idx_pos);
    v_stable = V_eig(:, idx_neg);
    
    % Integrate along UNSTABLE manifold (Forward in time)
    init_unstable_pos = [p_val; q_val] + eps * v_unstable;
    init_unstable_neg = [p_val; q_val] - eps * v_unstable;
    [~, Y_uni_pos] = ode45(@(t, y) replicator_dynamics(t, y, n), [0 2000], init_unstable_pos, options);
    [~, Y_uni_neg] = ode45(@(t, y) replicator_dynamics(t, y, n), [0 2000], init_unstable_neg, options);
    
    % Integrate along STABLE manifold (Backward in time)
    init_stable_pos = [p_val; q_val] + eps * v_stable;
    init_stable_neg = [p_val; q_val] - eps * v_stable;
    [~, Y_sta_pos] = ode45(@(t, y) replicator_dynamics(t, y, n), [0 -2000], init_stable_pos, options);
    [~, Y_sta_neg] = ode45(@(t, y) replicator_dynamics(t, y, n), [0 -2000], init_stable_neg, options);
    
    % Plot the UNSTABLE manifold arms
    h = plot(Y_uni_pos(:,1), Y_uni_pos(:,2), '--', 'Color', colors(i,:), 'LineWidth', 2, 'DisplayName', sprintf('$n=%.2f$', n));
    legend_handles(i) = h;
    plot(Y_uni_neg(:,1), Y_uni_neg(:,2), '--', 'Color', colors(i,:), 'LineWidth', 2, 'HandleVisibility','off');
    
    % Plot the STABLE manifold arms
    plot(Y_sta_pos(:,1), Y_sta_pos(:,2), '--', 'Color', colors(i,:), 'LineWidth', 2, 'HandleVisibility','off');
    plot(Y_sta_neg(:,1), Y_sta_neg(:,2), '--', 'Color', colors(i,:), 'LineWidth', 2, 'HandleVisibility','off');
    
    % Plot the saddle point marker
    plot(p_val, q_val, 'o', 'MarkerFaceColor', colors(i,:), 'MarkerEdgeColor','k', 'HandleVisibility','off');
end

% --- Plot Fixed Points and Annotations ---
plot(0, 0, 'ko', 'MarkerSize', 8, 'MarkerFaceColor', 'r');
text(0-0.02, 0-0.02, '$E_1$', 'FontSize', 14, 'VerticalAlignment', 'top', 'HorizontalAlignment', 'right', 'Interpreter', 'latex');
plot(0, 1, 'ko', 'MarkerSize', 8, 'MarkerFaceColor', 'g', 'MarkerEdgeColor', 'r');
text(0-0.02, 1+0.02, '$E_2$', 'FontSize', 14, 'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'right', 'Interpreter', 'latex');
plot(1, 0, 'ko', 'MarkerSize', 8, 'MarkerFaceColor', 'g', 'MarkerEdgeColor', 'r');
text(1+0.02, 0-0.02, '$E_3$', 'FontSize', 14, 'VerticalAlignment', 'top', 'HorizontalAlignment', 'left', 'Interpreter', 'latex');
plot(1, 1, 'ko', 'MarkerSize', 8, 'MarkerFaceColor', 'r');
text(1+0.02, 1+0.02, '$E_4$', 'FontSize', 14, 'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'left', 'Interpreter', 'latex');

% Plot representative initial points for analysis (D, E, F)
plot(0.2, 0.8, 'ks', 'MarkerSize', 8, 'MarkerFaceColor', 'y');
text(0.2, 0.8 - 0.02, '$D$', 'FontSize', 13, 'VerticalAlignment', 'top', 'HorizontalAlignment', 'center', 'Interpreter', 'latex');
plot(0.5, 0.4, 'ks', 'MarkerSize', 8, 'MarkerFaceColor', 'y');
text(0.5, 0.4 - 0.02, '$E$', 'FontSize', 13, 'VerticalAlignment', 'top', 'HorizontalAlignment', 'center', 'Interpreter', 'latex');
plot(0.8, 0.2, 'ks', 'MarkerSize', 8, 'MarkerFaceColor', 'y');
text(0.8, 0.2 - 0.02, '$F$', 'FontSize', 13, 'VerticalAlignment', 'top', 'HorizontalAlignment', 'center', 'Interpreter', 'latex');

% --- Final Figure Formatting ---
axis([0 1 0 1]); axis square;
xlabel('Founder''s Strategy, $p$ (Prob. Concentration)','Interpreter','latex');
ylabel('Successor''s Strategy, $q$ (Prob. Manager)','Interpreter','latex');
title('Sensitivity Analysis of Conflict Intensity $n$','Interpreter','latex');
grid on; box on; xticks(0:0.1:1); yticks(0:0.1:1);
set(gca, 'FontSize', 12, 'Layer', 'top');
lgd = legend(legend_handles, 'Interpreter', 'latex', 'Location', 'northeastoutside');
lgd.Title.String = '\hspace{.5em}Saddle point $E_5$ for varying $n$\hspace{.5em}';
lgd.Title.Interpreter = 'latex';
hold off;

% --- Save Figure ---
ax = gca; ax.LineWidth = 1.2;
outputFileName = 'Dynamics_n_2.pdf';
exportgraphics(ax, outputFileName, 'ContentType', 'vector');
disp(['Figure saved to: ' fullfile(pwd, outputFileName)]);


%%% --- LOCAL FUNCTIONS --- %%%
function dydt = replicator_dynamics(~, y, n)
    % Parameters for Regime Shift simulation
    gamma=0.5; omega=0.01; phi_FS=0.01; phi_NPM=0.25; m=0.1;
    V=100; W=120; C=6; K_m=5; f_gamma=45;
    
    p = y(1); q = y(2);
    dydt = zeros(2,1);
    
    dydt(1) = p*(1 - p)*( q*(phi_FS - n)*V - (1 - q)*n*W - f_gamma );
    
    X_FSM = ( p + (1 - p)*gamma*(1 - omega - phi_FS + n) )*V - C;
    X_FSS = ( p*(1 - omega - phi_NPM + m) + (1 - p)*gamma*(1 - omega - phi_NPM + m + n) )*W - K_m;
    dydt(2) = q*(1 - q)*( X_FSM - X_FSS );
end

function J = calculate_jacobian(p, q, n)
    gamma=0.5; omega=0.01; phi_FS=0.01; phi_NPM=0.25; m=0.1;
    V=100; W=120; C=6; K_m=5; f_gamma=45;
    
    A_q = q*(phi_FS - n)*V - (1 - q)*n*W - f_gamma;
    J11 = (1 - 2*p) * A_q;
    J12 = p*(1 - p) * ((phi_FS - n)*V + n*W);
    
    X_FSM = ( p + (1 - p)*gamma*(1 - omega - phi_FS + n) )*V - C;
    X_FSS = ( p*(1 - omega - phi_NPM + m) + (1 - p)*gamma*(1 - omega - phi_NPM + m + n) )*W - K_m;
    B_p = X_FSM - X_FSS;
    d_XFSM_dp = (1 - gamma*(1 - omega - phi_FS + n))*V;
    d_XFSS_dp = ((1 - omega - phi_NPM + m) - gamma*(1 - omega - phi_NPM + m + n))*W;
    B_p_prime = d_XFSM_dp - d_XFSS_dp;
    J21 = q*(1 - q) * B_p_prime;
    J22 = (1 - 2*q) * B_p;
    
    J = [J11, J12; J21, J22];
end