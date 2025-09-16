%%% === SCRIPT SETUP FOR F_GAMMA SENSITIVITY ANALYSIS (CONFLICTUAL CASE) === %%%
% --- Define System Parameters (from Table "params_extended") ---
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
% Define the range of 'f_gamma' to be analyzed.
% The theoretical range for a valid interior saddle point in this regime is (41, 48).
f_gamma_values = 42:0.5:47;
p_star = zeros(size(f_gamma_values));
q_star = zeros(size(f_gamma_values));

% --- Calculate Saddle Point Coordinates for each 'f_gamma' ---
for i = 1:length(f_gamma_values)
    f_gamma = f_gamma_values(i);
    
    % q* formula (depends on f_gamma)
    numerator_q = n*W + f_gamma;
    denominator_q = (phi_FS - n)*V + n*W;
    q_star(i) = numerator_q / denominator_q;
    
    % p* formula (is constant as it does not depend on f_gamma)
    numerator_p = C - K_m - gamma*(1-omega-phi_FS+n)*V + gamma*(1-omega-phi_NPM+m+n)*W;
    denominator_p = (1 - gamma*(1-omega-phi_FS+n))*V - ((1-omega-phi_NPM+m) - gamma*(1-omega-phi_NPM+m+n))*W;
    p_star(i) = numerator_p / denominator_p;
end

% --- Plotting Setup ---
eps = 1e-5;
options = odeset('RelTol', 1e-6, 'AbsTol', 1e-9);
figure('Units', 'inches', 'Position', [0 0 10 10]); 
hold on;
colors = cool(length(f_gamma_values)); 
legend_handles = gobjects(length(f_gamma_values), 1);

% --- Calculate and Plot Manifolds for Each 'f_gamma' ---
for i = 1:length(f_gamma_values)
    f_gamma = f_gamma_values(i);
    p_val = p_star(i);
    q_val = q_star(i);
    
    % Create a function handle for the ODE with the current f_gamma
    ode_function = @(t, y) replicator_dynamics_fgamma_conflict(t, y, f_gamma);
    
    % Get Jacobian and eigenvectors
    J_E = calculate_jacobian_fgamma_conflict(p_val, q_val, f_gamma);
    [V_eig, D_eig] = eig(J_E);
    lambda = diag(D_eig);
    [~, idx_pos] = max(real(lambda));
    [~, idx_neg] = min(real(lambda));
    v_unstable = V_eig(:, idx_pos);
    v_stable = V_eig(:, idx_neg);
    
    % Integrate along manifolds
    [~, Y_uni_pos] = ode45(ode_function, [0 2000], [p_val; q_val] + eps * v_unstable, options);
    [~, Y_uni_neg] = ode45(ode_function, [0 2000], [p_val; q_val] - eps * v_unstable, options);
    [~, Y_sta_pos] = ode45(ode_function, [0 -2000], [p_val; q_val] + eps * v_stable, options);
    [~, Y_sta_neg] = ode45(ode_function, [0 -2000], [p_val; q_val] - eps * v_stable, options);
    
    % Plot the manifold arms
    h = plot(Y_uni_pos(:,1), Y_uni_pos(:,2), '--', 'Color', colors(i,:), 'LineWidth', 2, 'DisplayName', sprintf('$f_{\\gamma}=%.1f$', f_gamma));
    legend_handles(i) = h;
    plot(Y_uni_neg(:,1), Y_uni_neg(:,2), '--', 'Color', colors(i,:), 'LineWidth', 2, 'HandleVisibility','off');
    plot(Y_sta_pos(:,1), Y_sta_pos(:,2), '--', 'Color', colors(i,:), 'LineWidth', 2, 'HandleVisibility','off');
    plot(Y_sta_neg(:,1), Y_sta_neg(:,2), '--', 'Color', colors(i,:), 'LineWidth', 2, 'HandleVisibility','off');
    
    % Plot the saddle point marker
    plot(p_val, q_val, 'o', 'MarkerFaceColor', colors(i,:), 'MarkerEdgeColor','k', 'HandleVisibility','off');
end

% --- Plot Fixed Points and Annotations (for Conflictual Case) ---
% E2 and E3 are stable (green), E1 and E4 are unstable (red)
plot(0, 0, 'ko', 'MarkerSize', 8, 'MarkerFaceColor', 'r', 'MarkerEdgeColor', 'r');
text(0-0.02, 0-0.02, '$E_1$', 'FontSize', 14, 'VerticalAlignment', 'top', 'HorizontalAlignment', 'right', 'Interpreter', 'latex');
plot(0, 1, 'ko', 'MarkerSize', 8, 'MarkerFaceColor', 'g', 'MarkerEdgeColor', 'g');
text(0-0.02, 1+0.02, '$E_2$', 'FontSize', 14, 'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'right', 'Interpreter', 'latex');
plot(1, 0, 'ko', 'MarkerSize', 8, 'MarkerFaceColor', 'g', 'MarkerEdgeColor', 'g');
text(1+0.02, 0-0.02, '$E_3$', 'FontSize', 14, 'VerticalAlignment', 'top', 'HorizontalAlignment', 'left', 'Interpreter', 'latex');
plot(1, 1, 'ko', 'MarkerSize', 8, 'MarkerFaceColor', 'r', 'MarkerEdgeColor', 'r');
text(1+0.02, 1+0.02, '$E_4$', 'FontSize', 14, 'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'left', 'Interpreter', 'latex');

% Plot representative initial points (optional for drift plots, but kept for context)
plot(0.2, 0.8, 'ks', 'MarkerSize', 8, 'MarkerFaceColor', 'y');
text(0.2, 0.8 - 0.04, '$P$', 'FontSize', 14, 'HorizontalAlignment', 'center', 'Interpreter', 'latex');
plot(0.5, 0.5, 'ks', 'MarkerSize', 8, 'MarkerFaceColor', 'y');
text(0.5, 0.5 - 0.04, '$Q$', 'FontSize', 14, 'HorizontalAlignment', 'center', 'Interpreter', 'latex');
plot(0.8, 0.2, 'ks', 'MarkerSize', 8, 'MarkerFaceColor', 'y');
text(0.8, 0.2 - 0.04, '$R$', 'FontSize', 14, 'HorizontalAlignment', 'center', 'Interpreter', 'latex');

% --- Final Figure Formatting ---
axis([0 1 0 1]); axis square;
xlabel('Founder''s Strategy, $p$ (Prob. Concentration)','Interpreter','latex');
ylabel('Successor''s Strategy, $q$ (Prob. Manager)','Interpreter','latex');
title('Sensitivity Analysis of Affective Benefit $f_{\gamma}$ (Conflictual Case)','Interpreter','latex');
grid on; box on; xticks(0:0.1:1); yticks(0:0.1:1);
set(gca, 'FontSize', 12, 'Layer', 'top');
hold off;

% --- Create, Format, and Add Padding to Legend ---
lgd = legend(legend_handles, 'Interpreter', 'latex', 'Location', 'northeastoutside');
lgd.Title.String = '\hspace{.5em}Saddle point $E_5$ for varying $f_{\gamma}$\hspace{.5em}';
lgd.Title.Interpreter = 'latex';

% --- Save Figure ---
ax = gca; ax.LineWidth = 1.2;
outputFileName = 'Dynamics_f_gamma_2.pdf';
exportgraphics(ax, outputFileName, 'ContentType', 'vector');
disp(['Figure saved to: ' fullfile(pwd, outputFileName)]);

%%% --- LOCAL FUNCTIONS (Conflictual Case) --- %%%
function dydt = replicator_dynamics_fgamma_conflict(~, y, f_gamma)
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
    
    p = y(1); q = y(2);
    dydt = zeros(2,1);
    
    % dp/dt equation
    dydt(1) = p*(1 - p)*( q*(phi_FS - n)*V - (1 - q)*n*W - f_gamma );
    
    % dq/dt equation
    X_FSM = ( p + (1 - p)*gamma*(1 - omega - phi_FS + n) )*V - C;
    X_FSS = ( p*(1 - omega - phi_NPM + m) + (1 - p)*gamma*(1 - omega - phi_NPM + m + n) )*W - K_m;
    dydt(2) = q*(1 - q)*( X_FSM - X_FSS );
end

function J = calculate_jacobian_fgamma_conflict(p, q, f_gamma)
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
    
    % Partial derivatives for Jacobian matrix (from manuscript Table 3)
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