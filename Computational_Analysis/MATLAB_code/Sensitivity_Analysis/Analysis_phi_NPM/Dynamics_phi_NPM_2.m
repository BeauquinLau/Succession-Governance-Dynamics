%%% === SCRIPT SETUP FOR PHI_NPM SENSITIVITY ANALYSIS (CONFLICTUAL CASE) === %%%
% --- Define System Parameters ---
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
% Define the range of 'phi_NPM' to be analyzed
% UPDATED based on the newly calculated valid range for the regime shift: (0.2233, 0.265)
phi_NPM_values = 0.23:0.002:0.26;
p_star = zeros(size(phi_NPM_values));
q_star = zeros(size(phi_NPM_values));

% --- Calculate Saddle Point Coordinates for each 'phi_NPM' ---
for i = 1:length(phi_NPM_values)
    phi_NPM = phi_NPM_values(i);
    
    % q* formula (is constant as it does not depend on phi_NPM)
    numerator_q = n*W + f_gamma;
    denominator_q = (phi_FS - n)*V + n*W;
    q_star(i) = numerator_q / denominator_q;
    
    % p* formula (Corrected based on manuscript equations)
    numerator_p = C - K_m - gamma*(1-omega-phi_FS+n)*V + gamma*(1-omega-phi_NPM+m+n)*W;
    denominator_p = (1 - gamma*(1-omega-phi_FS+n))*V - ((1-omega-phi_NPM+m) - gamma*(1-omega-phi_NPM+m+n))*W;
    p_star(i) = numerator_p / denominator_p;
end

% --- Plotting Setup ---
eps = 1e-5;
options = odeset('RelTol', 1e-6, 'AbsTol', 1e-9);
figure('Units', 'inches', 'Position', [0 0 10 10]); 
hold on;
colors = cool(length(phi_NPM_values)); 
legend_handles = gobjects(length(phi_NPM_values), 1);

% --- Calculate and Plot Manifolds for Each 'phi_NPM' ---
for i = 1:length(phi_NPM_values)
    phi_NPM = phi_NPM_values(i);
    p_val = p_star(i);
    q_val = q_star(i);
    
    % Get Jacobian and eigenvectors
    J_E = calculate_jacobian_conflict(p_val, q_val, phi_NPM);
    [V_eig, D_eig] = eig(J_E);
    lambda = diag(D_eig);
    [~, idx_pos] = max(real(lambda));
    [~, idx_neg] = min(real(lambda));
    v_unstable = V_eig(:, idx_pos);
    v_stable = V_eig(:, idx_neg);
    
    % Create a function handle for the ODE with the current phi_NPM
    ode_function = @(t, y) replicator_dynamics_conflict(t, y, phi_NPM);
    
    % Integrate along manifolds
    [~, Y_uni_pos] = ode45(ode_function, [0 2000], [p_val; q_val] + eps * v_unstable, options);
    [~, Y_uni_neg] = ode45(ode_function, [0 2000], [p_val; q_val] - eps * v_unstable, options);
    [~, Y_sta_pos] = ode45(ode_function, [0 -2000], [p_val; q_val] + eps * v_stable, options);
    [~, Y_sta_neg] = ode45(ode_function, [0 -2000], [p_val; q_val] - eps * v_stable, options);
    
    % Plot the manifold arms
    h = plot(Y_uni_pos(:,1), Y_uni_pos(:,2), '--', 'Color', colors(i,:), 'LineWidth', 2, 'DisplayName', sprintf('$\\phi_{NPM}=%.3f$', phi_NPM));
    legend_handles(i) = h;
    plot(Y_uni_neg(:,1), Y_uni_neg(:,2), '--', 'Color', colors(i,:), 'LineWidth', 2, 'HandleVisibility','off');
    plot(Y_sta_pos(:,1), Y_sta_pos(:,2), '--', 'Color', colors(i,:), 'LineWidth', 2, 'HandleVisibility','off');
    plot(Y_sta_neg(:,1), Y_sta_neg(:,2), '--', 'Color', colors(i,:), 'LineWidth', 2, 'HandleVisibility','off');
    
    % Plot the saddle point marker
    plot(p_val, q_val, 'o', 'MarkerFaceColor', colors(i,:), 'MarkerEdgeColor','k', 'HandleVisibility','off');
end

% --- Plot Fixed Points and Annotations (CONFLICTUAL CASE) ---
% E2 and E3 are stable (green), E1 and E4 are unstable (red)
plot(0, 0, 'ko', 'MarkerSize', 8, 'MarkerFaceColor', 'r'); % Unstable
text(0-0.02, 0-0.02, '$E_1$', 'FontSize', 14, 'VerticalAlignment', 'top', 'HorizontalAlignment', 'right', 'Interpreter', 'latex');
plot(0, 1, 'ko', 'MarkerSize', 8, 'MarkerFaceColor', 'g'); % Stable
text(0-0.02, 1+0.02, '$E_2$', 'FontSize', 14, 'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'right', 'Interpreter', 'latex');
plot(1, 0, 'ko', 'MarkerSize', 8, 'MarkerFaceColor', 'g'); % Stable
text(1+0.02, 0-0.02, '$E_3$', 'FontSize', 14, 'VerticalAlignment', 'top', 'HorizontalAlignment', 'left', 'Interpreter', 'latex');
plot(1, 1, 'ko', 'MarkerSize', 8, 'MarkerFaceColor', 'r'); % Unstable
text(1+0.02, 1+0.02, '$E_4$', 'FontSize', 14, 'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'left', 'Interpreter', 'latex');

% Plot representative initial points J, K, L
plot(0.2, 0.7, 'ks', 'MarkerSize', 8, 'MarkerFaceColor', 'y');
text(0.2, 0.7 - 0.04, '$J$', 'FontSize', 14, 'HorizontalAlignment', 'center', 'Interpreter', 'latex');
plot(0.6, 0.5, 'ks', 'MarkerSize', 8, 'MarkerFaceColor', 'y');
text(0.6, 0.5 - 0.04, '$K$', 'FontSize', 14, 'HorizontalAlignment', 'center', 'Interpreter', 'latex');
plot(0.9, 0.4, 'ks', 'MarkerSize', 8, 'MarkerFaceColor', 'y');
text(0.9, 0.4 + 0.04, '$L$', 'FontSize', 14, 'HorizontalAlignment', 'center', 'Interpreter', 'latex');

% --- Final Figure Formatting ---
axis([0 1 0 1]); axis square;
xlabel('Founder''s Strategy, $p$ (Prob. Concentration)','Interpreter','latex');
ylabel('Successor''s Strategy, $q$ (Prob. Manager)','Interpreter','latex');
title('Sensitivity Analysis of Agency Cost $\phi_{NPM}$ (Conflictual Case)','Interpreter','latex');
grid on; box on; xticks(0:0.1:1); yticks(0:0.1:1);
set(gca, 'FontSize', 12, 'Layer', 'top');
lgd = legend(legend_handles, 'Interpreter', 'latex', 'Location', 'northeastoutside');
lgd.Title.String = '\hspace{.5em}Saddle point $E_5$ for varying $\phi_{NPM}$\hspace{.5em}';
lgd.Title.Interpreter = 'latex';
hold off;

% --- Save Figure ---
ax = gca; ax.LineWidth = 1.2;
outputFileName = 'Dynamics_phi_NPM_2.pdf';
exportgraphics(ax, outputFileName, 'ContentType', 'vector');
disp(['Figure saved to: ' fullfile(pwd, outputFileName)]);

%%% --- LOCAL FUNCTIONS (Conflictual Case) --- %%%
function dydt = replicator_dynamics_conflict(~, y, phi_NPM)
    % System Parameters for Conflictual Case
    gamma = 0.5; omega = 0.01; phi_FS = 0.01; m = 0.1; n = -0.4;
    V = 100; W = 120; C = 6; K_m = 5; f_gamma = 45;
    
    p = y(1); q = y(2);
    dydt = zeros(2,1);
    
    % dp/dt equation
    dydt(1) = p*(1 - p)*( q*(phi_FS - n)*V - (1 - q)*n*W - f_gamma );
    
    % dq/dt equation
    X_FSM = ( p + (1 - p)*gamma*(1 - omega - phi_FS + n) )*V - C;
    X_FSS = ( p*(1 - omega - phi_NPM + m) + (1 - p)*gamma*(1 - omega - phi_NPM + m + n) )*W - K_m;
    dydt(2) = q*(1 - q)*( X_FSM - X_FSS );
end

function J = calculate_jacobian_conflict(p, q, phi_NPM)
    % System Parameters for Conflictual Case
    gamma = 0.5; omega = 0.01; phi_FS = 0.01; m = 0.1; n = -0.4;
    V = 100; W = 120; C = 6; K_m = 5; f_gamma = 45;
    
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