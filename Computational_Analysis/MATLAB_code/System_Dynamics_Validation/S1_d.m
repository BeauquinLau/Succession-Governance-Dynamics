%%% === SCRIPT SETUP === %%%
% --- Define System Parameters ---
% These are the fixed parameters for the simulation
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

% --- Calculate the Saddle Point (E5) using formulas ---
numerator_p = C - K_m - gamma*(1-omega-phi_FS+n)*V + gamma*(1-omega-phi_NPM+m+n)*W;
denominator_p = (1 - gamma*(1-omega-phi_FS+n))*V - ((1-omega-phi_NPM+m) - gamma*(1-omega-phi_NPM+m+n))*W;
p_star = numerator_p / denominator_p;
numerator_q = n * W + f_gamma;
denominator_q = (phi_FS - n) * V + n * W;
q_star = numerator_q / denominator_q;

% Display the new E5 coordinates
fprintf('New saddle point E5 is at (p*, q*) = (%.4f, %.4f)\n', p_star, q_star);

% --- Simulation and Plotting Setup ---
eps = 1e-5; 
options = odeset('RelTol', 1e-6, 'AbsTol', 1e-9);
figure('Units', 'inches', 'Position', [0 0 10 10]); 
hold on;

% --- Calculate and Plot Manifolds (Separatrix) ---
y1_val = p_star;
y2_val = q_star;
J_E = calculate_jacobian(y1_val, y2_val);
[V_eig, D_eig] = eig(J_E);
lambda = diag(D_eig);
[~, idx_pos] = max(real(lambda));
[~, idx_neg] = min(real(lambda));
v_unstable = V_eig(:, idx_pos);
v_stable = V_eig(:, idx_neg);

% Integrate along UNSTABLE manifold (forward in time)
init_unstable_pos = [y1_val; y2_val] + eps * v_unstable;
init_unstable_neg = [y1_val; y2_val] - eps * v_unstable;
[~, Y_uni_pos] = ode45(@(t, y) replicator_dynamics(t, y), [0 5000], init_unstable_pos, options);
[~, Y_uni_neg] = ode45(@(t, y) replicator_dynamics(t, y), [0 5000], init_unstable_neg, options);

% Integrate along STABLE manifold (backward in time)
init_stable_pos = [y1_val; y2_val] + eps * v_stable;
init_stable_neg = [y1_val; y2_val] - eps * v_stable;
[~, Y_sta_pos] = ode45(@(t, y) replicator_dynamics(t, y), [0 -5000], init_stable_pos, options);
[~, Y_sta_neg] = ode45(@(t, y) replicator_dynamics(t, y), [0 -5000], init_stable_neg, options);

% Plot the four manifold arms to form the complete separatrix
plot(Y_uni_pos(:,1), Y_uni_pos(:,2), '--', 'Color', [0.85 0.33 0.1], 'LineWidth', 2); 
plot(Y_uni_neg(:,1), Y_uni_neg(:,2), '--', 'Color', [0.85 0.33 0.1], 'LineWidth', 2, 'HandleVisibility','off'); 
plot(Y_sta_pos(:,1), Y_sta_pos(:,2), '--', 'Color', [0.85 0.33 0.1], 'LineWidth', 2, 'HandleVisibility','off');
plot(Y_sta_neg(:,1), Y_sta_neg(:,2), '--', 'Color', [0.85 0.33 0.1], 'LineWidth', 2, 'HandleVisibility','off');

% --- Plot Annotations ---
plot(0, 0, 'ko', 'MarkerSize', 8, 'MarkerFaceColor', 'g');
text(0-0.02, 0-0.02, '$E_1$', 'FontSize', 14, 'VerticalAlignment', 'top', 'HorizontalAlignment', 'right', 'Interpreter', 'latex');
plot(0, 1, 'ko', 'MarkerSize', 8, 'MarkerFaceColor', 'r', 'MarkerEdgeColor', 'r');
text(0-0.02, 1+0.02, '$E_2$', 'FontSize', 14, 'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'right', 'Interpreter', 'latex');
plot(1, 0, 'ko', 'MarkerSize', 8, 'MarkerFaceColor', 'r', 'MarkerEdgeColor', 'r');
text(1+0.02, 0-0.02, '$E_3$', 'FontSize', 14, 'VerticalAlignment', 'top', 'HorizontalAlignment', 'left', 'Interpreter', 'latex');
plot(1, 1, 'ko', 'MarkerSize', 8, 'MarkerFaceColor', 'g');
text(1+0.02, 1+0.02, '$E_4$', 'FontSize', 14, 'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'left', 'Interpreter', 'latex');
plot(y1_val, y2_val, 'x', 'MarkerSize', 10, 'MarkerEdgeColor','k', 'LineWidth', 2); 
text(y1_val, y2_val-0.03, '$E_5$', 'FontSize', 14, 'VerticalAlignment', 'top', 'HorizontalAlignment', 'center', 'Interpreter', 'latex');

% --- Final Figure Formatting ---
axis([0 1 0 1]);
axis square;
xlabel('Founder''s Strategy, $p$ (Prob. Concentration)','Interpreter','latex');
ylabel('Successor''s Strategy, $q$ (Prob. Manager)','Interpreter','latex');
title('Phase Portrait of Governance Dynamics','Interpreter','latex');
grid on;
box on;
xticks(0:0.1:1);
yticks(0:0.1:1);
set(gca, 'FontSize', 12, 'Layer', 'top');
hold off;

% --- Save Figure ---
ax = gca;
ax.LineWidth = 1.2;
outputFileName = 'S1_d.pdf';
exportgraphics(ax, outputFileName, 'ContentType', 'vector');
disp(['Figure saved to: ' fullfile(pwd, outputFileName)]);


%%% --- LOCAL FUNCTIONS --- %%%
function dydt = replicator_dynamics(~, y)
    gamma = 0.5; omega = 0.28; phi_FS = 0.2; phi_NPM = 0.25; m = 0.1; n = 0.1;
    V = 100; W = 90; C = 45; K_m = 20; f_gamma = 1;
    p = y(1); q = y(2);
    dydt = zeros(2,1);
    dydt(1) = p*(1 - p)*( q*(phi_FS - n)*V - (1 - q)*n*W - f_gamma );
    X_FSM = ( p + (1 - p)*gamma*(1 - omega - phi_FS + n) )*V - C;
    X_FSS = ( p*(1 - omega - phi_NPM + m) + (1 - p)*gamma*(1 - omega - phi_NPM + m + n) )*W - K_m;
    dydt(2) = q*(1 - q)*( X_FSM - X_FSS );
end

function J = calculate_jacobian(p, q)
    gamma = 0.5; omega = 0.28; phi_FS = 0.2; phi_NPM = 0.25; m = 0.1; n = 0.1;
    V = 100; W = 90; C = 45; K_m = 20; f_gamma = 1;
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