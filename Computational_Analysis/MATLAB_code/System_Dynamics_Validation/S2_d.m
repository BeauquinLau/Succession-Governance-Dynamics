%%% === SCRIPT SETUP === %%%
% --- Define System Parameters ---
gamma = 0.5;
omega = 0.2;
phi_FS = 0.05;
phi_NPM = 0.1;
m = 0.05;
n = -0.3;
V = 1;
W = 1.2;
C = 0.25;
K_m = 0.05;
f_gamma = 0.355;

% --- Calculate the Saddle Point (E5) ---
% Calculate p*
numerator_p = C - K_m - (gamma*(1 - omega) + (1 - gamma)*(phi_FS - n))*V + ...
              gamma*(1 - omega - phi_NPM + m + n)*W;
            
denominator_p = (1 - (gamma * (1 - omega) + (1 - gamma) * (phi_FS - n)))*V + ...
               (gamma*(1 - omega - phi_NPM + m + n) - (1 - omega - phi_NPM + m))*W;
p_star = numerator_p / denominator_p;

% Calculate q*
numerator_q = n*W + f_gamma;
denominator_q = (phi_FS - n)*V + n*W;
q_star = numerator_q / denominator_q;

% --- Simulation and Plotting Setup ---
eps = 1e-5; % Perturbation size for manifold calculation
options = odeset('RelTol', 1e-6, 'AbsTol', 1e-9);

% Initialize figure
figure;
hold on;

% --- Calculate and Plot Manifolds ---
y1_val = p_star;
y2_val = q_star;

% Calculate Jacobian matrix numerically at the saddle point
J_E = calculate_jacobian(y1_val, y2_val);

% Compute eigenvalues and eigenvectors to find stable/unstable directions
[V, D] = eig(J_E);
lambda = diag(D);
[~, idx_pos] = max(real(lambda));
[~, idx_neg] = min(real(lambda));
v_unstable = V(:, idx_pos);
v_stable = V(:, idx_neg);

% Integrate along the unstable manifold
init_unstable_pos = [y1_val; y2_val] + eps * v_unstable;
init_unstable_neg = [y1_val; y2_val] - eps * v_unstable;
[~, Y_uni_pos] = ode45(@f1, [0 5000], init_unstable_pos, options);
[~, Y_uni_neg] = ode45(@f1, [0 5000], init_unstable_neg, options);

% Integrate along the stable manifold (backward in time)
init_stable_pos = [y1_val; y2_val] + eps * v_stable;
init_stable_neg = [y1_val; y2_val] - eps * v_stable;
[~, Y_sta_pos] = ode45(@f1, [0 -5000], init_stable_pos, options);
[~, Y_sta_neg] = ode45(@f1, [0 -5000], init_stable_neg, options);

% Plot the phase trajectories (manifolds)
plot(Y_uni_pos(:,1), Y_uni_pos(:,2), 'r--', 'LineWidth', 2);
plot(Y_uni_neg(:,1), Y_uni_neg(:,2), 'r--', 'LineWidth', 2);
plot(Y_sta_pos(:,1), Y_sta_pos(:,2), 'r--', 'LineWidth', 2);
plot(Y_sta_neg(:,1), Y_sta_neg(:,2), 'r--', 'LineWidth', 2);

% Plot the saddle point (E5)
plot(y1_val, y2_val, 'o', 'MarkerFaceColor', 'r', 'MarkerEdgeColor','k', 'MarkerSize', 8);
text(y1_val + 0.03, y2_val - 0.02, 'E_5', 'FontSize', 12);

% Plot the pure-strategy equilibria (E1-E4)
plot(0, 0, 'ko', 'MarkerSize', 5, 'MarkerFaceColor', 'g');
text(0+0.02, 0-0.01, 'E_1', 'FontSize', 12, 'VerticalAlignment', 'top');
plot(0, 1, 'ko', 'MarkerSize', 5, 'MarkerFaceColor', 'g');
text(0, 1, 'E_2', 'FontSize', 12, 'VerticalAlignment', 'bottom');
plot(1, 0, 'ko', 'MarkerSize', 5, 'MarkerFaceColor', 'g');
text(1+0.02, 0, 'E_3', 'FontSize', 12, 'HorizontalAlignment', 'left');
plot(1, 1, 'ko', 'MarkerSize', 5, 'MarkerFaceColor', 'g');
text(1, 1, 'E_4', 'FontSize', 12, 'VerticalAlignment', 'bottom');

% --- Figure Formatting ---
axis([0 1 0 1]);
axis square;
xlabel('p');
ylabel('q');
grid on;
xticks(0:0.1:1);
yticks(0:0.1:1);
box on;
hold off;

% --- Save Figure ---
ax = gca;
ax.LineWidth = 1; % Set axis line width for better visibility in publication

% Define the output filename
outputFileName = 'S2_d.pdf'; 

% Save the figure using exportgraphics for a high-quality vector image
exportgraphics(gca, outputFileName, 'ContentType', 'vector');
disp(['Figure saved to: ' fullfile(pwd, outputFileName)]);


%%% --- LOCAL FUNCTIONS --- %%%

function dy = f1(~, y)
    % Define system parameters
    gamma = 0.5;
    omega = 0.2;
    phi_FS = 0.05;
    phi_NPM = 0.1;
    m = 0.05;
    n = -0.3;
    V = 1;
    W = 1.2;
    C = 0.25;
    K_m = 0.05;
    f_gamma = 0.355;
    
    % Replicator dynamics equations
    dy = zeros(2, 1);
    dy(1) = y(1)*(1 - y(1))*(y(2)*(phi_FS - n)*V - (1 - y(2))*n*W - f_gamma);
    
    term1 = (y(1) + (1 - y(1))*(gamma*(1 - omega) + (1 - gamma)*(phi_FS - n)))*V;
    term2 = ((y(1) + (1 - y(1))*gamma)*(1 - omega - phi_NPM + m) + (1 - y(1))*gamma*n)*W;
    dy(2) = y(2)*(1 - y(2))*(term1 - term2 - C + K_m);
end

function J = calculate_jacobian(y1, y2)
    % System parameters (must be identical to those in f1)
    gamma = 0.5;
    omega = 0.2;
    phi_FS = 0.05;
    phi_NPM = 0.1;
    m = 0.05;
    n = -0.3;
    V = 1;
    W = 1.2;
    C = 0.25;
    K_m = 0.05;
    f_gamma = 0.355;

    % Analytically derived partial derivatives
    
    % Partial derivatives of the first ODE (for p)
    J11 = (1 - 2*y1) * (y2*(phi_FS - n)*V - (1 - y2)*n*W - f_gamma);
    J12 = y1*(1 - y1) * ((phi_FS - n)*V + n*W);
    
    % Partial derivatives of the second ODE (for q)
    d_term1_dy1 = (1 - (gamma*(1 - omega) + (1 - gamma)*(phi_FS - n)))*V;
    d_term2_dy1 = ((1 - gamma)*(1 - omega - phi_NPM + m) - gamma*n)*W;
    J21 = y2*(1 - y2) * (d_term1_dy1 - d_term2_dy1);
    
    term1 = (y1 + (1 - y1)*(gamma*(1 - omega) + (1 - gamma)*(phi_FS - n)))*V;
    term2 = ((y1 + (1 - y1)*gamma)*(1 - omega - phi_NPM + m) + (1 - y1)*gamma*n)*W;
    J22 = (1 - 2*y2) * (term1 - term2 - C + K_m);
    
    % Assemble the Jacobian matrix
    J = [J11, J12; J21, J22];
end