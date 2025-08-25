%%% === SCRIPT SETUP === %%%
% --- Define System Parameters ---
gamma = 0.5;        
omega = 0.28;
phi_FS = 0.2;
phi_NPM = 0.25;
m = 0.1;
n = 0.1;
V = 100;          
W = 90;     
C = 50;
K_m = 20;
f_gamma = 4;

% --- Calculate the Saddle Point (E5) ---
% Calculate p*
numerator_p = C - K_m - (gamma * (1 - omega) + (1 - gamma) * (phi_FS - n)) * V + gamma * (1 - omega - phi_NPM + m + n) * W;
denominator_p = (1 - (gamma * (1 - omega) + (1 - gamma) * (phi_FS - n))) * V + (gamma * (1 - omega - phi_NPM + m + n) - (1 - omega - phi_NPM + m)) * W;
p_star = numerator_p / denominator_p;

% Calculate q*
numerator_q = n * W + f_gamma;
denominator_q = (phi_FS - n) * V + n * W;
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
J_E = calculate_jacobian(y1_val, y2_val, n);

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
[~, Y_uni_pos] = ode45(@(t, y) f1(t, y, n), [0 5000], init_unstable_pos, options);
[~, Y_uni_neg] = ode45(@(t, y) f1(t, y, n), [0 5000], init_unstable_neg, options);

% Integrate along the stable manifold (backward in time)
init_stable_pos = [y1_val; y2_val] + eps * v_stable;
init_stable_neg = [y1_val; y2_val] - eps * v_stable;
[~, Y_sta_pos] = ode45(@(t, y) f1(t, y, n), [0 -5000], init_stable_pos, options);
[~, Y_sta_neg] = ode45(@(t, y) f1(t, y, n), [0 -5000], init_stable_neg, options);

% Plot the four manifold arms
plot(Y_uni_pos(:,1), Y_uni_pos(:,2), '--', 'Color', 'r', 'LineWidth', 2); 
plot(Y_uni_neg(:,1), Y_uni_neg(:,2), '--', 'Color', 'r', 'LineWidth', 2, 'HandleVisibility','off'); 
plot(Y_sta_pos(:,1), Y_sta_pos(:,2), '--', 'Color', 'r', 'LineWidth', 2, 'HandleVisibility','off');
plot(Y_sta_neg(:,1), Y_sta_neg(:,2), '--', 'Color', 'r', 'LineWidth', 2, 'HandleVisibility','off');

% Plot the saddle point (E5)
plot(y1_val, y2_val, 'o', 'MarkerFaceColor', 'r', 'MarkerEdgeColor','k', 'HandleVisibility','off'); 
text(y1_val - 0.05, y2_val - 0.01, 'E_5', 'FontSize', 12, 'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'left');

% Plot the four pure-strategy equilibria (E1-E4)
plot(0, 0, 'ko', 'MarkerSize', 5, 'MarkerFaceColor', 'g');
text(0+0.02, 0-0.01, 'E_1', 'FontSize', 12, 'VerticalAlignment', 'top');
plot(0, 1, 'ko', 'MarkerSize', 5, 'MarkerFaceColor', 'g');
text(0, 1, 'E_2', 'FontSize', 12, 'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'left');
plot(1, 0, 'ko', 'MarkerSize', 5, 'MarkerFaceColor', 'g');
text(1+0.02, 0, 'E_3', 'FontSize', 12, 'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'left');
plot(1, 1, 'ko', 'MarkerSize', 5, 'MarkerFaceColor', 'g');
text(1, 1, 'E_4', 'FontSize', 12, 'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'left');

% --- Figure Formatting ---
axis([0 1 0 1])
axis('square')
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
outputFileName = 'S1_d.pdf';

% Save the figure using exportgraphics for a high-quality vector image
exportgraphics(gca, outputFileName, 'ContentType', 'vector');
disp(['Figure saved to: ' fullfile(pwd, outputFileName)]);


%%% --- LOCAL FUNCTIONS --- %%%

function dydt = f1(~, y, n)
    % Define system parameters
    gamma = 0.5;        
    omega = 0.28;
    phi_FS = 0.2;
    phi_NPM = 0.25;
    m = 0.1;
    V = 100;          
    W = 90;     
    C = 50;
    K_m = 20;
    f_gamma = 4;
    
    y1 = y(1);
    y2 = y(2);
    
    dydt = zeros(2,1);

    % Replicator dynamics for p
    dydt(1) = y1*(1 - y1)*( y2*(phi_FS - n)*V - (1 - y2)*n*W - f_gamma ); % <-- 直接赋值
    
    % Replicator dynamics for q
    term1 = (y1 + (1 - y1)*(gamma*(1 - omega) + (1 - gamma)*(phi_FS - n)))*V;
    term2 = ( (y1 + (1 - y1)*gamma)*(1 - omega - phi_NPM + m) + (1 - y1)*gamma*n )*W;
    dydt(2) = y2*(1 - y2)*( term1 - C - term2 + K_m );
end

function J = calculate_jacobian(y1, y2, n)
    % System parameters (must be identical to those in f1)
    gamma = 0.5;        
    omega = 0.28;
    phi_FS = 0.2;
    phi_NPM = 0.25;
    m = 0.1;
    V = 100;          
    W = 90;     
    C = 50;
    K_m = 20;
    f_gamma = 4;

    % Analytically derived partial derivatives
    
    % Partial derivatives of the first ODE (dy1dt)
    J11 = (1 - 2*y1) * (y2*(phi_FS - n)*V - (1 - y2)*n*W - f_gamma);
    J12 = y1*(1 - y1) * ((phi_FS - n)*V + n*W);
    
    % Partial derivatives of the second ODE (dy2dt)
    d_term1_dy1 = (1 - (gamma*(1 - omega) + (1 - gamma)*(phi_FS - n)))*V;
    d_term2_dy1 = ((1 - gamma)*(1 - omega - phi_NPM + m) - gamma*n)*W;
    J21 = y2*(1 - y2) * (d_term1_dy1 - d_term2_dy1);
    
    term1 = (y1 + (1 - y1)*(gamma*(1 - omega) + (1 - gamma)*(phi_FS - n)))*V;
    term2 = ( (y1 + (1 - y1)*gamma)*(1 - omega - phi_NPM + m) + (1 - y1)*gamma*n )*W;
    J22 = (1 - 2*y2) * (term1 - C - term2 + K_m);
    
    % Assemble the Jacobian matrix
    J = [J11, J12; J21, J22];
end