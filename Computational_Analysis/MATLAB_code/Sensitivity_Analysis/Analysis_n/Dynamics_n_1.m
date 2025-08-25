%%% === SCRIPT SETUP === %%%
% --- Define System Parameters ---
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

% --- Simulation Configuration ---
% Define the range of the parameter 'n' to be analyzed
n_values = 0.02:0.01:0.11; 

% Pre-allocate arrays to store the results for the saddle point coordinates
p_star = zeros(size(n_values));
q_star = zeros(size(n_values));

% --- Calculate Saddle Point Coordinates for each 'n' ---
for i = 1:length(n_values)
    n = n_values(i);
    
    % Calculate q*
    numerator_q = n*W + f_gamma;
    denominator_q = (phi_FS - n)*V + n*W;
    q_star(i) = numerator_q / denominator_q;
    
    % Calculate p*
    term1 = C - K_m;
    term2 = -(gamma*(1-omega) + (1-gamma)*(phi_FS - n)) * V;
    term3 = gamma*(1 - omega - phi_NPM + m + n)*W;
    numerator_p = term1 + term2 + term3;
    
    denom_term1 = (1 - (gamma * (1 - omega) + (1 - gamma) * (phi_FS - n))) * V;
    denom_term2 = (gamma*(1 - omega - phi_NPM + m + n) - (1 - omega - phi_NPM + m)) * W;
    denominator_p = denom_term1 + denom_term2;
    
    p_star(i) = numerator_p / denominator_p;
end

% --- Plotting Setup ---
eps = 1e-5; % Perturbation size for numerical stability
options = odeset('RelTol', 1e-6, 'AbsTol', 1e-9);

% Initialize figure and colormap
figure;
hold on;
colors = jet(length(n_values));
legend_handles = gobjects(length(n_values), 1); % Pre-allocate for legend handles

% --- Calculate and Plot Manifolds for Each 'n' ---
for i = 1:length(n_values)
    n = n_values(i);
    p_val = p_star(i);
    q_val = q_star(i);
    
    % Calculate the Jacobian matrix numerically at the saddle point
    J_E = calculate_jacobian(p_val, q_val, n);
    
    % Eigenvalue analysis to find stable and unstable directions
    [V, D] = eig(J_E);
    lambda = diag(D);
    [~, idx_pos] = max(real(lambda));
    [~, idx_neg] = min(real(lambda));
    v_unstable = V(:, idx_pos);
    v_stable = V(:, idx_neg);
    
    % Integrate along the unstable manifold
    init_unstable_pos = [p_val; q_val] + eps * v_unstable;
    init_unstable_neg = [p_val; q_val] - eps * v_unstable;
    [~, Y_uni_pos] = ode45(@(t, y) replicator_dynamics(t, y, n), [0 1000], init_unstable_pos, options);
    [~, Y_uni_neg] = ode45(@(t, y) replicator_dynamics(t, y, n), [0 1000], init_unstable_neg, options);
    
    % Integrate along the stable manifold (backward in time)
    init_stable_pos = [p_val; q_val] + eps * v_stable;
    init_stable_neg = [p_val; q_val] - eps * v_stable;
    [~, Y_sta_pos] = ode45(@(t, y) replicator_dynamics(t, y, n), [0 -1000], init_stable_pos, options);
    [~, Y_sta_neg] = ode45(@(t, y) replicator_dynamics(t, y, n), [0 -1000], init_stable_neg, options);
    
    % Plot the primary trajectory and store its handle for the legend
    h = plot(Y_uni_pos(:,1), Y_uni_pos(:,2), '--', 'Color', colors(i,:), 'LineWidth', 2, 'DisplayName', sprintf('$n=%.2f$', n));
    legend_handles(i) = h;
    
    % Plot the other three manifold arms without legend entries
    plot(Y_uni_neg(:,1), Y_uni_neg(:,2), '--', 'Color', colors(i,:), 'LineWidth', 2, 'HandleVisibility','off');
    plot(Y_sta_pos(:,1), Y_sta_pos(:,2), '--', 'Color', colors(i,:), 'LineWidth', 2, 'HandleVisibility','off');
    plot(Y_sta_neg(:,1), Y_sta_neg(:,2), '--', 'Color', colors(i,:), 'LineWidth', 2, 'HandleVisibility','off');
    
    % Plot the saddle point marker
    plot(p_val, q_val, 'o', 'MarkerFaceColor', colors(i,:), 'MarkerEdgeColor','k', 'HandleVisibility','off');
end

% --- Plot Fixed Points and Annotations ---
% Plot the pure-strategy equilibria (E1-E4)
plot(0, 0, 'ko', 'MarkerSize', 5, 'MarkerFaceColor', 'g');
text(0+0.02, 0-0.01, 'E_1', 'FontSize', 12, 'VerticalAlignment', 'top');
plot(0, 1, 'ko', 'MarkerSize', 5, 'MarkerFaceColor', 'g');
text(0, 1, 'E_2', 'FontSize', 12, 'VerticalAlignment', 'bottom');
plot(1, 0, 'ko', 'MarkerSize', 5, 'MarkerFaceColor', 'g');
text(1+0.02, 0, 'E_3', 'FontSize', 12, 'HorizontalAlignment', 'left');
plot(1, 1, 'ko', 'MarkerSize', 5, 'MarkerFaceColor', 'g');
text(1, 1, 'E_4', 'FontSize', 12, 'VerticalAlignment', 'bottom');

% Plot representative initial points for analysis
plot(0.2, 0.4, 'ko', 'MarkerSize', 6, 'MarkerFaceColor', 'r');
text(0.2 - 0.03, 0.4 - 0.05, 'A', 'FontSize', 13, 'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'left');
plot(0.6, 0.4, 'ko', 'MarkerSize', 6, 'MarkerFaceColor', 'r');
text(0.602, 0.4-0.05, 'B', 'FontSize', 13, 'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'left');
plot(0.9, 0.5, 'ko', 'MarkerSize', 6, 'MarkerFaceColor', 'r');
text(0.902, 0.502, 'C', 'FontSize', 13, 'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'left');

% --- Final Figure Formatting ---
axis([0 1 0 1]);
axis square;
xlabel('p');
ylabel('q');
grid on;
box on;
xticks(0:0.1:1);
yticks(0:0.1:1);
hold off;

% --- Create, Format, and Add Padding to Legend ---
% Artificially lengthen the longest entry with invisible space
% to trick the layout engine into creating a wider box.
last_handle = legend_handles(end);
last_handle.DisplayName = [last_handle.DisplayName, '\hspace{1em}']; % Adjust the em value as needed

% Create the legend. It will now be wide enough.
lgd = legend(legend_handles, 'Interpreter', 'latex', 'AutoUpdate', 'off');

% Set the single-line, padded title, which will now fit perfectly.
lgd.Title.String = '\hspace{0.5em}Varying Parameter $n$\hspace{0.5em}';
lgd.Title.Interpreter = 'latex';
set(lgd, 'Location', 'northeastoutside');

% Adjust figure window size for better layout with legend
set(gcf, 'Position', [100 100 800 600]);

% --- Save Figure ---
ax = gca;
ax.LineWidth = 1; % Set axis line width for publication quality

% Define the output filename
outputFileName = 'Dynamics_n_1.pdf';

% Save the figure using exportgraphics, which tightly crops the output
exportgraphics(ax, outputFileName, 'ContentType', 'vector');
disp(['Figure saved to: ' fullfile(pwd, outputFileName)]);

%%% --- LOCAL FUNCTIONS --- %%%

function dydt = replicator_dynamics(~, y, n)
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
    
    p = y(1);
    q = y(2);
    
    dydt = zeros(2,1);
    
    % Replicator dynamics equation for p
    dydt(1) = p*(1 - p)*( q*(phi_FS - n)*V - (1 - q)*n*W - f_gamma );
    
    % Replicator dynamics equation for q
    term1 = (p + (1 - p)*(gamma*(1 - omega) + (1 - gamma)*(phi_FS - n)))*V;
    term2 = ( (p + (1 - p)*gamma)*(1 - omega - phi_NPM + m) + (1 - p)*gamma*n )*W;
    dydt(2) = q*(1 - q)*( term1 - C - term2 + K_m );
end

function J = calculate_jacobian(p, q, n)
    % System parameters (must be identical to those in the ODE function)
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
    
    % Partial derivatives of the first ODE (dpdt)
    J11 = (1 - 2*p) * (q*(phi_FS - n)*V - (1 - q)*n*W - f_gamma);
    J12 = p*(1 - p) * ((phi_FS - n)*V + n*W);
    
    % Partial derivatives of the second ODE (dqdt)
    d_term1_dp = (1 - (gamma*(1 - omega) + (1 - gamma)*(phi_FS - n)))*V;
    d_term2_dp = ((1 - gamma)*(1 - omega - phi_NPM + m) - gamma*n)*W;
    J21 = q*(1 - q) * (d_term1_dp - d_term2_dp);
    
    term1 = (p + (1 - p)*(gamma*(1 - omega) + (1 - gamma)*(phi_FS - n)))*V;
    term2 = ( (p + (1 - p)*gamma)*(1 - omega - phi_NPM + m) + (1 - p)*gamma*n )*W;
    J22 = (1 - 2*q) * (term1 - C - term2 + K_m);
    
    % Assemble the Jacobian matrix
    J = [J11, J12; J21, J22];
end