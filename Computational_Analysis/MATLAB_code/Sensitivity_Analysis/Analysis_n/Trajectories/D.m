%%% === SCRIPT SETUP === %%%
% --- Define System Parameters ---
% These are the new parameters for the 'D' simulation.
gamma = 0.5;
omega = 0.2;
phi_FS = 0.05;
phi_NPM = 0.1;
m = 0.05;
V = 1;
W = 1.2;
C = 0.25;
K_m = 0.05;
f_gamma = 0.355;

% --- Simulation Configuration ---
% Define the range of the parameter 'n' to be analyzed
n_values = -0.304:0.001:-0.296;
% Define the single initial condition for all trajectories
initial_condition = [0.2 0.7]; % This is the new starting point 'D'

% --- Plotting Setup ---
options = odeset('RelTol', 1e-6, 'AbsTol', 1e-9);
% Initialize figure and colormap
figure;
hold on;
colors = jet(length(n_values));

%%% --- Plot Trajectories for Each 'n' --- %%%
legend_handles = gobjects(length(n_values), 1); % Pre-allocate for legend handles
for i = 1:length(n_values)
    n = n_values(i);
    
    % Create a function handle for the ODE with the current 'n'
    ode_function = @(t, y) replicator_dynamics(t, y, n);
    
    % Solve the differential equation from the new initial condition
    [~, Y] = ode45(ode_function, [0 1000], initial_condition, options);
    
    % Create the legend label using LaTeX interpreter (formatted to 3 decimal places)
    legend_label = sprintf('$n = %.3f$', n);
    
    % Plot the phase trajectory and store its handle for the legend
    h = plot(Y(:,1), Y(:,2), 'o-', 'Color', colors(i,:), 'LineWidth', 2, 'DisplayName', legend_label);
    legend_handles(i) = h;
end

%%% --- Plot Annotations --- %%%
% Plot the initial starting point 'D' with a marker and a label
plot(initial_condition(1), initial_condition(2), 'ko', 'MarkerSize', 6, 'MarkerFaceColor', 'r');
text(initial_condition(1) + 0.01, initial_condition(2) - 0.03, 'D', 'FontSize', 13, ...
    'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'left');

% Plot the pure-strategy equilibria (E1-E4)
plot(0, 0, 'ko', 'MarkerSize', 5, 'MarkerFaceColor', 'g');
text(0+0.02, 0-0.01, 'E_1', 'FontSize', 12, 'VerticalAlignment', 'top');
plot(0, 1, 'ko', 'MarkerSize', 5, 'MarkerFaceColor', 'g');
text(0, 1, 'E_2', 'FontSize', 12, 'VerticalAlignment', 'bottom');
plot(1, 0, 'ko', 'MarkerSize', 5, 'MarkerFaceColor', 'g');
text(1+0.02, 0, 'E_3', 'FontSize', 12, 'HorizontalAlignment', 'left');
plot(1, 1, 'ko', 'MarkerSize', 5, 'MarkerFaceColor', 'g');
text(1, 1, 'E_4', 'FontSize', 12, 'VerticalAlignment', 'bottom');

%%% --- Final Figure Formatting --- %%%
axis([0 1 0 1]);
axis square;
xlabel('p');
ylabel('q');
grid on;
box on;
xticks(0:0.1:1);
yticks(0:0.1:1);
hold off;

%%% --- Create, Format, and Add Padding to Legend --- %%%
% Artificially lengthen the longest entry with invisible space
% to trick the layout engine into creating a wider box.
last_handle = legend_handles(end);
last_handle.DisplayName = [last_handle.DisplayName, '\hspace{1em}']; % Adjust space as needed

% Create the legend. It will now be wide enough.
lgd = legend(legend_handles, 'Interpreter', 'latex', 'AutoUpdate', 'off');

% Set the single-line, padded title, which will now fit perfectly.
lgd.Title.String = '\hspace{.5em}Varying Parameter $n$\hspace{.5em}';
lgd.Title.Interpreter = 'latex';
set(lgd, 'Location', 'northeastoutside');

% Adjust figure window size for better layout with legend
set(gcf, 'Position', [100 100 800 600]);

%%% --- Save Figure --- %%%
ax = gca;
ax.LineWidth = 1; % Set axis line width for publication quality

% Define the output filename
outputFileName = 'D.pdf'; % Changed filename to reflect the new plot

% Save the figure using exportgraphics, which tightly crops the output
exportgraphics(ax, outputFileName, 'ContentType', 'vector');
disp(['Figure saved to: ' fullfile(pwd, outputFileName)]);


%%% --- LOCAL FUNCTIONS --- %%%
function dydt = replicator_dynamics(~, y, n)
    % Define system parameters (kept here for function encapsulation)
    gamma = 0.5;
    omega = 0.2;
    phi_FS = 0.05;
    phi_NPM = 0.1;
    m = 0.05;
    V = 1;
    W = 1.2;
    C = 0.25;
    K_m = 0.05;
    f_gamma = 0.355;
   
    % Use p and q for better readability
    p = y(1);
    q = y(2);
    dydt = zeros(2,1);
    
    % Replicator dynamics equation for p (identical to original)
    dydt(1) = p*(1 - p)*( q*(phi_FS - n)*V - (1 - q)*n*W - f_gamma );
    
    % Replicator dynamics equation for q (identical to original)
    term1 = (p + (1 - p)*(gamma*(1 - omega) + (1 - gamma)*(phi_FS - n)))*V;
    term2 = ( (p + (1 - p)*gamma)*(1 - omega - phi_NPM + m) + (1 - p)*gamma*n )*W;
    dydt(2) = q*(1 - q)*( term1 - term2 - C + K_m );
end