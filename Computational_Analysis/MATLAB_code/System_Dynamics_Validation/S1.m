%%% === SCRIPT SETUP === %%%
% --- Define System Parameters ---
% These are the fixed parameters for the 'S1' simulation.
gamma = 0.5;        
omega = 0.28;
phi_FS = 0.2;       % Renamed from phi_FS for consistency
phi_NPM = 0.25;
m = 0.1;
n = 0.1;
V = 100;          % Renamed from V for consistency
W = 90;     % Renamed from W for consistency
C = 50;           % Renamed from C for consistency
K_m = 20;
f_gamma = 4;

% --- Plotting Setup ---
options = odeset('RelTol', 1e-6, 'AbsTol', 1e-9);
% Initialize figure
figure;
hold on;

%%% --- Plot Trajectories from Various Initial Conditions --- %%%
% Define the ODE function handle once
ode_function = @(t, y) replicator_dynamics(t, y);

% Group initial conditions by color for cleaner code
initial_conditions = {
    % Red trajectories
    {[0.1 0.95], [0.2 0.9], [0.2 0.95], [0.3 0.9]}, 'r-o';
    % Cyan trajectories
    {[0.3 0.95], [0.4 0.82], [0.4 0.85], [0.4 0.95], [0.6 0.95], [0.65 0.95]}, 'c-o';
    % Magenta trajectories
    {[0.7 0.4], [0.8 0.05], [0.8 0.15], [0.8 0.2], [0.9 0.05]}, 'm-o';
    % Blue trajectories
    {[0.8 0.4], [0.85 0.2], [0.85 0.3], [0.9 0.2], [0.95 0.1]}, 'b-o'
};

% Loop through the groups and plot trajectories
for i = 1:size(initial_conditions, 1)
    conditions = initial_conditions{i, 1};
    style = initial_conditions{i, 2};
    for j = 1:length(conditions)
        ic = conditions{j};
        [~, Y] = ode45(ode_function, [0 1000], ic, options);
        plot(Y(:,1), Y(:,2), style, 'LineWidth', 2);
    end
end

%%% --- Plot Annotations --- %%%
% Plot the four pure-strategy equilibria
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

%%% --- Save Figure --- %%%
ax = gca;
ax.LineWidth = 1; % Set axis line width for publication quality

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
    C = 50;
    K_m = 20;
    f_gamma = 4;
   
    % Use p and q for better readability
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