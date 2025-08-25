%%% === SCRIPT SETUP === %%%
% --- Define System Parameters ---
% These are the fixed parameters for the 'S2' simulation.
gamma = 0.5;
omega = 0.2;
phi_FS = 0.05;      % Renamed from phi_FS for consistency
phi_NPM = 0.1;
m = 0.05;
n = -0.3;
V = 1;            % Renamed from V for consistency
W = 1.2;    % Renamed from W for consistency
C = 0.25;         % Renamed from C for consistency
K_m = 0.05;
f_gamma = 0.355;

% --- Plotting Setup ---
options = odeset('RelTol', 1e-6, 'AbsTol', 1e-9);
% Initialize figure
figure;
hold on;

%%% --- Plot Trajectories from Various Initial Conditions --- %%%
% Define the ODE function handle once
ode_function = @(t, y) replicator_dynamics(t, y);

% Plot trajectories
[~,Y]=ode45(ode_function,[0 1000],[0.05 0.05]);
plot(Y(:,1),Y(:,2),'r-o','LineWidth', 2)
[~,Y]=ode45(ode_function,[0 1000],[0.1 0.05]);
plot(Y(:,1),Y(:,2),'r-o','LineWidth', 2)
[~,Y]=ode45(ode_function,[0 1000],[0.15 0.05]);
plot(Y(:,1),Y(:,2),'r-o','LineWidth', 2)
[~,Y]=ode45(ode_function,[0 1000],[0.2 0.05]);
plot(Y(:,1),Y(:,2),'r-o','LineWidth', 2)
[~,Y]=ode45(ode_function,[0 1000],[0.25 0.05]);
plot(Y(:,1),Y(:,2),'r-o','LineWidth', 2)
[~,Y]=ode45(ode_function,[0 1000],[0.3 0.05]);
plot(Y(:,1),Y(:,2),'r-o','LineWidth', 2)
[~,Y]=ode45(ode_function,[0 1000],[0.35 0.05]);
plot(Y(:,1),Y(:,2),'r-o','LineWidth', 2)
[~,Y]=ode45(ode_function,[0 1000],[0.4 0.05]);
plot(Y(:,1),Y(:,2),'c-o','LineWidth', 2)
[~,Y]=ode45(ode_function,[0 1000],[0.45 0.05]);
plot(Y(:,1),Y(:,2),'c-o','LineWidth', 2)
[~,Y]=ode45(ode_function,[0 1000],[0.45 0.2]);
plot(Y(:,1),Y(:,2),'c-o','LineWidth', 2)
[~,Y]=ode45(ode_function,[0 1000],[0.6 0.8]);
plot(Y(:,1),Y(:,2),'b-o','LineWidth', 2)
[~,Y]=ode45(ode_function,[0 1000],[0.6 0.85]);
plot(Y(:,1),Y(:,2),'m-o','LineWidth', 2)
[~,Y]=ode45(ode_function,[0 1000],[0.6 0.9]);
plot(Y(:,1),Y(:,2),'m-o','LineWidth', 2)
[~,Y]=ode45(ode_function,[0 1000],[0.6 0.95]);
plot(Y(:,1),Y(:,2),'m-o','LineWidth', 2)
[~,Y]=ode45(ode_function,[0 1000],[0.7 0.95]);
plot(Y(:,1),Y(:,2),'b-o','LineWidth', 2)
[~,Y]=ode45(ode_function,[0 1000],[0.75 0.95]);
plot(Y(:,1),Y(:,2),'b-o','LineWidth', 2)
[~,Y]=ode45(ode_function,[0 1000],[0.8 0.95]);
plot(Y(:,1),Y(:,2),'b-o','LineWidth', 2)
[~,Y]=ode45(ode_function,[0 1000],[0.85 0.95]);
plot(Y(:,1),Y(:,2),'b-o','LineWidth', 2)
[~,Y]=ode45(ode_function,[0 1000],[0.9 0.95]);
plot(Y(:,1),Y(:,2),'b-o','LineWidth', 2)
[~,Y]=ode45(ode_function,[0 1000],[0.95 0.95]);
plot(Y(:,1),Y(:,2),'b-o','LineWidth', 2)

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
outputFileName = 'S2.pdf'; % Changed filename to reflect the new plot

% Save the figure using exportgraphics, which tightly crops the output
exportgraphics(ax, outputFileName, 'ContentType', 'vector');
disp(['Figure saved to: ' fullfile(pwd, outputFileName)]);


%%% --- LOCAL FUNCTIONS --- %%%
function dydt = replicator_dynamics(~, y)
    % Define system parameters (now defined globally in the main script)
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
   
    % Use p and q for better readability
    p = y(1);
    q = y(2);
    dydt = zeros(2,1);
    
    % Replicator dynamics equation for p
    dydt(1) = p*(1 - p)*( q*(phi_FS - n)*V - (1 - q)*n*W - f_gamma );
    
    % Replicator dynamics equation for q
    term1 = (p + (1 - p)*(gamma*(1 - omega) + (1 - gamma)*(phi_FS - n)))*V;
    term2 = ( (p + (1 - p)*gamma)*(1 - omega - phi_NPM + m) + (1 - p)*gamma*n )*W;
    dydt(2) = q*(1 - q)*( term1 - term2 - C + K_m );
end