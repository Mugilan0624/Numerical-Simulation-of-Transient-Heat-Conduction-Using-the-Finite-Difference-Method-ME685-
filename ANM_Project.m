% ANM Project

% Domain: 3mm x 3mm, Grid: 120 x 120 nodes
% Laser power: 1500 W and the Beam radius: 50 um.

clear; clc; close all;

%% Material and Physical Properties %% Titanium alloy Grade 5 (Ti-6Al-4V)
rho = 4430;              % Density (kg/m^3)
cp  = 560;               % Specific heat (J/kg.K)
k   = 6.7;               % Thermal conductivity (W/m.K)
alpha = k / (rho*cp);    % Thermal diffusivity (m^2/s)
T_ambient = 300;         % Ambient temperature (K)

%% Convective Heat Loss Parameter
h = 50;                  % Convective heat transfer coefficient (W/m^2.K)

%% Domain Setup
Nx = 120;                % Number of nodes in x-direction
Ny = 120;                % Number of nodes in y-direction
L = 0.003;               % Domain length (meters)
dx = L / Nx;             % Grid spacing in x
dy = dx;                 % Uniform grid spacing in y

x = linspace(0, L, Nx);  % x-coordinate vector
y = linspace(0, L, Ny);  % y-coordinate vector
[X, Y] = meshgrid(x, y); % 2D grid for spatial domain

%% Time Setup
t_end = 1;                             % Simulation end time in seconds
dt = 0.25 * dx^2 / alpha;                % Time step (CFL condition for stability)
Nt = floor(t_end / dt);                  % Total number of time steps

%% Snapshot Times Setup
% Snapshots for 0, 0.25, 0.5, 0.75, and 1.0 seconds
snapshot_times = [0, 0.25, 0.5, 0.75, 1.0];  
snapshot_indices = round(snapshot_times / dt);  % Convert times to indices
T_snapshots = cell(length(snapshot_times), 1);  % Cell array to store temperature fields
snapshot_counter = 1;

%% Laser (Gaussian heat source) Parameters
P = 1500;              % Laser power (W)
r0 = 0.00005;          % Laser beam radius (m) [50 um]
x0 = L / 2;            % Laser center x-coordinate
y0 = L / 2;            % Laser center y-coordinate

%% Initial Condition
T = T_ambient * ones(Nx, Ny);  % Initialize temperature field to ambient temp

% Save the initial snapshot (t = 0 s)
T_snapshots{snapshot_counter} = T;
snapshot_counter = snapshot_counter + 1;

%% MP4 Video Setup (Optional)
v = VideoWriter('Temperature_Evolution.mp4', 'MPEG-4'); % Create video file
v.FrameRate = 30;
v.Quality = 100;
open(v);  % Open video file for writing

fig = figure('Position', [200, 0, 1000, 800]); % Figure for simulation video frames

%% Main Time-Stepping Loop
for n = 1:Nt
    Q = gaussianHeatSource(X, Y, x0, y0, P, r0);  % Compute laser heat source
    T = computeTemp(T, Q, alpha, dt, dx, rho, cp, T_ambient, h);  % Update temperature field

    % Save snapshot if the current time index matches one of the snapshot indices
    if snapshot_counter <= length(snapshot_indices) && n == snapshot_indices(snapshot_counter)
        T_snapshots{snapshot_counter} = T;
        snapshot_counter = snapshot_counter + 1;
    end

    % Render and record a video frame every 30 time steps (optional)
    if mod(n, 30) == 0
        imagesc(x, y, T');                  
        set(gca, 'YDir', 'normal');         
        colorbar;
        title(['Temperature distribution at t = ', num2str(n*dt, '%.3f'), ' s']);
        xlabel('x (m)');
        ylabel('y (m)');
        caxis([300 500]);
        axis square;
        drawnow;
        frame = getframe(fig);
        writeVideo(v, frame);
    end
end

% Ensure the final time (t = 1 s) snapshot is saved if not already captured in loop
if snapshot_counter <= length(snapshot_times)
    T_snapshots{snapshot_counter} = T;
end

close(v);  % Finalize video file
disp('Simulation complete. Video (if enabled) and snapshots are created.');

%% Combined Mid-line (x = L/2) Temperature Line Plot for All Snapshot Times
mid_x_index = round(Nx/2);  % Index corresponding to x = L/2

figure('Position', [250, 100, 800, 400]);
hold on;  % Hold on to overlay multiple line plots

for i = 1:length(snapshot_times)
    T_cross = T_snapshots{i}(mid_x_index, :);  % Temperature along the mid-line
    plot(y, T_cross, 'LineWidth', 2, 'DisplayName', sprintf('t = %.2f s', snapshot_times(i)));
end

title('Temperature Distribution Along x = L/2');
xlabel('y (m)');
ylabel('Temperature (K)');
grid on;
legend('show', 'Location', 'best');
axis([min(y) max(y) 300 500]);  % Fixed temperature range
hold off;
saveas(gcf, 'Combined_TempLine_xL2.png');

%% Enhanced Separate Discrete Contour Plots for Each Selected Time
nLevels = 15;                           
levels = linspace(300, 500, nLevels);  % Define levels from 300 to 500 K

for i = 1:length(snapshot_times)
    figure('Position',[200, 0, 1000, 800]);
    
    % Create a filled contour plot with discrete levels
    [Cs_filled, hContourF] = contourf(x, y, T_snapshots{i}', levels, 'LineStyle','none');
    colormap(parula(nLevels-1));   % Use a discrete colormap with nLevels-1 bands
    caxis([levels(1), levels(end)]);   % Fix the color axis
    
    % Add a colorbar with ticks matching the discrete levels
    cb = colorbar;
    cb.Ticks = levels;
    cb.TickLabels = round(levels, 1);
    cb.Label.String = 'Temperature (K)';
    cb.Label.FontSize = 12;
    
    hold on;
    % Overlay contour lines with the same discrete levels and add labels
    [C, hContourL] = contour(x, y, T_snapshots{i}', levels, 'LineColor','k', 'LineWidth',1);
    clabel(C, hContourL, 'FontSize', 8, 'Color', 'k');  % Label contour lines
    
    title(sprintf('Temperature Distribution at t = %.2f s', snapshot_times(i)));
    xlabel('x (m)');
    ylabel('y (m)');
    axis square;
    grid on;
    hold off;
    
    saveas(gcf, sprintf('Enhanced_Discrete_Contour_Temperature_t_%.2f_s.png', snapshot_times(i)));
end

%% Functions Used

% Gaussian heat source function
function Q = gaussianHeatSource(X, Y, x0, y0, P, r0)
    % Returns a 2D Gaussian heat source centered at (x0, y0)
    Q = (2*P/(pi*r0^2)) .* exp(-2*((X - x0).^2 + (Y - y0).^2) / r0^2);
end

% Temperature field update function
function T_new = computeTemp(T, Q, alpha, dt, dx, rho, cp, T_ambient, h)
    % Updates temperature field using an explicit finite difference method
    [Nx, Ny] = size(T);
    T_new = T;
    for i = 2:Nx-1
        for j = 2:Ny-1
            d2Tdx2 = (T(i+1,j) - 2*T(i,j) + T(i-1,j)) / dx^2;
            d2Tdy2 = (T(i,j+1) - 2*T(i,j) + T(i,j-1)) / dx^2;
            T_new(i,j) = T(i,j) + dt * ( ...
                alpha*(d2Tdx2 + d2Tdy2) + ...         % Thermal diffusion
                Q(i,j)/(rho*cp) - ...                 % Laser heating
                h*(T(i,j) - T_ambient)/(rho*cp) );       % Convective heat loss
        end
    end

    % Apply Neumann boundary conditions (zero flux)
    T_new(1,:)   = T_new(2,:);
    T_new(end,:) = T_new(end-1,:);
    T_new(:,1)   = T_new(:,2);
    T_new(:,end) = T_new(:,end-1);
end