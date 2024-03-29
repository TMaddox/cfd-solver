clc; clearvars; close all; % clear console; clear varibles; figures schliessen
fprintf('start\n')

% solver settings
max_iter = 10000;       % max. # Iterationen
CFL = 0.1;              % Courant-Friedrichs-Lewy Kriterium
N_i = 15;               % # of nodes in R-Richtung
N_j = 4 * N_i;          % # of nodes in Phi-Richtung
alpha = 0.2;            % Relaxation für velocity
beta = 0.2;             % Relaxation für pressure
GS_max_iter = 100;      % max. # Iterationen für Gauss-Seidel solver

% settings
T_h_absolut = 400;      % max. Temperatur (Beispielhaft, normierung erfolgt später)
T_c_absolut = 300;      % min. Temperatur (Beispielhaft, normierung erfolgt später)
R = 1;                  % Radius (Beispielhaft)
R_a = R;                % Radius aussen (lt. Angabe)
R_i = R / 2;            % Radius innen (lt. Angabe)
Phi_start = pi;         % Startwinkel (lt. Angabe)
Phi_end = 2*pi;         % Endwinkel (lt. Angabe)
Re = 60;                % Reynolds Zahl (lt. Angabe)
Pr = 4;                 % Prandtl Zahl (lt. Angabe)
display_ghost_cells = false;    % specify if ghost cells should be displayed

% add ghost cells
i_max = N_i + 2;
j_max = N_j + 2;
i_fixed = 3; % one p_corr needs to be set to a constant fixed value because only Neumann boundaries are present
j_fixed = 3;

% wall temperatures
T_h = (T_h_absolut - T_c_absolut) / (T_h_absolut - T_c_absolut);    % normierung Temperatur
T_c = (T_c_absolut - T_c_absolut) / (T_h_absolut - T_c_absolut);    % normierung Temperatur
% wall velocity
v_wall_max = 1;         % Geschwindigkeit bewegte Wand (Beispielhaft)

% Grid Parameter
dr = (R_a - R_i) / N_i;                             % Grid size in R-Richtung
rp = linspace(R_i - dr/2, R_a + dr/2, i_max);       % linspace von allen rp's (Radius aller Points)
ru = linspace(R_i, R_a + dr, i_max);                % linspace von allen ru's (Radius aller Zellwände)
dphi = (Phi_end - Phi_start) / N_j;                 % Grid size in phi-Richtung
phip = linspace(Phi_start - dphi/2, Phi_end + dphi/2, j_max);   % linspace von allen phip's (Winkel aller Points)
phiv = linspace(Phi_start, Phi_end + dphi, j_max);              % linspace von allen phiv's (Winkel aller Zellwände)
[R_coord, Phi_coord] = meshgrid(linspace(R_i-dr, R_a+dr, i_max), linspace(Phi_start-dphi, Phi_end+dphi, j_max)); % Meshgrid aller Nodes (incl. ghost cells)
R_coord_with_ghost = transpose(R_coord);                        % Korrektur
Phi_coord_with_ghost = transpose(Phi_coord);                    % Korrektur
if ~display_ghost_cells
    R_coord = R_coord_with_ghost(2:end-1, 2:end-1);             % cut out Ghost cells
    Phi_coord = Phi_coord_with_ghost(2:end-1, 2:end-1);         % cut out ghost cells
end

% small functions
avg = @(a, b) (a+b)/2;          % average between two values

% other variables for plotting
fconti = NaN(max_iter, 1);      % residual_der_kontinuitätsgleichung(iter)
fp_corr = NaN(max_iter, 1);     % p_corr(iter)

% time difference per cell
dt = NaN(i_max, j_max);         % Zeitschritt per cell

% n + 1 timestep
u_np1 = NaN(i_max, j_max);      % Geschwindigkeit u (in r-Richtung)
v_np1 = NaN(i_max, j_max);      % Geschwindigkeit v (in phi-Richtung)
p_np1 = NaN(i_max, j_max);      % Druck
T_np1 = NaN(i_max, j_max);      % Temperatur

% current (n) timestep
u_n = ones(i_max, j_max) * 0;   % Geschwindigkeit u (in r-Richtung)
v_n = ones(i_max, j_max) * 0;   % Geschwindigkeit v (in phi-Richtung)
p_n = ones(i_max, j_max) * 0;   % Druck
T_n = ones(i_max, j_max) * 0.5; % Temperatur
Fi_n = zeros(i_max, j_max);     % RHS der Impulserhaltungsgleichung in R-Richtung zum Zeitpunkt n
Fj_n = zeros(i_max, j_max);     % RHS der Impulserhaltungsgleichung in phi-Richtung zum Zeitpunkt n
Fe_n = zeros(i_max, j_max);     % RHS der Energiegleichung zum Zeitpunkt n

% n - 1 timestep
u_nm1 = ones(i_max, j_max) * 0; % Geschwindigkeit u (in r-Richtung)
v_nm1 = ones(i_max, j_max) * 0; % Geschwindigkeit v (in phi-Richtung)
p_nm1 = ones(i_max, j_max) * 0; % Druck
T_nm1 = ones(i_max, j_max) * 0.5;   % Temperatur
Fi_nm1 = zeros(i_max, j_max);   % RHS der Impulserhaltungsgleichung in R-Richtung zum Zeitpunkt n-1
Fj_nm1 = zeros(i_max, j_max);   % RHS der Impulserhaltungsgleichung in phi-Richtung zum Zeitpunkt n-1
Fe_nm1 = zeros(i_max, j_max);   % RHS der Energiegleichung zum Zeitpunkt n-1

% projection step
u_star = NaN(i_max, j_max);     % Geschwindigkeit u* (in r-Richtung) im Zeitpunkt *
v_star = NaN(i_max, j_max);     % Geschwindigkeit v* (in phi-Richtung) im Zeitpunkt *
p_corr = NaN(i_max, j_max);     % Druck für Druckkorrektur (für * --> n+1)

% apply initial boundaries
u_n = apply_u_boundary(u_n, i_max, j_max);
v_n = apply_v_boundary(v_n, i_max, j_max, 0);
T_n = apply_T_boundary(T_n, i_max, j_max, T_c, T_h);
p_n = apply_p_boundary(p_n, i_max, j_max);
u_nm1 = apply_u_boundary(u_nm1, i_max, j_max);
v_nm1 = apply_v_boundary(v_nm1, i_max, j_max, 0);
T_nm1 = apply_T_boundary(T_nm1, i_max, j_max, T_c, T_h);
p_nm1 = apply_p_boundary(p_nm1, i_max, j_max);

% static Eq. for poisson Eq.
A = eye(i_max*j_max, i_max*j_max);
for i = 2:i_max-1
    for j = 2:j_max-1
        if i ~= i_fixed || j ~= j_fixed
            idx = j + (i-1) * j_max;
            A(idx, idx) = - ru(i)*dphi/dr - ru(i-1)*dphi/dr - 2*dr/(rp(i)*dphi);  % i,j
            A(idx, idx+1) = ru(i)*dphi/dr;          % i+1,j
            A(idx, idx-1) = ru(i-1)*dphi/dr;        % i-1,j
            A(idx, idx+j_max) = dr/(rp(i)*dphi);    % i,j+1
            A(idx, idx-j_max) = dr/(rp(i)*dphi);    % i,j-1
        end
    end
end

% pressure boundaries in A matrix
for i = 2:i_max-1
    for j = 1:j_max
        idx = j + (i-1) * j_max;
        if j == 1 % no pressure gradient - B1
            A(idx, idx+1) = -1;
        elseif j == j_max % no pressure gradient - B3
            A(idx, idx-1) = -1;
        end
    end
end
for i = 1:i_max
    for j = 2:j_max-1
        idx = j + (i-1) * j_max;
        if i == 1 % no pressure gradient - B2
            A(idx, idx+j_max) = -1;
        elseif i == i_max % no pressure gradient - B4
            A(idx, idx-j_max) = -1;
        end
    end
end

% %%%%%%%% ITERATIVE SOLVING %%%%%%%%
for itr = 1:max_iter
    if mod(itr, 100) == 0 || itr == 1
        fprintf('start iter: %i\n', itr)

        % %%%%%%%% PLOTS %%%%%%%%
        % velocity
        figure(1)
        n_plot = 100;
        phi_plot = linspace(Phi_start, Phi_end, n_plot);
        x_plot = zeros(2*n_plot+1,1);
        x_plot(1:n_plot) = R_i * cos(phi_plot);
        x_plot(n_plot+1:2*n_plot) = - R_a * cos(phi_plot);
        x_plot(2*n_plot+1) = x_plot(1);
        y_plot = zeros(2*n_plot+1,1);
        y_plot(1:n_plot) = R_i * sin(phi_plot);
        y_plot(n_plot+1:2*n_plot) = R_a * sin(phi_plot);
        y_plot(2*n_plot+1) = y_plot(1);
        plot(x_plot, y_plot, 'k');
        
        hold on;
        % combined
        X = R_coord .* cos(Phi_coord);
        Y = R_coord .* sin(Phi_coord);
        if ~display_ghost_cells
            u_x = u_n(2:end-1, 2:end-1) .* cos(Phi_coord);
            u_y = u_n(2:end-1, 2:end-1) .* sin(Phi_coord);
            v_x = -v_n(2:end-1, 2:end-1) .* sin(Phi_coord);
            v_y = v_n(2:end-1, 2:end-1) .* cos(Phi_coord);
        else
            u_x = u_n .* cos(Phi_coord);
            u_y = u_n .* sin(Phi_coord);
            v_x = -v_n .* sin(Phi_coord);
            v_y = v_n .* cos(Phi_coord);
        end
        plot_vel_x = u_x + v_x;
        plot_vel_y = u_y + v_y;
        quiver(X(1:2:end, 1:2:end), Y(1:2:end, 1:2:end), plot_vel_x(1:2:end, 1:2:end), plot_vel_y(1:2:end, 1:2:end), 'r', 'linewidth', 1);

        % u
%         X = (R_coord + dr/2) .* cos(Phi_coord);
%         Y = (R_coord + dr/2) .* sin(Phi_coord);
%         u_x = u_n .* cos(Phi_coord);
%         u_y = u_n .* sin(Phi_coord);
%         quiver(X, Y, u_x, u_y);

        % v
%         X = R_coord .* cos(Phi_coord + dphi/2);
%         Y = R_coord .* sin(Phi_coord + dphi/2);
%         v_x = -v_n .* sin(Phi_coord + dphi/2);
%         v_y = v_n .* cos(Phi_coord + dphi/2);
%         quiver(X, Y, v_x, v_y);
        title("Velocity");
        xlabel ('x');
        ylabel('y');
        hold off;
        axis equal;
        
        % temperature
        figure(2)
        X = R_coord .* cos(Phi_coord);
        Y = R_coord .* sin(Phi_coord);
        if ~display_ghost_cells
            contourf(X,Y,T_n(2:end-1, 2:end-1),30,'LineColor','none');
        else
            contourf(X,Y,T_n,30,'LineColor','none');
        end
        colorbar;
        colormap('jet');
        
        hold on;
        n_plot = 100;
        phi_plot = linspace(Phi_start, Phi_end, n_plot);
        x_plot = zeros(2*n_plot+1,1);
        x_plot(1:n_plot) = R_i * cos(phi_plot);
        x_plot(n_plot+1:2*n_plot) = - R_a * cos(phi_plot);
        x_plot(2*n_plot+1) = x_plot(1);
        y_plot = zeros(2*n_plot+1,1);
        y_plot(1:n_plot) = R_i * sin(phi_plot);
        y_plot(n_plot+1:2*n_plot) = R_a * sin(phi_plot);
        y_plot(2*n_plot+1) = y_plot(1);
        plot(x_plot, y_plot, 'k');
        title("Temperature");
        xlabel ('x');
        ylabel('y');
        hold off;
        axis equal;
        
        % pressure
        figure(3)
        X = R_coord .* cos(Phi_coord);
        Y = R_coord .* sin(Phi_coord);
        if ~display_ghost_cells
            contourf(X,Y,p_n(2:end-1, 2:end-1),100,'LineColor','none');
        else
            contourf(X,Y,p_n,30,'LineColor','none');
        end
        colorbar;
        colormap('cool');
        
        hold on;
        n_plot = 100;
        phi_plot = linspace(Phi_start, Phi_end, n_plot);
        x_plot = zeros(2*n_plot+1,1);
        x_plot(1:n_plot) = R_i * cos(phi_plot);
        x_plot(n_plot+1:2*n_plot) = - R_a * cos(phi_plot);
        x_plot(2*n_plot+1) = x_plot(1);
        y_plot = zeros(2*n_plot+1,1);
        y_plot(1:n_plot) = R_i * sin(phi_plot);
        y_plot(n_plot+1:2*n_plot) = R_a * sin(phi_plot);
        y_plot(2*n_plot+1) = y_plot(1);
        plot(x_plot, y_plot, 'k');
        title("Pressure");
        xlabel ('x');
        ylabel('y');
        hold off;
        axis equal;

        if itr ~= 1
            % streamline
            figure(4)
            clf(4)
            X = R_coord_with_ghost .* cos(Phi_coord_with_ghost);
            Y = R_coord_with_ghost .* sin(Phi_coord_with_ghost);
%             u_x = u_n .* cos(Phi_coord);
%             u_y = u_n .* sin(Phi_coord);
%             v_x = -v_n .* sin(Phi_coord);
%             v_y = v_n .* cos(Phi_coord);
            if ~display_ghost_cells
                vel_x = NaN(i_max,j_max);
                vel_y = NaN(i_max,j_max);
                vel_x(2:i_max-1, 2:j_max-1) = u_x + v_x;
                vel_y(2:i_max-1, 2:j_max-1) = u_y + v_y;
            else
                vel_x = NaN(i_max,j_max);
                vel_y = NaN(i_max,j_max);
                vel_x(2:i_max-1, 2:j_max-1) = u_x(2:i_max-1, 2:j_max-1) + v_x(2:i_max-1, 2:j_max-1);
                vel_y(2:i_max-1, 2:j_max-1) = u_y(2:i_max-1, 2:j_max-1) + v_y(2:i_max-1, 2:j_max-1);
            end

            num_x = 1000;
            num_y = num_x/2;
            [X_cart, Y_cart] = meshgrid(linspace(min(min(X)), max(max(X)), num_x), linspace(max(max(Y)), min(min(Y)), num_y));
            vel_x_cart = griddata(X(:), Y(:), vel_x(:), X_cart(:), Y_cart(:));
            vel_y_cart = griddata(X(:), Y(:), vel_y(:), X_cart(:), Y_cart(:));
            vel_x_cart = reshape(vel_x_cart, num_y, num_x);
            vel_y_cart = reshape(vel_y_cart, num_y, num_x);
            streamslice(X_cart, Y_cart, vel_x_cart, vel_y_cart);

            hold on;
            n_plot = 100;
            phi_plot = linspace(Phi_start, Phi_end, n_plot);
            x_plot = zeros(2*n_plot+1,1);
            x_plot(1:n_plot) = R_i * cos(phi_plot);
            x_plot(n_plot+1:2*n_plot) = - R_a * cos(phi_plot);
            x_plot(2*n_plot+1) = x_plot(1);
            y_plot = zeros(2*n_plot+1,1);
            y_plot(1:n_plot) = R_i * sin(phi_plot);
            y_plot(n_plot+1:2*n_plot) = R_a * sin(phi_plot);
            y_plot(2*n_plot+1) = y_plot(1);
            plot(x_plot, y_plot, 'k');
            title("Streamlines");
            xlabel ('x');
            ylabel('y');
            hold off;
            axis equal;
        end

        % plot conti
        figure(5)
        semilogy(linspace(1,itr-1,itr-1), fconti(1:itr-1));
        xlabel ('itr');
        ylabel('max(abs(conti))');
        title("max(abs(conti))(itr)");

        % plot conti
        figure(6)
        semilogy(linspace(1,itr-1,itr-1), fp_corr(1:itr-1));
        xlabel ('itr');
        ylabel('max(abs(p_corr))');
        title("max(abs(p corr))(itr)");

        drawnow;
    end

    % calculate v_wall as a function of iter for smooth start
    v_wall = v_wall_max * (1 - exp(-itr/(max_iter/20)));

    % calculate dt per cell (local timestep)
    for i = 2:i_max-1
        for j = 2:j_max-1
            dt(i,j) = CFL/(abs(u_n(i,j))/dr + abs(v_n(i,j))/(dphi*rp(i)) + 2/dr^2/Re + 2/(dphi*rp(i))^2/Re);
        end
    end

    % %%%%%%%% PROJECTION %%%%%%%%
    % RHS for r-impulse & projection u*
    for i = 2:i_max-2
        for j = 2:j_max-1
            % RHS for r-impulse
            i_t1 = (rp(i+1) * avg(u_n(i+1,j), u_n(i,j))^2 ...
                - rp(i) * avg(u_n(i,j), u_n(i-1,j))^2) * dphi;
            i_t2 = (avg(u_n(i,j+1), u_n(i,j)) * avg(v_n(i+1,j), v_n(i,j)) ...
                - avg(u_n(i,j), u_n(i,j-1)) * avg(v_n(i+1,j-1), v_n(i,j-1))) * dr;
            i_t3 = ((v_n(i,j) + v_n(i+1,j) + v_n(i+1,j-1) + v_n(i,j-1))/4)^2 * dr * dphi;
            i_t4 = (p_n(i+1,j) - p_n(i,j)) * ru(i) * dphi;
            i_t5 = 2/Re * (rp(i+1) * (u_n(i+1,j) - u_n(i,j))/dr ...
                - rp(i) * (u_n(i,j) - u_n(i-1,j))/dr) * dphi;
            i_t6 = 1/Re * (ru(i) * (v_n(i+1,j)/rp(i+1) - v_n(i,j)/rp(i))/dr ...
                - ru(i) * (v_n(i+1,j-1)/rp(i+1) - v_n(i,j-1)/rp(i))/dr ...
                + (u_n(i,j+1) - u_n(i,j)) / (ru(i)*dphi) ...
                - (u_n(i,j) - u_n(i,j-1)) / (ru(i)*dphi)) * dr;
            i_t7 = 2/(Re*ru(i)) * ...
                ((avg(v_n(i+1,j), v_n(i,j)) - avg(v_n(i+1,j-1), v_n(i,j-1))) / dphi ...
                + u_n(i,j)) * dr * dphi;
            Fi_n(i,j) = - i_t1 - i_t2 + i_t3 - i_t4 + i_t5 + i_t6 - i_t7;

            % projection u*
            u_star(i,j) = dt(i,j)/(dr*dphi*ru(i)) * (3/2 * Fi_n(i,j) - 1/2 * Fi_nm1(i,j)) + u_n(i,j);
        end
    end
    % RHS fpr phi-impulse & projection v*
    for i = 2:i_max-1
        for j = 2:j_max-2
            % RHS for phi-impulse
            j_t1 = (ru(i) * avg(u_n(i,j+1), u_n(i,j)) * avg(v_n(i+1,j), v_n(i,j)) ...
                - ru(i-1) * avg(u_n(i-1,j+1), u_n(i-1,j)) * avg(v_n(i,j), v_n(i-1,j))) * dphi;
            j_t2 = (avg(v_n(i,j+1), v_n(i,j))^2 ...
                - avg(v_n(i,j), v_n(i,j-1))^2) * dr;
            j_t3 = ((u_n(i,j) + u_n(i,j+1) + u_n(i-1,j) + u_n(i-1,j+1))/4) * v_n(i,j) * dr * dphi;
            j_t4 = (p_n(i,j+1) - p_n(i,j)) * dr;
            j_t5 = 1/(Re * rp(i)) * (ru(i)^3 * (v_n(i+1,j)/rp(i+1) - v_n(i,j)/rp(i))/dr ...
                - ru(i-1)^3 * (v_n(i,j)/rp(i) - v_n(i-1,j)/rp(i-1))/dr ...
                + ru(i) * (u_n(i,j+1) - u_n(i,j))/dphi ...
                - ru(i-1) * (u_n(i-1,j+1) - u_n(i-1,j))/dphi) * dphi;
            j_t6 = 2/(Re * rp(i)) * ((v_n(i,j+1) - v_n(i,j))/dphi ...
                - (v_n(i,j) - v_n(i,j-1))/dphi ...
                + avg(u_n(i,j+1), u_n(i-1,j+1)) ...
                - avg(u_n(i,j), u_n(i-1,j))) * dr;
            Fj_n(i,j) = - j_t1 - j_t2 - j_t3 - j_t4 + j_t5 + j_t6;

            % projection v*
            v_star(i,j) = dt(i,j)/(dr*dphi*rp(i)) * (3/2 * Fj_n(i,j) - 1/2 * Fj_nm1(i,j)) + v_n(i,j);
        end
    end

    % apply boundaries for projected velocities
    u_star = apply_u_boundary(u_star, i_max, j_max);
    v_star = apply_v_boundary(v_star, i_max, j_max, v_wall);

    % %%%%%%%% POISSON %%%%%%%%
    b = zeros(i_max*j_max, 1);
    for i = 2:i_max-1
        for j = 2:j_max-1
            if i ~= i_fixed || j ~= j_fixed
                idx = j + (i-1) * j_max;
                b(idx) = ((ru(i) * u_star(i,j) - ru(i-1) * u_star(i-1,j)) * dphi ...
                    + (v_star(i,j) - v_star(i,j-1)) * dr) / dt(i,j);
            end
        end
    end

    % gauss seidl solver for solving the A*x=b Eq. system
    p_corr_tmp = zeros(i_max*j_max, 1);
    p_corr_tmp = GS(A, b, p_corr_tmp, GS_max_iter);

    % reshape solution into correct dimensions
    p_corr = reshape(p_corr_tmp, [j_max, i_max])';
    p_corr(1, 1) = nan;
    p_corr(i_max, 1) = nan;
    p_corr(1, j_max) = nan;
    p_corr(i_max, j_max) = nan;

    % %%%%%%%% CORRECTION %%%%%%%%
    % correct u
    for i = 2:i_max-2
        for j = 2:j_max-1
            u_np1(i,j) = u_star(i,j) - dt(i,j)/dr * (p_corr(i+1,j) - p_corr(i,j)) * alpha;
        end
    end
    % correct v
    for i = 2:i_max-1
        for j = 2:j_max-2
            v_np1(i,j) = v_star(i,j) - dt(i,j)/(dphi * rp(i)) * (p_corr(i,j+1) - p_corr(i,j)) * alpha;
        end
    end
    % correct pressure
    p_np1 = p_n + p_corr * beta;

    % %%%%%%%% CALCULATE ENERGY EQUATION %%%%%%%%
    for i = 2:i_max-1
        for j = 2:j_max-1
            % RHS for energy
            e_t1 = dphi * (rp(i) * u_n(i,j) * avg(T_n(i+1,j), T_n(i,j)) ...
                - rp(i-1) * u_n(i-1,j) * avg(T_n(i,j), T_n(i-1,j)));
            e_t2 = dr * (v_n(i,j) * avg(T_n(i,j+1), T_n(i,j)) ...
                - v_n(i,j-1) * avg(T_n(i,j), T_n(i,j-1)));
            e_t3 = 1/(Pr*Re) * ((ru(i) * (T_n(i+1,j)-T_n(i,j))/dr - ru(i-1) * (T_n(i,j)-T_n(i-1,j))/dr) * dphi ...
                + dr/rp(i) * ((T_n(i,j+1) - T_n(i,j))/dphi - (T_n(i,j) - T_n(i,j-1))/dphi));
            Fe_n(i,j) = - e_t1 - e_t2 + e_t3;

            T_np1(i,j) = dt(i,j)/(dr*dphi*ru(i)) * (3/2 * Fe_n(i,j) - 1/2 * Fe_nm1(i,j)) + T_n(i,j);
        end
    end

    % %%%%%%%% STORE VALUES %%%%%%%%
    % copy values from timestep n to n-1 for next iter
    u_nm1 = u_n;
    v_nm1 = v_n;
    p_nm1 = p_n;
    T_nm1 = T_n;
    Fi_nm1 = Fi_n;
    Fj_nm1 = Fj_n;
    Fe_nm1 = Fe_n;

    % copy values from timestep n+1 to n for next iter
    u_n = u_np1;
    v_n = v_np1;
    p_n = p_np1;
    T_n = T_np1;

    % apply boundaries
    u_n = apply_u_boundary(u_n, i_max, j_max);
    v_n = apply_v_boundary(v_n, i_max, j_max, v_wall);
    T_n = apply_T_boundary(T_n, i_max, j_max, T_c, T_h);
    p_n = apply_p_boundary(p_n, i_max, j_max);

    % %%%%%%%% check conti %%%%%%%%
    conti = NaN(i_max, j_max);
    for i = 2:i_max-1
        for j = 2:j_max-1
            conti(i, j) = ((ru(i) * u_n(i,j) - ru(i-1) * u_n(i-1,j)) * dphi ...
                + (v_n(i,j) - v_n(i,j-1)) * dr);
        end
    end
    fprintf('max(abs(conti)): %i\n', max(max(abs(conti))))
    %fprintf('max(abs(p_corr)): %i\n', max(max(abs(p_corr))))

    % track convergence metrics over iter
    fconti(itr) = max(max(abs(conti)));
    fp_corr(itr) = max(max(abs(p_corr)));
end

fprintf('\n--------------------\n')
fprintf('max(abs(conti)): %i\n', max(max(abs(conti))))
fprintf('max(abs(p_corr)): %i\n', max(max(abs(p_corr))))
fprintf('completed\n')

% %%%%%%%% calculate Nusselt number %%%%%%%%
Nu = 0;
for j = 2:j_max-1
    Nu = Nu + (T_n(i_max,j) - T_n(i_max-1,j))/dr * R * dphi; 
end
fprintf('Nu: %i\n', Nu)


% %%%%%%%% BOUNDARIES %%%%%%%%
% Boundary 1 - adiabat left
% Boundary 2 - inner radius
% Boundary 3 - adiabat right
% Boundary 4 - outer radius
function u = apply_u_boundary(u, i_max, j_max)
    u(1, 2:j_max-1) = 0; % B2
    u(i_max-1, 2:j_max-1) = 0; % B4
    u(2:i_max-2, 1) = - u(2:i_max-2, 2); % B1
    u(2:i_max-2, j_max) = - u(2:i_max-2, j_max - 1); % B3
end

function v = apply_v_boundary(v, i_max, j_max, v_wall)
    v(2:i_max-1, 1) = 0; % B1
    v(2:i_max-1, j_max-1) = 0; % B3
    v(1, 2:j_max-2) = 2 * v_wall - v(2, 2:j_max-2); % B2
    v(i_max, 2:j_max-2) = - v(i_max - 1, 2:j_max-2); % B4
end

function T = apply_T_boundary(T, i_max, j_max, T_c, T_h)
    T(2:i_max-1, 1) = T(2:i_max-1, 2); % B1
    T(1, 2:j_max-1) = 2 * T_c - T(2, 2:j_max-1); % B2
    T(2:i_max-1, j_max) = T(2:i_max-1, j_max - 1); % B3
    T(i_max, 2:j_max-1) = 2 * T_h - T(i_max - 1, 2:j_max-1); % B4
end

function p = apply_p_boundary(p, i_max, j_max)
    p(2:i_max-1, 1) = p(2:i_max-1, 2); % B1
    p(1, 2:j_max-1) = p(2, 2:j_max-1); % B2
    p(2:i_max-1, j_max) = p(2:i_max-1, j_max - 1); % B3
    p(i_max, 2:j_max-1) = p(i_max - 1, 2:j_max-1); % B4
end

function x = GS(A, b, x_init, max_iter)
    % gauss seidl solver for solving the A*x=b Eq. system
    x = x_init;
    for idx = 1:max_iter
        for jdx = 1:size(A,1)
            x(jdx) = (b(jdx) - sum(A(jdx,:)'.*x) + A(jdx,jdx)*x(jdx)) / A(jdx,jdx);
        end
    end
end