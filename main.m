clc; clearvars; close all;
fprintf('start\n')

% solver settings
max_iter = 100;
CFL = 0.9;
N_i = 10; % number of nodes in R direction
N_j = 10; % number of nodes in Phi direction
GS_max_iter = 1000;
GS_tol = 1e-6;

% settings
U0 = 1;
T_h_absolut = 400;
T_c_absolut = 300;
R = 1;
R_a = R;
R_i = R / 2;
Phi_start = pi;
Phi_end = 2*pi;
Re = 60;
Pr = 4;
nu = U0 * R / Re;

% add ghost cells
i_max = N_i + 2;
j_max = N_j + 2;
% wall temperatures
T_h = (T_h_absolut - 300) / (400 - 300);
T_c = (T_c_absolut - 300) / (400 - 300);
% wall velocity
omega = 2 * pi;
v_wall = omega * R_i;

% grid
dr = (R_a - R_i) / N_i;
rp = linspace(R_i - dr/2, R_a + dr/2, i_max);
ru = linspace(R_i, R_a + dr, i_max);
dphi = (Phi_end - Phi_start) / N_j;
phip = linspace(Phi_start - dphi/2, Phi_end + dphi/2, j_max);
phiv = linspace(Phi_start, Phi_end + dphi, j_max);
[R_coord, Phi_coord] = meshgrid(linspace(R_i-dr/2, R_a+dr/2, i_max), linspace(Phi_start-dphi/2, Phi_end+dphi/2, j_max));
R_coord = transpose(R_coord);
Phi_coord = transpose(Phi_coord);

% small functions
avg = @(a, b) (a+b)/2;

% define variables
dt = zeros(i_max, j_max);

% n + 1 timestep
u_np1 = NaN(i_max, j_max);
v_np1 = NaN(i_max, j_max);
p_np1 = NaN(i_max, j_max);
T_np1 = NaN(i_max, j_max);

% current timestep
u_n = zeros(i_max, j_max);
v_n = zeros(i_max, j_max);
p_n = ones(i_max, j_max);
T_n = zeros(i_max, j_max);
Fi_n = zeros(i_max, j_max);
Fj_n = zeros(i_max, j_max);
Fe_n = zeros(i_max, j_max);

% n - 1 timestep
u_nm1 = zeros(i_max, j_max);
v_nm1 = zeros(i_max, j_max);
p_nm1 = zeros(i_max, j_max);
T_nm1 = zeros(i_max, j_max);
Fi_nm1 = zeros(i_max, j_max);
Fj_nm1 = zeros(i_max, j_max);
Fe_nm1 = zeros(i_max, j_max);

% projection step
u_star = NaN(i_max, j_max);
v_star = NaN(i_max, j_max);
p_corr = NaN(i_max, j_max);

% apply initial boundaries
u_n = apply_u_boundary(u_n, i_max, j_max);
v_n = apply_v_boundary(v_n, i_max, j_max, v_wall);
T_n = apply_T_boundary(T_n, i_max, j_max, T_c, T_h);
p_n = apply_p_boundary(p_n, i_max, j_max);
u_nm1 = apply_u_boundary(u_nm1, i_max, j_max);
v_nm1 = apply_v_boundary(v_nm1, i_max, j_max, v_wall);
T_nm1 = apply_T_boundary(T_nm1, i_max, j_max, T_c, T_h);
p_nm1 = apply_p_boundary(p_nm1, i_max, j_max);

% static Eq. for poisson Eq.
A = eye(i_max*j_max, i_max*j_max);
for i = 2:i_max-1
    for j = 2:j_max-1
        idx = j + (i-1) * j_max;
        A(idx, idx) = - ru(i)*dphi/dr - ru(i-1)*dphi/dr - 2 * dr/(rp(i)*dphi);
        A(idx, idx+1) = ru(i)*dphi/dr;
        A(idx, idx-1) = ru(i-1)*dphi/dr;
        A(idx, idx+j_max) = dr/(rp(i)*dphi);
        A(idx, idx-j_max) = dr/(rp(i)*dphi);
    end
end

% pressure boundaries in A matrix
for i = 2:i_max-1
    for j = 1:j_max
        idx = j + (i-1) * j_max;
        if j == 1 && i ~= 2 % no pressure gradient - B1
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
    % calculate dt per cell
    for i = 2:i_max-1
        for j = 2:j_max-1
            dt(i,j) = CFL/(abs(u_n(i,j))/dr + abs(v_n(i,j))/(dphi*rp(i)) + 2*nu/dr^2 + 2*nu/(dphi*rp(i))^2);
        end
    end

    % %%%%%%%% PROJECTION %%%%%%%%
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

            % projection
            u_star(i,j) = dt(i,j)/(dr*dphi*ru(i)) * (3/2 * Fi_n(i,j) - 1/2 * Fi_nm1(i,j)) + u_n(i,j);
        end
    end
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

            % projection
            v_star(i,j) = dt(i,j)/(dr*dphi*ru(i)) * (3/2 * Fj_n(i,j) - 1/2 * Fj_nm1(i,j)) + v_n(i,j);
        end
    end

    % apply boundaries
    u_star = apply_u_boundary(u_star, i_max, j_max);
    v_star = apply_v_boundary(v_star, i_max, j_max, v_wall);

    % %%%%%%%% POISSON %%%%%%%%
    p_corr_tmp = zeros(i_max*j_max, 1);
    b = zeros(i_max*j_max, 1);
    for i = 2:i_max-1
        for j = 2:j_max-1
            idx = j + (i-1) * j_max;
            b(idx) = ((ru(i) * u_star(i,j) - ru(i-1) * u_star(i-1,j)) * dphi ...
                + (v_star(i,j) - v_star(i,j-1)) * dr) / dt(i,j);
        end
    end

    % gauss seidl solver for solving the A*x=b Eq. system
%     x = p_corr_tmp;
%     for k = 1:GS_max_iter
%         err = 0;
%         for idx = 1:i_max*j_max
%             sum = 0;
%             for jdx = 1:i_max*j_max
%                 if jdx ~= idx
%                     sum = sum + A(idx,jdx) * x(jdx);
%                 end
%             end
%             x(idx) = (b(idx) - sum)/A(idx,idx);
%             err = max(err, abs(x(idx) - p_corr_tmp(idx)));
%         end
% 
%         if err <= GS_tol
%             break;
%         end
%         p_corr_tmp = x;
%     end
    p_corr_tmp = inv(A) * b;

    % reshape solution into correct dimensions
    p_corr = reshape(p_corr_tmp, [j_max, i_max])';

    % %%%%%%%% CORRECTION %%%%%%%%
    for i = 2:i_max-1
        for j = 2:j_max-1
            % RHS for energy
            e_t1 = dphi * (rp(i) * u_n(i,j) * (T_n(i+1,j) + T_n(i,j))/2 ...
                - rp(i-1) * u_n(i-1,j) * (T_n(i,j) + T_n(i-1,j))/2);
            e_t2 = dr * (v_n(i,j) * (T_n(i,j+1) + T_n(i,j))/2 ...
                - v_n(i,j-1) * (T_n(i,j) + T_n(i,j-1))/2);
            e_t3 = 1/(Pr*Re) * ((ru(i) * (T_n(i+1,j)-T_n(i,j))/dr - ru(i-1) * (T_n(i,j)-T_n(i-1,j))/dr) * dphi ...
                + dr/rp(i) * ((T_n(i,j+1) - T_n(i,j))/dphi - (T_n(i,j) - T_n(i,j-1))/dphi));
            Fe_n(i,j) = - e_t1 - e_t2 + e_t3;

            T_np1(i,j) = dt(i,j)/(dr*dphi*ru(i)) * (3/2 * Fe_n(i,j) - 1/2 * Fe_nm1(i,j)) + T_n(i,j);
        end
    end
    for i = 2:i_max-2
        for j = 2:j_max-1
            % correction
            u_np1(i,j) = u_star(i,j) - dt(i,j)/dr * (p_corr(i+1,j) - p_corr(i,j));
        end
    end
    for i = 2:i_max-1
        for j = 2:j_max-2
            % correction
            v_np1(i,j) = v_star(i,j) - dt(i,j)/dphi * (p_corr(i,j+1) - p_corr(i,j));
        end
    end
    p_np1 = p_n + p_corr;

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

    if mod(itr, 10) == 0
        fprintf('iter: %i\n', itr)
    end
end

fprintf('completed\n')

% %%%%%%%% PLOTS %%%%%%%%
figure
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
X = R_coord .* cos(Phi_coord);
Y = R_coord .* sin(Phi_coord);
u_x = u_n .* cos(Phi_coord);
u_y = u_n .* sin(Phi_coord);
v_x = -v_n .* sin(Phi_coord);
v_y = v_n .* cos(Phi_coord);
quiver(X, Y, u_x + v_x, u_y + v_y);
axis equal;

figure
contourf(X,Y,T_n,100,'LineColor','none');
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
axis equal;

figure
contourf(X,Y,p_n,100,'LineColor','none');
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
axis equal;


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