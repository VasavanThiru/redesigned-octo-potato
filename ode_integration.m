function []=ode_integration(R_0, v_0, M_0, theta)
    steps = [15, 2600; 10, 3000; 10, 4400];
    k = 0;
    R_0 = [R_t; 0];
    hold on;
    while k < length(steps)
        f = @(t, q) F(q, k + 1, step(1:3, 1), step(1:3, 2), theta, c_x, rho_0, R_t);
        [t, Q] = ode45(f, [R_0, V(k + 1), M(k + 1)], 
    end
    hold off;
end

% rho_0 = 1.225
% H = 7000
% R_t = 6378137
% mu = 3.986e14
% c_x = 0.1

% parameters:
% R_0 = [R_t; 0];
% V_0 = v_0 * [cos(theta); sin(theta)];
% M_0 = 1500;
% v_e = [2600, 3000, 4400];
% v = @(theta_0, theta_1, theta_2, theta_3)
