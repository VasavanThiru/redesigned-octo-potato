function []=ode_integration(R_0, v_0, M_0, theta)
    steps = [15, 2600; 10, 3000; 10, 4400];
    k = 0;
    hold on;
    while k < length(steps)
        [t, Q] = ode45(@(t, Q) [Q(4); Q(5); Q(6); 0
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
% v = @(theta_0, theta_1, theta_2, theta_3)
