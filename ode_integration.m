function []=ode_integration(m_0, theta, R)
    mu = 3.986e14;
    R_t = 6378137;
    steps = [15, 2600, 0.1; 10, 3000, 0.15; 10, 4400, 0.20]; % acceleration, vitesse, indice

    k = 0;
    R_0 = [R_t; 0];
    x_m
    V_p = sqrt(mu / R);
    mass = @(x) [-prod(steps(1, 1:3) ./ x - step(1, 1:3)), dot(step(2, 1:3), log(x)) - V_p];
    [x, ~, ~, ~] = SQP(x_0, mass, 0.001, domain, 1000, 1);
    M(4) = m_0;
    for i = 3:-1:1
        M(i) = M(i+1) / x(i);
    end
    t(0) = 0;
    hold on;
    while k < length(steps)
        f = @(t, q) F(q, k + 1, step(1:3, 1), step(1:3, 2), theta, c_x, rho_0, R_t);
        [t, Q] = ode45(f, [R_0, V(k + 1), M(k + 1)], [t(k), t(k+1)]);
        plot(t, Q);
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
%
