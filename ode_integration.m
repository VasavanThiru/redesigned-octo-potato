function VC = ode_integration(m_0, m_e, theta, R, draw)
    % TODO: comments and presentation/formatting
    R_t = 6378137;
    % acceleration, vitesse, indice
    steps = [15, 2600, 0.1; 10, 3000, 0.15; 10, 4400, 0.20];

    R_i = zeros(4, 2);
    R_i(1, :) = [R_t, 0];
    % m_e = [145349, 31215, 7933, 1700]; % TODO: m_e = 1500
    M_i(4) = m_0;
    % calcul des masses
    for j = 3:-1:1
        M_f(j) = M_i(j + 1) + steps(j, 3) * m_e(j);
        M_i(j) = M_f(j) + m_e(j);
    end
    % calcul des durees
    t = zeros(1, 3);
    for j = 2:4
        t(j) = steps(j - 1, 2) / (steps(j - 1, 1) * M_i(j)) * m_e(j);
    end
    if draw == 1
        hold on;
        axis equal;
    end
    k = 0;
    v = zeros(4, 2);
    v(1, :) = 100 * [cos(theta(k + 1)), sin(theta(k + 1))];
    while k < length(steps)
        f = @(t, q) F(q, k + 1, steps(1:3, 1), steps(1:3, 2), M_i(k + 1), theta(k + 1));
        %f(0, [R_0', v(k + 1) * [cos(theta(k + 1)), sin(theta(k + 1))], M_i(k + 1)])
        [T, Q] = ode45(f, t(k + 1):0.01:t(k + 2), [R_i(k + 1, :), v(k + 1, :), M_i(k + 1)]);
        n = length(T);
        R_i(k + 2, :) = Q(n, 1:2);
        v(k + 2, :) = Q(n, 3:4);
        if draw == 1
            plot(Q(:, 1), Q(:, 2));
        end
        k = k + 1;
    end
    if draw == 1
        T = 0:0.001:1;
        plot(R_t * T, R_t * sqrt(1 - T.^2));
        plot(R_t * T, -R_t * sqrt(1 - T.^2));
        R_t = R_t + R;
        plot(R_t * T, R_t * sqrt(1 - T.^2));
        plot(R_t * T, -R_t * sqrt(1 - T.^2));
        hold off;
    end
    N = length(t);
    V = -dot(Q(N, 3:4), Q(N, 3:4)); % valeur
    C = zeros(2, 1);
    C(1) = dot(Q(N, 1:2), Q(N, 3:4)); % contraintes
    C(2) = R - norm(Q(N, 1:2), 2);
    VC = [V; C];
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
