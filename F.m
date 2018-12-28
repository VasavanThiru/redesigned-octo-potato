function Q_prime = F(Q, i, alpha, v_e, M_i, theta)
    % i: etape
    % alpha: acceleration initiale
    % v_e: vitesse d'ejection
    % M_i: masse initiale a l'etape i
    % theta: angles
    % c_x: coefficient de trainee du lanceur
    % rho_0: densite de l'air au sol
    % R_t: rayon terrestre
    R_t = 6378137;
    mu = 3.986e14;
    c_x = 0.1;
    rho_0 = 1.225;
    H = 7000; % facteur d’échelle
    e_r = 1 / norm(Q(1:2), 2) * Q(1:2);
    e_h = [-e_r(2); e_r(1)];
    gamma = asin(dot(Q(1:2), Q(3:4)) / (norm(Q(1:2), 2) * norm(Q(3:4), 2)));
    u = cos(gamma + theta) * e_h + sin(gamma + theta) * e_r;
    Q_prime = zeros(size(Q));
    Q_prime(1:2) = Q(3:4);
    Q_prime(3:4) = -mu * Q(1:2) / norm(Q(1:2), 2)^3 + alpha(i) * u - c_x * rho_0 * exp(-(norm(Q(1:2), 2) - R_t) / H) / Q(5) * norm(Q(3:4), 2) * Q(3:4);
    Q_prime(5) = -alpha(i) / v_e(i) * M_i;
    % Q_prime = Q_prime';
end
