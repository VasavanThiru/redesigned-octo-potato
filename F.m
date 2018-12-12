function Q_prime = F(Q, i, alpha, v_e, M, theta, c_x, rho_0, R_t)
    % i: etape
    % alpha: acceleration initiale
    % v_e: vitesse d'ejection
    % M: masse initiale a l'etape i
    % theta: angles
    % c_x: coefficient de trainee du lanceur
    % rho_0: densite de l'air au sol
    % R_t: rayon terrestre
    Q_prime(1:2) = Q(3:4);
    Q_prime(3:4) = -mu * Q(1:2) / norm(Q(1:2), 2)^3 + alpha(i) * u - c_x * rho_0 * exp(-(norm(Q(1:2), 2) - R_t) / Q(5) * norm(Q(3, 4) * Q(3, 4)
    Q_prime(5) = -alpha(i) / v_e(i) * M(i);
end
