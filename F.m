function Q_prime = F(Q, alpha, theta)
    Q_prime(1:3) = Q(4:6);
    Q_prime(4:6) = -mu * Q(1:3) / norm(Q(1:3), 2)^3 + alpha(1) * u +
end
