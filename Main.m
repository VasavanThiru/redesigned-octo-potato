% Test de ode_integration
clear;
clear all;

% Constants
mu = 3.986e14;

% Parametres
m_u = 1500;
R = 200e3;
V = sqrt(mu / R);
m_0 = [1.8e5; 5e4; 1.5e4];
% acceleration, vitesse, indice
steps = [15, 2600, 0.1; 10, 3000, 0.15; 10, 4400, 0.20];
m_3 = @(m) steps(3, 2) * log((m_u+(1+steps(3, 1))*m(3))/(m_u+steps(3, 1)*m(3)));
m_2 = @(m) steps(2, 2) * log((m_u+(1+steps(3, 1))*m(3)+(1+steps(2, 1))*m(2))/(m_u+(1+steps(3, 1))*m(3)+steps(2, 1)*m(2)));
m_1 = @(m) steps(1, 2) * log((m_u+(1+steps(3, 1))*m(3)+(1+steps(2, 1))*m(2)+(1+steps(1, 1))*m(1))/(m_u+(1+steps(3, 1))*m(3)+(1+steps(2, 1))*m(2)+steps(1, 1)*m(1)));
F = @(m) [m_u + dot((1 + steps(1:3, 1)), m); -V + m_1(m) + m_2(m) + m_3(m)];
domain = [[1e5, 2e5]; [1e4, 5e4]; [5e3, 1e4]];

printf("Calcul des masses d'ergols\n");
[m_e, ~, ~, k] = SQP(m_0, F, 1e-6, domain, 30, 1, 0);
theta_0 = [1, 1, 1, 1];
domain_theta = [[-pi / 2, pi / 2]; [-pi / 2, pi / 2]; [-pi / 2, pi / 2]; [-pi / 2, pi / 2]];
G = @(theta) ode_integration(m_u, [m_e', m_u], theta, R, 0);
G(theta_0)
printf("Calcul des angles\n");
[theta, ~, ~, k] = SQP(theta_0, G, 1e-6, domain_theta, 30, 1, 0);
ode_integration(m_u, m_e, theta, R, 1)
