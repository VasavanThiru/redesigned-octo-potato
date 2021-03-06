% Test de ode_integration
clear;
clear all;

% Constants
mu = 3.986e14;
R_t = 6378137;

% Parametres
m_u = 1500;
R = 200e3;
V = sqrt(mu / (R_t + R));
m_0 = [1.8e5; 5e4; 1.5e4];
% acceleration, vitesse, indice
steps = [15, 2600, 0.1; 10, 3000, 0.15; 10, 4400, 0.20];

m_3 = @(m) steps(3, 2) * log((m_u+(1+steps(3, 3))*m(3))/(m_u+steps(3, 3)*m(3)));
m_2 = @(m) steps(2, 2) * log((m_u+(1+steps(3, 3))*m(3)+(1+steps(2, 3))*m(2))/(m_u+(1+steps(3, 3))*m(3)+steps(2, 3)*m(2)));
m_1 = @(m) steps(1, 2) * log((m_u+(1+steps(3, 3))*m(3)+(1+steps(2, 3))*m(2)+(1+steps(1, 3))*m(1))/(m_u+(1+steps(3, 3))*m(3)+(1+steps(2, 3))*m(2)+steps(1, 3)*m(1)));
F = @(m) [m_u + dot((1 + steps(1:3, 1)), m); -V + m_1(m) + m_2(m) + m_3(m)];
domain = [[1e5, 2e5]; [1e4, 5e4]; [5e3, 1e4]];

fprintf("Calcul des masses d'ergols\n");
[m_e, ~, ~, k] = SQP(m_0, F, 1e-6, domain, 30, 1, 0);

theta_0 = [1; 1; 1; 1];
%theta_0 = [1; 1];
domain_theta = [[-pi / 2, pi / 2]; [-pi / 2, pi / 2]; [-pi / 2, pi / 2]; [-pi / 2, pi / 2]];
G = @(theta) ode_integration(m_u, [m_e', m_u], [1, 1, theta'], R, 0);
% G(theta_0)
fprintf("Calcul des angles\n");
[theta, ~, ~, k] = SQP(theta_0, G, 1e-6, domain_theta, 30, 1, 1);

ode_integration(m_u, [m_e', m_u], [1, 1, theta'], R, 1)
