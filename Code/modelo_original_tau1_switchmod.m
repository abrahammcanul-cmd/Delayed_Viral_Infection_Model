function modelo_original_tau1_switchmod
%figure 4
%format long g

% Define the tau values for the two step intervals
tau_value = [0.1:0.01:2.2, 2.3:0.1:18];%[2.3:0.1:18] 

global tau
for i = 1:length(tau_value)
    tau_value(i);
    G = solve_nonlinear_system(tau_value(i));
    g1(i) = G(1);
    g2(i) = G(2);
    g3(i) = G(3);
end

for i = 1:length(tau_value)
    tau = tau_value(i);
    x0 = [g1(i) + 1, g2(i) + 1, g3(i) + 1];
    intervalo = [0, 1200];
    sol = dde23(@sistema, tau, x0, intervalo); 
    
    % Last 10 points for analysis
    last_points = sol.y(:, end-10:end);
    
    
    plot3(tau_value(i), g1(i), g2(i), '*r', 'LineWidth', 2);
    plot3(sol.y(2, end-100:end) .* 0 + tau_value(i), sol.y(1, end-100:end), sol.y(2, end-100:end), 'b');
    grid on;
    hold on;
end

xlabel('\tau_1', 'Fontsize', 10);
ylabel('x(t)', 'Fontsize', 10);
zlabel('y(t)', 'Fontsize', 10);

end

function dydt = sistema(t, y, z)
global tau

% Parameters
s = 10;
beta = 0.023;
d = 0.1;
a = 0.08;
b = 0.75;
c = 0.04;
m = 0.01;
p = 9;
r = 0.6;
xm = 500;
delta = 0.06;
q = 0.001;

yretardo1 = z(:,1);
dydt = [
    s + r * y(1) * (1 - y(1) / xm) - d * y(1) - (beta * y(1) * y(2)) / (1 + q * y(3));
    (beta * exp(-delta * tau) * yretardo1(1) * yretardo1(2)) / (1 + q * yretardo1(3)) - a * y(2) - p * y(2) * y(3);
    (c - m) * y(2) * y(3) - b * y(3)
];
end

function positive_solutions = solve_nonlinear_system(T)
syms x y v

% Parameters
s = 10;
beta = 0.023;
d = 0.1;
a = 0.08;
b = 0.75;
c = 0.04;
m = 0.01;
p = 9;
r = 0.6;
xm = 500;
delta = 0.06;
q = 0.001;
Tau = T;

% System Equations
F1 = s - d * x + r * x * (1 - x / xm) - (beta * x * y) / (1 + q * v) == 0;
F2 = (beta * exp(-delta * Tau) * x * y) / (1 + q * v) - a * y - p * y * v == 0;
F3 = (c - m) * y * v - b * v == 0;

% Symbolic Solution of the System
solutions = solve([F1, F2, F3], [x, y, v]);
x_sol = double(solutions.x);
y_sol = double(solutions.y);
v_sol = double(solutions.v);

% Filter Positive Solutions
tolerance = 1e-6;
positive_solutions = [];
for i = 1:length(x_sol)
    if real(x_sol(i)) > 0 && abs(imag(x_sol(i))) < tolerance && ...
       real(y_sol(i)) > 0 && abs(imag(y_sol(i))) < tolerance && ...
       real(v_sol(i)) > 0 && abs(imag(v_sol(i))) < tolerance
        positive_solutions = [positive_solutions; real(x_sol(i)), real(y_sol(i)), real(v_sol(i))];
    end
end

% Select the First Positive Solution
positive_solutions = [positive_solutions(1), positive_solutions(2), positive_solutions(3)];
end
