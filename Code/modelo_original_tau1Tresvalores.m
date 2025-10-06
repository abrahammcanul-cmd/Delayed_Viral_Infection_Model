function modelo_original_tau1Tresvalores
%figure5
    global retardos tau
    clc;

    % Define the tau values to evaluate
    valores_tau = [0.1, 0.3, 2.4];

   % Create a main figure for all subplots
set(0, 'DefaultAxesFontSize', 15);
set(0, 'DefaultTextFontSize', 15); 
    figure;

    for i = 1:length(valores_tau)
        % Define general parameters
        retardos = valores_tau(i);
        a = 0;
        b = 600;
        intervalo = [a, b];
        tau = valores_tau(i);
        
        % Solve the system
        sol = dde23(@sistema, retardos, @historiasistema, intervalo);
        
       % Solution plots for each tau value in the same figure
        subplot(3, 3, (i - 1) * 3 + 1); % Solution plots for each tau value in the same figure
        plot(sol.x, sol.y(1,:), 'LineWidth', 1.5);
        xlabel('t');
        ylabel('x(t)');
        title([' \tau = ', num2str(tau)]);
        
        subplot(3, 3, (i - 1) * 3 + 2);
        plot(sol.x, sol.y(2,:), 'LineWidth', 1.5);
        xlabel('t');
        ylabel('y(t)');
        
        subplot(3, 3, (i - 1) * 3 + 3);
        plot(sol.x, sol.y(3,:), 'LineWidth', 1.5);
        xlabel('t');
        ylabel('z(t)');
    end
end

function S = historiasistema(t)
    S = [32.63; 13; 0.078];
end

function dydt = sistema(t, y, z)
    global tau
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

    yretardo1 = z(:, 1);

    dydt = [s + r * y(1) * (1 - y(1) / xm) - d * y(1) - (beta * y(1) * y(2)) / (1 + q * y(3));
            (beta * exp(-delta * tau) * yretardo1(1) * yretardo1(2)) / (1 + q * yretardo1(3)) - a * y(2) - p * y(2) * y(3);
            (c - m) * y(2) * y(3) - b * y(3)];
end


