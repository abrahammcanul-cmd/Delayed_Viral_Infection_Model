function simulacionvaria_q1
%figure 10

    zg =[887; 2; 9];
    q=0:0.01:0.12
    for j=1:length(q)
        f = fsolve(@(z) myfunctionq(z, 0.1,q(j)), zg)
        f1 = fsolve(@(z) myfunctionp(z, 0.1,q(j)), zg)
        X(j)=f(1)
        X1(j)=f1(1)
        Z(j)=f(3)
        Z1(j)=f1(3)
    end
    set(0, 'DefaultAxesFontSize', 15);
    set(0, 'DefaultTextFontSize', 15); 
    figure(1)
    hold on
    plot(q,X,'-.b','LineWidth', 1.5)
    plot(q,X1,'-.r','LineWidth', 1.5)
    ylabel('x_2')
    legend({'q','p'},'Location','southwest')
    
    figure(2)
    hold on
    plot(q,Z,'-.b','LineWidth', 1.5)
    plot(q,Z1,'-.r','LineWidth', 1.5)
    ylabel('z_2')
    legend({'q','p'},'Location','southwest')
end

function [F] = myfunctionq(z, T,q1)
    tau = T;
    q=q1;
s=500;
beta=0.003;
d=0.5;
a=1.9;
b=0.003;
c=0.003;
m=0.002;
p=0.05;
r=0.5;
K=800;
delta=0.0001;

    x = z(1);
    y = z(2);
    v = z(3);

    % Define the equations as residuals
    F(1) = s - d * x + r * x * (1 - x / K) - beta * x * y / (1 + q1 * v);
    F(2) = (beta * exp(-delta * tau) * x * y) / (1 + q1 * v) - a * y - p * y * v;
    F(3) = (c - m) * y * v - b * v;
end
function [F1] = myfunctionp(z, T,p)
    tau = T;
    p1=p;
s=500;
beta=0.003;
d=0.5;
a=1.9;
b=0.003;
c=0.003;
m=0.002;
r=0.5;
K=800;
delta=0.0001;
q=0.01;

    x = z(1);
    y = z(2);
    v = z(3);

    % Define the equations as residuals
    F1(1) = s - d * x + r * x * (1 - x / K) - beta * x * y / (1 + q * v);
    F1(2) = (beta * exp(-delta * tau) * x * y) / (1 + q * v) - a * y - p1 * y * v;
    F1(3) = (c - m) * y * v - b * v;
end

