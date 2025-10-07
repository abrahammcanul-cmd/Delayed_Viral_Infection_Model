function smAEEEE
%figure 8
%0.1
%clear all
%clc
 global tau tau_2
a=0;
b=2000;
%c=800;
%intervalo1=linspace(c,b,200);
intervalo=[a,b];
    tau =0.1;
    tau_2=0.1;
    sol=dde23(@sistema,[tau,tau_2],@historiasistema,intervalo);
   
    
figure(1)
plot(sol.x,sol.y(1,:),'b-','LineWidth', 1.5)
xlabel('t')
ylabel('x(t)')
%ylabel('x(t), y(t), z(t)')
%xlim([0 10])    
  
        
figure(2)
plot(sol.x,sol.y(2,:),'g-','LineWidth', 1.5)
xlabel('t')
ylabel('y(t)')


figure(3)
plot(sol.x,sol.y(3,:),'r-','LineWidth', 1.5)
xlabel('t')
ylabel('z(t)')


figure(4)
plot3(sol.y(1,:),sol.y(2,:),sol.y(3,:),'b','LineWidth', 1.5)
title('phase space')
xlabel('x(t)');
ylabel('y(t)');
zlabel('z(t)');
axis square;grid on










function S= historiasistema(t)
S=[850;18;11];

function dydt = sistema(t,y,z)
yretardo1=z(:,1);
yretardo2=z(:,2);

%%%Parameters
s=500;
beta=0.003;
d=0.5;
a=1.9;
b=0.003;
c=0.003;
m=0.002;
p=0.05;
r=0.5;
xm=800;
delta=0.0001;
q=0.01;
tau=0.1;

dydt=[s-d*y(1)+r*y(1)*(1-y(1)/xm)-(beta*y(1)*y(2))/(1+q*y(3))
      (beta*yretardo1(1)*yretardo1(2)*exp(-delta*tau))/(1+q*yretardo1(3))-a*y(2)-p*y(2)*y(3)
      c*yretardo2(2)*yretardo2(3)-b*y(3)-m*y(2)*y(3)];



   
%Steady-State Calculation
function positive_solutions = solve_nonlinear_system(T)
    % Define symbolic variables
    syms x y v
 global s, global beta,global d, global a, global b, global c, global m
global p, global r, global xm, global delta, global q  
%     % Define parameters

s=500;
beta=0.003;
d=0.5;
a=1.9;
b=0.003;
c=0.003;
m=0.002;
p=0.05;
r=0.5;
xm=800;
delta=0.0001;
q=0.01;
tau=0.1;
Tau=T;


    % Define the equations
    F1 = s-d*x+r*x*(1-x/xm)-(beta*x*y)/(1+q*v) == 0;
    F2 = (beta*exp(-delta*Tau)*x*y)/(1+q*v)-a*y-p*y*v == 0;
    F3 = (c-m)*y*v-b*v == 0;
    
    % Solve the system of equations symbolically
    solutions = solve([F1, F2, F3], [x, y, v]);
    
    % Extract numeric solutions
    x_sol = double(solutions.x);
    y_sol = double(solutions.y);
    v_sol = double(solutions.v);
    
    % Filter positive real solutions
    valid_indices = find(imag(x_sol) == 0 & x_sol > 0 & imag(y_sol) == 0 & y_sol > 0 & imag(v_sol) == 0 & v_sol > 0);
    
    % Store positive real solutions in a matrix
    positive_solutions = [x_sol(valid_indices), y_sol(valid_indices), v_sol(valid_indices)];

    

