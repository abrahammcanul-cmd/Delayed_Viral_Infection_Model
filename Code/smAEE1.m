function smAEE1
%figure2
%clear all
%clc
 global tau tau_2
a=0;
b=2000;
intervalo=[a,b];
    tau =0.01;
    tau_2=2.9;
    sol=dde23(@sistema,[tau,tau_2],@historiasistema,intervalo);
    sol1=dde23(@sistema,[tau,tau_2],@historiasistema1,intervalo);
    sol2=dde23(@sistema,[tau,tau_2],@historiasistema2,intervalo);
    sol3=dde23(@sistema,[tau,tau_2],@historiasistema3,intervalo);
   
%%%%Axis font size   
set(0, 'DefaultAxesFontSize', 15);
set(0, 'DefaultTextFontSize', 15);  
%%%%%%%%%%%%%%   
figure(1)
hold on
plot(sol.x,sol.y(1,:),'b-','LineWidth', 1.5)
hold on
plot(sol1.x,sol1.y(1,:),'g:','LineWidth', 1.5)
hold on
plot(sol2.x,sol2.y(1,:),'r-.','LineWidth', 1.5)
hold on
plot(sol3.x,sol3.y(1,:),'k--','LineWidth', 1.5)
xlabel('t')
ylabel('x(t)')
xlim([0 350])    
h=legend({'$\tilde{\varphi}_{1}(\theta)$','$\tilde{\varphi}_{2}(\theta)$', '$\tilde{\varphi}_{3}(\theta)$','$\tilde{\varphi}_{4}(\theta)$'},'Interpreter', 'latex');   
h.FontSize = 14; % Set the legend font size directly   
        
figure(2)
hold on
plot(sol.x,sol.y(2,:),'b-','LineWidth', 1.5)
hold on
plot(sol1.x,sol1.y(2,:),'g:','LineWidth', 1.5)
hold on
plot(sol2.x,sol2.y(2,:),'r-..','LineWidth', 1.5)
hold on
plot(sol3.x,sol3.y(2,:),'k--','LineWidth', 1.5)
xlabel('t')
ylabel('y(t)')
xlim([0 350])
h=legend({'$\tilde{\varphi}_{1}(\theta)$','$\tilde{\varphi}_{2}(\theta)$', '$\tilde{\varphi}_{3}(\theta)$','$\tilde{\varphi}_{4}(\theta)$'},'Interpreter', 'latex');   
h.FontSize = 14; % Set the legend font size directly 
figure(3)
hold on
plot(sol.x,sol.y(3,:),'b-','LineWidth', 1.5)
hold on
plot(sol1.x,sol1.y(3,:),'g:','LineWidth', 1.5)
hold on
plot(sol2.x,sol2.y(3,:),'r-.','LineWidth', 1.5)
hold on
plot(sol3.x,sol3.y(3,:),'k--','LineWidth', 1.5)
xlabel('t')
ylabel('z(t)')
xlim([0 1200])
h=legend({'$\tilde{\varphi}_{1}(\theta)$','$\tilde{\varphi}_{2}(\theta)$', '$\tilde{\varphi}_{3}(\theta)$','$\tilde{\varphi}_{4}(\theta)$'},'Interpreter', 'latex');   
h.FontSize = 14; % Set the legend font size directly 

figure(4)
%Parameters to estimate x0
s=255;
d=0.1;
r=0.1;
xm=500;
tau=0.01;

x0=(xm/(2*r))*((r-d)+sqrt((r-d)^2+(4*r*s)/xm))
%%%%%%%%%%%%%%%%%%%
plot3(sol.y(1,:),sol.y(2,:),sol.y(3,:),'b','LineWidth', 1.5)
hold on
plot3(sol1.y(1,:),sol1.y(2,:),sol1.y(3,:),'g','LineWidth', 1.5)
hold on
plot3(sol2.y(1,:),sol2.y(2,:),sol2.y(3,:),'r','LineWidth', 1.5)
hold on
plot3(sol3.y(1,:),sol3.y(2,:),sol3.y(3,:),'k','LineWidth', 1.5)
hold on
plot3(x0,0,0,'ro', 'MarkerSize', 10, 'MarkerFaceColor', 'r');
xlabel('x(t)');
ylabel('y(t)');
zlabel('z(t)');
h=legend({'$\tilde{\varphi}_{1}(\theta)$','$\tilde{\varphi}_{2}(\theta)$', '$\tilde{\varphi}_{3}(\theta)$','$\tilde{\varphi}_{4}(\theta)$'},'Interpreter', 'latex');   
h.FontSize = 14; % Set the legend font size directly 
axis square;grid on




 function S= historiasistema(t)
S=[110;40;15];
function S1= historiasistema1(t)
S1=[300;60;20];
function S2= historiasistema2(t)
S2=[500;80;30];
 function S3= historiasistema3(t)
 S3=[1000;50;0.1];


function dydt = sistema(t,y,z)
yretardo1=z(:,1);
yretardo2=z(:,2);
%%%Parameters
s=255;
d=0.1;
r=0.1;
beta=0.02;
q=0.001;
delta=0.001;
a=5;
b=0.1;
c=0.003;
m=0.001;
p=0.1;
xm=500;
tau=0.01;

x0=(xm/(2*r))*((r-d)+sqrt((r-d)^2+(4*r*s)/xm))
R0=(beta*exp(-delta*tau)*x0)/a
x1=a/(beta*exp(-delta*tau))
y1=(s*xm+(r-d)*xm*x1-r*x1^2)/(beta*x1*xm)
RCTL=((c-m)*y1)/b
Con=d-r+((r*x1)/xm)

dydt=[s-d*y(1)+r*y(1)*(1-y(1)/xm)-(beta*y(1)*y(2))/(1+q*y(3))
      (beta*yretardo1(1)*yretardo1(2)*exp(-delta*tau))/(1+q*yretardo1(3))-a*y(2)-p*y(2)*y(3)
      c*yretardo2(2)*yretardo2(3)-b*y(3)-m*y(2)*y(3)];



   
%Equilibrium Point Computation
function positive_solutions = solve_nonlinear_system(T)
    % Define symbolic variables
    syms x y v
 global s, global beta,global d, global a, global b, global c, global m
global p, global r, global xm, global delta, global q  

s=255;
d=0.1;
r=0.1;
beta=0.02;
q=0.5;
delta=0.001;
a=5;
b=0.1;
c=0.003;
m=0.001;
p=0.1;
xm=500;
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

    
