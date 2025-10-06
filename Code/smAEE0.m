function smAEE0
%%%%%%%%%%%%%%%%%%Figure1%%%%%%%%%%%%
%clear all
%clc
 global tau tau_2
a=0;
b=1000;
intervalo=[a,b];
    tau =0.01;
    tau_2=1;
    sol=dde23(@sistema,[tau,tau_2],@historiasistema,intervalo);
    sol1=dde23(@sistema,[tau,tau_2],@historiasistema1,intervalo);
    sol2=dde23(@sistema,[tau,tau_2],@historiasistema2,intervalo);
    sol3=dde23(@sistema,[tau,tau_2],@historiasistema3,intervalo);
%%%%Tamaño de letra del eje   
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
xlim([0 100]) 
h=legend({'$\varphi_{1}(\theta)$','$\varphi_{2}(\theta)$', '$\varphi_{3}(\theta)$','$\varphi_{4}(\theta)$'},'Interpreter', 'latex');   
h.FontSize = 14; % Set the legend font size directly  
        
figure(2)
hold on
plot(sol.x,sol.y(2,:),'b-','LineWidth', 1.5)
hold on
plot(sol1.x,sol1.y(2,:),'g:','LineWidth', 1.5)
hold on
plot(sol2.x,sol2.y(2,:),'r-.','LineWidth', 1.5)
hold on
plot(sol3.x,sol3.y(2,:),'k--','LineWidth', 1.5)
xlabel('t')
ylabel('y(t)')
ylim([0 90])
xlim([0 1.5])
h=legend({'$\varphi_{1}(\theta)$','$\varphi_{2}(\theta)$', '$\varphi_{3}(\theta)$','$\varphi_{4}(\theta)$'},'Interpreter', 'latex'); 
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
xlim([0 100]) 
h=legend({'$\varphi_{1}(\theta)$','$\varphi_{2}(\theta)$', '$\varphi_{3}(\theta)$','$\varphi_{4}(\theta)$'},'Interpreter', 'latex'); 
h.FontSize = 14; % Set the legend font size directly 

figure(4)
grid on
plot3(sol.y(1,:),sol.y(2,:),sol.y(3,:),'b','LineWidth', 1.5)
hold on
plot3(sol1.y(1,:),sol1.y(2,:),sol1.y(3,:),'g','LineWidth', 1.5)
hold on
plot3(sol2.y(1,:),sol2.y(2,:),sol2.y(3,:),'r','LineWidth', 1.5)
hold on
plot3(sol3.y(1,:),sol3.y(2,:),sol3.y(3,:),'k','LineWidth', 1.5)
%%%Punto E
plot3(1.0607e+03, 0, 0, 'ro', 'MarkerSize', 10, 'MarkerFaceColor', 'r');
text(1100, 3, 100, 'E0', 'FontSize', 12, 'FontWeight', 'bold', 'HorizontalAlignment','left');
%%%flechas%%%%%
%title('phase space')
xlabel('x(t)');
ylabel('y(t)');
zlabel('z(t)');
h=legend({'$\varphi_{1}(\theta)$','$\varphi_{2}(\theta)$', '$\varphi_{3}(\theta)$','$\varphi_{4}(\theta)$'},'Interpreter', 'latex'); 
h.FontSize = 14; % Set the legend font size directly 
axis square;grid on











function S= historiasistema(t)
S=[50;30;11];
function S1= historiasistema1(t)
S1=[110;40;15];
function S2= historiasistema2(t)
S2=[300;60;20];
function S3= historiasistema3(t)
S3=[500;80;30];


function dydt = sistema(t,y,z)
yretardo1=z(:,1);
yretardo2=z(:,2);
%%%Parametros adecuados
s=225;
d=0.1;
r=0.1;
beta=0.002;
q=0.5;
delta=0.001;
a=5;
b=0.1;
c=0.3;
m=0.1;
p=0.1;
xm=500;
tau=0.01;

x0=(xm/(2*r))*((r-d)+sqrt((r-d)^2+(4*r*s)/xm))
R0=(beta*exp(-delta*tau)*x0)/a

dydt=[s-d*y(1)+r*y(1)*(1-y(1)/xm)-(beta*y(1)*y(2))/(1+q*y(3))
      (beta*yretardo1(1)*yretardo1(2)*exp(-delta*tau))/(1+q*yretardo1(3))-a*y(2)-p*y(2)*y(3)
      c*yretardo2(2)*yretardo2(3)-b*y(3)-m*y(2)*y(3)];



   
