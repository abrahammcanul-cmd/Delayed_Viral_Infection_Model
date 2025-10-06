function modelo_original_tau1mod
%figure 3
% Plot local stability when a(T), c(T), b(T) > 0
 global retardos tau
%%%%%%%max 0.4
clc;

retardos=0.1;
tau=retardos
a=0;
b=2000;
c=-retardos;


d=0;
intervalo=[a,b];
sol=dde23(@sistema,retardos,@historiasistema,intervalo);
sol.y(1,end)
sol.y(2,end)
sol.y(3,end)

%%%%%%%%% Solution graph%%%%
set(0, 'DefaultAxesFontSize', 15);
set(0, 'DefaultTextFontSize', 15); 
figure(1)
w=linspace(c,d,1000);
S=historiasistema(w);
S1=historiasistema(w);
S2=historiasistema(w);
hold on
plot(sol.x,sol.y(1,:),'b','LineWidth',2);
xlabel('t');
ylabel('x(t)');

%%%%%phase space%%%%%%%
figure(2)
w=linspace(c,d,1000);
S=historiasistema(w);
S1=historiasistema(w);
S2=historiasistema(w);
hold on
plot(sol.x,sol.y(2,:),'g','LineWidth',2);
xlabel('t');
ylabel('y(t)');


%%%%%phase space%%%%%%%
figure(3)
w=linspace(c,d,1000);
S=historiasistema(w);
S1=historiasistema(w);
S2=historiasistema(w);
hold on
plot(sol.x,sol.y(3,:),'r','LineWidth',2)
xlabel('t');
ylabel('z(t)');


figure(4)
plot3(sol.y(1,:),sol.y(2,:),sol.y(3,:),'b','LineWidth',1.5)
hold on
plot3(887.9234,3,10.3,'*r','LineWidth',2)
xlabel('x(t)');
ylabel('y(t)');
zlabel('z(t)');
axis square; grid on







function S= historiasistema(t)
S=[868;11;10];%figure 3;


function S1= historiasistema1(t)
S1=[210;120;25];%figure 3;

function S2= historiasistema2(t)

S2=[330;200;30];%figure 3;



function dydt = sistema(t,y,z)
global tau
s=500; beta=0.003; d=0.5; a=1.9; b=0.003;
c=0.003; m=0.002; p=0.05; r=0.5; xm=800;
delta=0.0001;q=0.01; tau=0.1;
yretardo1=z(:,1);

x0=(xm/(2*r))*((r-d)+sqrt((r-d)^2+(4*r*s)/xm))
R0=(beta*exp(-delta*tau)*x0)/a
x1=a/(beta*exp(-delta*tau))
y1=(s*xm+(r-d)*xm*x1-r*x1^2)/(beta*x1*xm)
RCTL=((c-m)*y1)/b
CON=d-r+(2*r*887.9220)/xm




dydt=[s+r*y(1)*(1-y(1)/xm)-d*y(1)-(beta*y(1)*y(2))/(1+q*y(3))
      (beta*exp(-delta*tau)*yretardo1(1)*yretardo1(2))/(1+q*yretardo1(3))-a*y(2)-p*y(2)*y(3)
      (c-m)*y(2)*y(3)-b*y(3)];