function smAEqceroqmayorcero
%figure 11
%%Model Comparison
%clear all
%clc
 global tau tau_2
a=0;
b=1600;
intervalo=[a,b];
    tau =0.1;
    tau_2=0.7;
    sol=dde23(@sistema,[tau,tau_2],@historiasistema,intervalo);
    solsq=dde23(@sistemasq,[tau,tau_2],@historiasistema,intervalo);
   
   
set(0, 'DefaultAxesFontSize', 15);
set(0, 'DefaultTextFontSize', 15);     
figure(1)
hold on
plot(sol.x,sol.y(1,:),'b-','LineWidth', 1.5)
plot(solsq.x,solsq.y(1,:),'r-','LineWidth', 1.5)
legend({'q>0','q=0'},'Location', 'northeast')
xlabel('t')
ylabel('x(t)')
xlim([0 1500])    
  
        
figure(2)
hold on
plot(sol.x,sol.y(2,:),'b-','LineWidth', 1.5)
plot(solsq.x,solsq.y(2,:),'r-','LineWidth', 1.5)
legend({'q>0','q=0'},'Location', 'northeast')
xlabel('t')
ylabel('y(t)')
xlim([0 1500])  

figure(3)
hold on
plot(sol.x,sol.y(3,:),'b-','LineWidth', 1.5)
plot(solsq.x,solsq.y(3,:),'r-','LineWidth', 1.5)
legend({'q>0','q=0'},'Location', 'northeast')
xlabel('t')
ylabel('z(t)')
xlim([0 1500])  

Solucion= solve_nonlinear_system(tau)







function S= historiasistema(t)
S=[850;18;11];

function dydt = sistema(t,y,z)
yretardo1=z(:,1);
yretardo2=z(:,2);
%%%Define parameters
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

  function dydt = sistemasq(t,y,z)
yretardo1=z(:,1);
yretardo2=z(:,2);
%%%DEfine parameters
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
tau=0.1


dydt=[s-d*y(1)+r*y(1)*(1-y(1)/xm)-beta*y(1)*y(2)
      beta*yretardo1(1)*yretardo1(2)*exp(-delta*tau)-a*y(2)-p*y(2)*y(3)
      c*yretardo2(2)*yretardo2(3)-b*y(3)-m*y(2)*y(3)];

 


   
%Equilibrium Points Calculation
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

    