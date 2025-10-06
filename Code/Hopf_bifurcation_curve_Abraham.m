
%figure 6

% Define the range for a
a = linspace(0.9,2.3 , 100); % Avoid omega = 0 to prevent division by zero
%%%%%First branches
for j=1:length(a)
    T=a(j)
    f=solve_nonlinear_system(T)
    f
    f=solve_nonlinear_system(T)
    fdun= solve_nonlinear_system(1.9)
    [a2 a1 a0 b2 b1 b0]=GP(f(1),f(2),f(3))

    u=rt(f(1),f(2),f(3))
    w(1)=sqrt(u(1))
    w(2)=sqrt(u(2))
% Calculate sigma
sigma1 = b2*w(1).^4 + (a2*b1 - a1*b2 - b0).*w(1).^2 + a1*b0 - a0*b1;

% Calculate theta
theta1 = ((b1 - a2*b2).*w(1).^4 + (a0*b2 + a2*b0 - a1*b1).*w(1).^2 - a0*b0) ./ ((b2.*w(1).^2 - b0).^2 + (b1^2).*w(1).^2);
% Calculate sigma
sigma2 = b2*w(2).^4 + (a2*b1 - a1*b2 - b0).*w(2).^2 + a1*b0 - a0*b1;

% Calculate theta
theta2 = ((b1 - a2*b2).*w(2).^4 + (a0*b2 + a2*b0 - a1*b1).*w(2).^2 - a0*b0) ./ ((b2.*w(2).^2 - b0).^2 + (b1^2).*w(2).^2);
   if sigma1 >= 0
        tau0(j) = 1./w(1) * (2*0*pi + acos(theta1));
        tau1(j) = 1./w(1) * (2*1*pi + acos(theta1));
        tau2(j) = 1./w(1) * (2*2*pi + acos(theta1));
        tau3(j) = 1./w(1) * (2*3*pi + acos(theta1));
    else
        tau0(j) = 1./w(1) * (2*0*pi + 2*pi - acos(theta1));
        tau1(j) = 1./w(1) * (2*1*pi + 2*pi - acos(theta1));
        tau2(j) = 1./w(1) * (2*2*pi + 2*pi - acos(theta1));
        tau3(j) = 1./w(1) * (2*3*pi + 2*pi - acos(theta1));
   end
    if sigma2 >= 0
        tau20(j) = 1./w(2) * (2*0*pi + acos(theta2));
        tau21(j) = 1./w(2) * (2*1*pi + acos(theta2));
        tau22(j) = 1./w(2) * (2*2*pi + acos(theta2));
        tau23(j) = 1./w(2) * (2*3*pi + acos(theta2));
    else
        tau20(j) = 1./w(2) * (2*0*pi + 2*pi - acos(theta2));
        tau21(j) = 1./w(2) * (2*1*pi + 2*pi - acos(theta2));
        tau22(j) = 1./w(2) * (2*2*pi + 2*pi - acos(theta2));
        tau23(j) = 1./w(2) * (2*3*pi + 2*pi - acos(theta2));
   end
end
tau0(1)
tau1(1)
%a
%tau
%Plot the piecewise function
figure;
hold on
plot(a, tau0,'b-', 'LineWidth', 1.5);
plot(a, tau22,'r--', 'LineWidth', 1.5);
plot(a, tau1,'b-', 'LineWidth', 1.5);
plot(a, tau2,'b-', 'LineWidth', 1.5);
plot(a, tau3,'b-', 'LineWidth', 1.5);
plot(a, tau20,'r--', 'LineWidth', 1.5);
plot(a, tau21,'r--', 'LineWidth', 1.5);
plot(a, tau23,'r--', 'LineWidth', 1.5);
plot(1.557,4.712,'k*', 'LineWidth', 1.5);
text(1.5,3.9, 'HH_1 (1.557,4.712)', 'FontSize', 12, 'FontWeight', 'bold', 'HorizontalAlignment','left');
plot(2.157,15.558,'k*', 'LineWidth', 1.5);
text(2,14, 'HH_4 (2.157,15.558)', 'FontSize', 12, 'FontWeight', 'bold', 'HorizontalAlignment','left');
plot(2.016,10.057,'k*', 'LineWidth', 1.5);
text(1.9,8.5, 'HH_3 (2.016,10.057)', 'FontSize', 12, 'FontWeight', 'bold', 'HorizontalAlignment','left');
plot(1.640,13.545,'k*', 'LineWidth', 1.5);
text(1.65,14.5, 'HH_2 (1.640,13.545)', 'FontSize', 12, 'FontWeight', 'bold', 'HorizontalAlignment','left');
 xlim([1.07 2.3])
 ylim([0 16])
xlabel('$a$');
ylabel('$\tau_2$');
box on;
grid on;
legend('\tau^{1}_2_{n}','\tau^{2}_2_{n}', 'Location', 'northwest')


function r1=rt(aa,bb,cc)
x=aa;
y=bb;
z=cc;
%%%%%%%%%%%Parameters%%%%%%%%%%%%%%%%%%%
s=20;
beta=0.15;
d=0.35;
b=0.45;
c=0.2;
p=1.6;
r=0.01;
xm=150;
q=0.05;
m=0.01;

a2=b+(s/x)+(r*x/xm)+m*y;
a1=(b+m*y)*((s/x)+(r*x/xm))+((beta^2)*x*y)/(1+q*z)^2 -p*m*y*z+ (beta*m*q*x*y*z)/(1+q*z)^2;
a0=((beta^2)*x*y*(b+m*y))/(1+q*z)^2+((beta^2)*x*y^2*q*m*z)/(1+q*z)^3 -((s/x)+(r*x/xm))*(p*m*y*z+(beta*m*q*x*y*z)/(1+q*z)^2);
b2=-(b+m*y);
b1=p*(b+m*y)*z-(b+m*y)*((s/x)+(r*x/xm))+(beta*q*x*z*(b+m*y))/(1+q*z)^2;
b0=(p*(b+m*y)*z+(beta*q*x*z*(b+m*y))/(1+q*z)^2)*((s/x)+(r*x/xm))-((beta^2)*(b+m*y)*x*y)/(1+q*z)^2-((beta^2)*q*(b+m*y)*x*y*z)/(1+q*z)^3;

C1=a2^2-2*a1-b2^2
C2=a1^2-2*a0*a2-b1^2+2*b0*b2
C3=a0^2-b0^2
PD3=[1 C1 C2 C3];

r11=roots(PD3)
    valid_r= find(imag(r11) == 0 & r11 > 0);
    
    % Store positive real solutions in a matrix
    r1 =r11(valid_r)
end

function [a2 a1 a0 b2 b1 b0]=GP(aa,bb,cc)
x=aa;
y=bb;
z=cc;
%%%%%%%%%%%Parameters%%%%%%%%%%%%%%%%%%%
s=20;
beta=0.15;
d=0.35;
b=0.45;
c=0.2;
p=1.6;
r=0.01;
xm=150;
q=0.05;
m=0.01;

a2=b+(s/x)+(r*x/xm)+m*y;
a1=(b+m*y)*((s/x)+(r*x/xm))+((beta^2)*x*y)/(1+q*z)^2 -p*m*y*z+ (beta*m*q*x*y*z)/(1+q*z)^2;
a0=((beta^2)*x*y*(b+m*y))/(1+q*z)^2+((beta^2)*x*y^2*q*m*z)/(1+q*z)^3 +((s/x)+(r*x/xm))*(p*m*y*z+(beta*m*q*x*y*z)/(1+q*z)^2);
b2=-(b+m*y);
b1=p*(b+m*y)*z-(b+m*y)*((s/x)+(r*x/xm))+(beta*q*x*z*(b+m*y))/(1+q*z)^2;
b0=(p*(b+m*y)*z+(beta*q*x*z*(b+m*y))/(1+q*z)^2)*((s/x)+(r*x/xm))-((beta^2)*(b+m*y)*x*y)/(1+q*z)^2-((beta^2)*q*(b+m*y)*x*y*z)/(1+q*z)^3;
end

function positive_solutions = solve_nonlinear_system(T)
    % Define symbolic variables
    syms x y v
    
    % Define parameters
a = T;
s=20;
beta=0.15;
d=0.35;
b=0.45;
c=0.2;
p=1.6;
r=0.01;
xm=150;
q=0.05;
m=0.01;



    % Define the equations
    F1 = s+r*x*(1-x/xm)-d*x-(beta*x*y)/(1+q*v) == 0;
    F2 = (beta*x*y)/(1+q*v) -a*y-p*y*v == 0;
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
end