%% Add paths and load sym package if GNU Octave is used
clear
ddebiftoolpath='../../';
addpath(strcat(ddebiftoolpath,'ddebiftool'),...
strcat(ddebiftoolpath,'ddebiftool_extra_symbolic'));
if dde_isoctave()
pkg load symbolic
end
%% Create parameter names as strings and define fixed parameters
% The demo has the parameters | a| and |tau| parametros que varian
parnames ={'a','tau'};
s= sym(20,'r');%lambda
beta= sym(0.15,'r');
d=sym(0.35, 'r');%d1
%a=1.9;%d2
b=sym(0.45, 'r');%d3
c=sym(0.2, 'r');%gamma
p=sym(1.6, 'r');
r1=sym(0.01, 'r');%0.5;
xm=sym(150, 'r');%K
q=sym(0.05,'r');
m=sym(0.01,'r');


%% Create symbols for parameters , states and delays states
% | par | is the array of symbols in the same order as parnames .
% Due to the following two lines we may , for example ,
% use either beta or par (1) to refer to the delay .
syms(parnames{:}); % create symbols for beta , alpha and tua
par=cell2sym(parnames); % now beta is par (1) etc
%% Define system using symbolic algebra
% create symbols for u1(t) u1(t-tau), u2(t), u2(t-tau)
syms u1 u1t u2 u2t u3 u3t
du1_dt=s-d*u1+r1*u1*(1-u1/xm)-(beta*u1*u2)/(1+q*u3);
du2_dt=(beta*u1*u2)/(1+q*u3)-a*u2-p*u2*u3;
du3_dt=c*u2t*u3t-b*u3-m*u2*u3;

%% Differentiate and generate code (multi - linear forms )
[fstr, derivs ]=dde_sym2funcs(...
[du1_dt; du2_dt; du3_dt] ,... % n x 1 array of derivative symbolic expressions
[u1, u1t ; u2, u2t; u3, u3t] ,... % n x ( ntau +1) array of symbols for states ( current & delayed )
par ,... % 1 x np (or np x 1) array of symbols used for parameters
'filename','sym_canul_mf' ,... % optional argument specifying output file
'directional_derivative ',false ) ;
%% Differentiate and generate code ( directional derivatives )
[fstr, derivs]=dde_sym2funcs (...
[du1_dt; du2_dt; du3_dt] ,... % n x 1 array of derivative symbolic expressions
[u1, u1t ; u2, u2t; u3, u3t] ,... % n x ( ntau +1) array of symbols for states ( current &delayed )
par ,... % 1 x np (or np x 1) array of symbols used for parameters
'filename','sym_canul' ,... % optional argument specifying output file
'directional_derivative',true ) ;



