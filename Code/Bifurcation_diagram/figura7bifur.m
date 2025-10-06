%% Load DDE-BifTool into MATLAB path
clear
close all
ddebiftoolpath='../../';
addpath(strcat(ddebiftoolpath,'ddebiftool'),...
    strcat(ddebiftoolpath,'ddebiftool_extra_psol'),...
    strcat(ddebiftoolpath,'ddebiftool_extra_nmfm'),...
    strcat(ddebiftoolpath,'ddebiftool_utilities'));

format compact
set(groot,'defaultTextInterpreter','LaTeX');

%% Initial parameters and state
parnames={'a','tau'};
cind=[parnames;num2cell(1:length(parnames))];
ind=struct(cind{:});
bounds = {'max_bound',[ind.tau 17],'max_step',[0,0.5],...
    'min_bound',[ind.tau 0]};

%% Set user-defined functions
funcs=set_symfuncs(@sym_canul_mf,'sys_tau',@()ind.tau,'x_vectorized',true,'p_vectorized',true);

%% Construct steady-state point
stst = dde_stst_create('x',[29.6816;2.3684;1.4117]);
stst.parameter(ind.tau)=0;
stst.parameter(ind.a)=1.9;

%% Hopf - Hopf point
stst.kind = 'stst';
stst.x =[29.6816;2.3684;1.4117];
stst.parameter(ind.tau)=0;
stst.parameter(ind.a)=1.9;
method = df_mthod (funcs,'stst');
method_stst.stability.minimal_real_part=-25;
stst.stability = p_stabil(funcs,stst,method.stability);

%% Initialization of branch of steady-states
contpar=ind.tau;
steadystate_br=SetupStst(funcs,'x',stst.x,'parameter',stst.parameter,...
    'step',0.1,'contpar',contpar,'max_step',[contpar,0.3],bounds{:});

%% Continue steady-state branch
n_steps=200;
steadystate_br=br_contn(funcs,steadystate_br,n_steps);

%% Detect bifurcations on steady-state branch
[steadystate_br,~,ind_hopf,bif1types]=LocateSpecialPoints(funcs,...
    steadystate_br);
nunst_eqs=GetStability(steadystate_br);
fprintf('Hopf bifurcation near point %d\n',ind_hopf);

%% Initialize limit-cycle branches from Hopf bifurcations
psol_branches = {};
for i = 1:length(ind_hopf)
    [psol_br,success] = SetupPsol(funcs,steadystate_br,ind_hopf(i),...
        'contpar',contpar,'degree',4,'intervals',40,'max_bound',[ind.tau,16],...
        'max_step',[0,0.5],'print_residual_info',1);
    if success
        psol_branches{end+1} = br_contn(funcs,psol_br,100);
    end
end

%% Plot the amplitude y(t) vs tau - FIGURE 3
figure(3); clf; hold on;
for i = 1:length(psol_branches)
    [amplitude_y, param_tau] = get_cycle_amplitude(psol_branches{i}, ind.tau);
    plot(param_tau, amplitude_y, '-b', 'LineWidth', 2, 'DisplayName', sprintf('Branch %d', i));
end

xlabel('$\tau$', 'FontSize', 12);
ylabel('Amplitude of $y(t)$', 'FontSize', 12);
legend('show');
title('Amplitude of $y(t)$ vs $\tau$', 'FontSize', 14);
grid on;

%% Function to compute the amplitude of y(t) for limit-cycle branches
function [amplitude_y, param_tau] = get_cycle_amplitude(branch, tau_ind)
    num_points = length(branch.point);
    amplitude_y = zeros(1, num_points);
    param_tau = zeros(1, num_points);
    for j = 1:num_points
        cycle = branch.point(j);
        param_tau(j) = cycle.parameter(tau_ind);
        y = cycle.profile(2,:); % Assuming y is the second component
        amplitude_y(j) = max(y) - min(y); % Compute amplitude based on y(t)
    end
end