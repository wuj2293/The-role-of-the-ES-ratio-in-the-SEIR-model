close all
clear all

% Load Data
load('SierraLeone_originaldata.mat')
load('SierraLeone_week.mat')
ind_used_data = ind(1:1:150);
Cases_used_data = Cases(1:1:150);
Deaths_used_data = Deaths(1:1:150);

% Curve Fitting
Curvefit_cases = createFit(ind_used_data, Cases_used_data);
Curvefit_Deaths = createFit(ind_used_data, Deaths_used_data);

x_coeff = coeffvalues(Curvefit_cases);
a = x_coeff(1);
b = x_coeff(2);

y_coeff = coeffvalues(Curvefit_Deaths);
c = y_coeff(1);
d = y_coeff(2);

%% Construction Beta(t)
f = 0.304969310893238;
gamma = 0.098614519051379;
R0 = 0;

dx_dt = @(t) a*Curvefit_cases(t).*(1-b*Curvefit_cases(t));
dy_dt = @(t) c*Curvefit_Deaths(t).*(1-d*Curvefit_Deaths(t));

d2x_dt2 = @(t) (a^2)*Curvefit_cases(t).*(1-b*Curvefit_cases(t)).*(1-2*b*Curvefit_cases(t));
d2y_dt2 = @(t) (c^2)*Curvefit_Deaths(t).*(1-d*Curvefit_Deaths(t)).*(1-2*d*Curvefit_Deaths(t));

I = @(t) (1/(f*gamma))*dy_dt(t);

% sig_arr=[0.1887,0.0877,0.04762];
Nt = linspace(20000,80000,100); % N = 6092075;

for k = 1:100; 
    N = Nt(k);

    sigma_1 = 0.1887;
    E_1 = @(t) (1/sigma_1)*dx_dt(t);
    dE_1 = @(t) (1/sigma_1)*d2x_dt2(t);
    phi_1 = @(t) dE_1(t) + dx_dt(t);
    S0_1 = N-(E_1(1)+I(1)+R0);
    S_1 = @(t) S0_1 - integral(phi_1,0, t, 'ArrayValued', true);
    beta_1 = @(t) N*(1./(S_1(t).*I(t))).*phi_1(t);
    
    sigma_2 = 1/11.4;
    E_2 = @(t) (1/sigma_2)*dx_dt(t);
    dE_2 = @(t) (1/sigma_2)*d2x_dt2(t);
    phi_2 = @(t) dE_2(t) + dx_dt(t);
    S0_2 = N-(E_2(1)+I(1)+R0);
    S_2 = @(t) S0_2 - integral(phi_2,0, t, 'ArrayValued', true);
    beta_2 = @(t) N*(1./(S_2(t).*I(t))).*phi_2(t);
    
    sigma_3 = 0.04762;
    E_3 = @(t) (1/sigma_3)*dx_dt(t);
    dE_3 = @(t) (1/sigma_3)*d2x_dt2(t);
    phi_3 = @(t) dE_3(t) + dx_dt(t);
    S0_3 = N-(E_3(1)+I(1)+R0);
    S_3 = @(t) S0_3 - integral(phi_3,0, t, 'ArrayValued', true);
    beta_3 = @(t) N*(1./(S_3(t).*I(t))).*phi_3(t);
    
    root_fun = @(t) beta_1(t) - beta_3(t);
    x_root(k) = fzero(root_fun, 183);
end

%% Plot
close all
figure(1)

Nt = linspace(20000,80000,100);
% load('S_2mto20m.mat');

hold on
plot(Nt, x_root, 'LineWidth',2)
yline(189,'-','Rt = 1', 'LineWidth',2);
xlabel('S(0)')
ylabel('Days from May 27, 2014')
title('The day of Crosspoint at S(0) in SierraLeone','Interpreter', 'latex')
ylim([185 210])
% ylim([186 194])
hold off
