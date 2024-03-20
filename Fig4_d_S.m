clear all

%% Load Data
load('SierraLeone_originaldata.mat')
load('SierraLeone_week.mat')
ind_used_data = ind(1:1:150);
Cases_used_data = Cases(1:1:150);
Deaths_used_data = Deaths(1:1:150);

%% Curve Fitting
Curvefit_cases = createFit(ind_used_data, Cases_used_data);
Curvefit_Deaths = createFit(ind_used_data, Deaths_used_data);
t = (1:1:500);
x = Curvefit_cases(t);
y = Curvefit_Deaths(t);

%% Construction Beta(t)

x_coeff = coeffvalues(Curvefit_cases);
a = x_coeff(1);
b = x_coeff(2);

y_coeff = coeffvalues(Curvefit_Deaths);
c = y_coeff(1);
d = y_coeff(2);

dx_dt = a*x.*(1-b*x);
dy_dt = c*y.*(1-d*y);

d2x_dt2 = (a^2)*x.*(1-b*x).*(1-2*b*x);
d2y_dt2 = (c^2)*y.*(1-d*y).*(1-2*d*y);

A1 = [dx_dt -dy_dt];
pA1=pinv(A1,1e-2);
B1 = d2y_dt2;
sol1 = pA1*B1;
f = sol1(1)/sol1(2);
gamma = sol1(2);
gamma_1 = 1/gamma;

sigma=1/11.4;

I = (1/(f*gamma))*dy_dt;
dI_dt=(1/(f*gamma))*d2y_dt2;
E = (1/sigma)*dx_dt;
dE_dt = (1/sigma)*d2x_dt2;
phi = dE_dt + dx_dt+5;

E0 = E(1);
I0 = I(1);
R0 = 0;

% St = [20000 30000 40000 50000 60000 70000 80000];
% St = [40000 50000 60000 70000 80000 90000 100000 110000 120000];
St = 2e4:1e4:20e4;
% St = linspace(20000,200000,200);

for i=1:length(St);
    N=St(i)+74;
    S0 = N-(E0+I0+R0);
    S = S0 - cumtrapz(t,phi);
    sigind=E./S;
    test(i)=find(sigind==max(sigind));
end

%% Plot
close all
figure(1)

hold on
plot(St,test,'k.', 'MarkerSize',20)
yline(189,'-','Rt = 1', 'LineWidth',2);
xlabel('S(0)')
ylabel('Days from May 27, 2014')
title({'The extreme point($\frac{d}{dt}$(E/S)=0) in SierraLeone', '(Incubation period = 1/11.4)'}, 'Interpreter', 'latex');
ylim([185 210])
% ylim([186 194])
hold off
