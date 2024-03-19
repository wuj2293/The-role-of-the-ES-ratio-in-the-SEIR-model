close all
clear all
%%%%% Figure3 %%%%%
%% Load Data
load('Guinea_originaldata.mat')
load('Guinea_week.mat')
ind_used_data = ind(1:1:260);
Cases_used_data = Cases(1:1:260);
Deaths_used_data = Deaths(1:1:260);

%% Curve Fitting
Curvefit_cases = createFit(ind_used_data, Cases_used_data);
Curvefit_Deaths = createFit(ind_used_data, Deaths_used_data);
t = (1:1:700);
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

N = 11745189;
sigma=1/11.4;

A1 = [dx_dt -dy_dt];
B1 = d2y_dt2;
sol1 = A1\B1;
f = sol1(1)/sol1(2);
gamma = sol1(2);
gamma_1 = 1/gamma;

I = (1/(f*gamma))*dy_dt;
dI_dt=(1/(f*gamma))*d2y_dt2;
E = (1/sigma)*dx_dt;
dE_dt = (1/sigma)*d2x_dt2;
phi = dE_dt + dx_dt+5;

E0 = E(1);
I0 = I(1);
R0 = 0;
S0 = N-(E0+I0+R0);
S = S0 - cumtrapz(t,phi);
R = (1-f)*gamma.*cumtrapz(t,I);

beta = N*(dE_dt+sigma*E)./(S.*I); 
Re = beta/gamma;

%% Plot
figure(3)

subplot(2,1,1)
hold on
plot(t, beta,'linewidth',2);
plot([237 237],[0.0 1],':r','linewidth',2);
hold off

title('Guinea')
legend('1/\sigma=11.4','Location','northeast')
text(238,0.11,'CP: (237, 0.1076)')
ylim([0.01 0.17])
ylabel('\beta (t)')
xticks([1 167.3333 333.6667 500])
xticklabels({'25 Mar 2014','07 Sep 2014','20 Feb 2015','06 Aug 2015'})

subplot(2,1,2)
hold on
plot(t, Re,'linewidth',2);
plot([241 241],[0.0 1],':r','linewidth',2);
plot([0 241],[1 1],':r','linewidth',2);
hold off

title('Guinea')
legend('1/\sigma=11.4','Location','northeast')
text(242,1.05,'Rt: (241, 1)')
ylabel('Rt')
xticks([1 167.3333 333.6667 500])
xticklabels({'25 Mar 2014','07 Sep 2014','20 Feb 2015','06 Aug 2015'})
