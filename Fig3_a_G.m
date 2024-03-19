close all
clear all

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

A1 = [dx_dt -dy_dt];
B1 = d2y_dt2;
sol1 = A1\B1;
f = sol1(1)/sol1(2);
gamma = sol1(2);
gamma_1 = 1/gamma;

NN = length(t);
BETA = zeros(NN,3);
N = 11745189;
sig_arr=[0.1887,0.0877,0.04762];
for i=1:3
    sigma=sig_arr(i);
    I = (1/(f*gamma))*dy_dt;
    E = (1/sigma)*dx_dt;
    dE = (1/sigma)*d2x_dt2;
    phi = dE + dx_dt;
    E0 = E(1);
    I0 = I(1);
    D0 = y(1);
    R0 = 0;
    S0 = N-(E0+I0+R0);
    S = S0 - cumtrapz(t,phi);
    R = (1-f)*gamma.*cumtrapz(t,I);
    beta = N*(1./(S.*I)).*phi;
    for j=1:NN
        BETA(j,i) = beta(j);
    end
end
Cross_Point = BETA(237,1);

figure(1)
hold on
plot(t,BETA(:,1),'-.', t,BETA(:,2),'--',t,BETA(:,3),'linewidth',001);
plot([237 237],[0. 0.18],':r','linewidth',2);
hold off

title('Guinea')
legend('1/\sigma = 5.3','1/\sigma = 11.4','1/\sigma = 21')
text(244,Cross_Point,'CP: (237, 0.1076)')
xticks([1 167.3333 333.6667 500])
xticklabels({'25 Mar 2014','07 Sep 2014','20 Feb 2015','06 Aug 2015'})
ylabel('\beta (t)')
