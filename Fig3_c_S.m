close all
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

% A1 = [dx_dt -dy_dt];
% pA1=pinv(A1,1e-2);
% B1 = d2y_dt2;
% sol1 = pA1*B1;
% f = sol1(1)/sol1(2);
% gamma = sol1(2);
% gamma_1 = 1/gamma;

A1 = [dx_dt -dy_dt];
B1 = d2y_dt2;
sol1 = A1\B1;
f = sol1(1)/sol1(2);
gamma = sol1(2);
gamma_1 = 1/gamma;

NN = length(t);
BETA = zeros(NN,3);
N = 6092075;
sig_arr=[0.1887,0.0877,0.04762];
for i=1:3
    %sigma = 0.05+i*0.05;
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
cr=BETA(186,1)

figure(1)
hold on
plot(t,BETA(:,1),'-.', t,BETA(:,2),'--',t,BETA(:,3),'linewidth',001);
plot([185 185],[0.0 0.18],':r','linewidth',2);
hold off

title('Sierra Leone')
legend('1/\sigma = 5.3','1/\sigma = 11.4','1/\sigma = 21')
text(187,0.105,'CP: (185, 0.1003)')
%text(244,Cross_Point,'CP: (237, 0.1076)')
xticks([1 167.3333 333.6667 500])
xticklabels({'27 May 2014','09 Nov 2014','24 Apr 2015','08 Oct 2015'})
ylabel('\beta (t)')
