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

sigma=1/11.4;

I = (1/(f*gamma))*dy_dt;
dI_dt=(1/(f*gamma))*d2y_dt2;
E = (1/sigma)*dx_dt;
dE_dt = (1/sigma)*d2x_dt2;
phi = dE_dt + sigma*E;

E0 = E(1);
I0 = I(1);
R0 = 0;

Nt = [20024 40024 60024 80024]; % linspace(20024,70024,6); %인구 범위
St = [20000 40000 60000 80000]; % linspace(20000,70000,6);

for i=1:length(Nt);
    N=Nt(i);
    S0 = N-(E0+I0+R0);
    S = S0 - cumtrapz(t,phi);
    sigind(i,:) = E./S;
end

%%
close all
figure(1)

hold on
S2m = plot(t,sigind(1,:),'k-','linewidth',2)
S4m = plot(t,sigind(2,:),'k--','linewidth',2);
S6m = plot(t,sigind(3,:),'k:','linewidth',2);
S8m = plot(t,sigind(4,:),'k-.','linewidth',2);

plot(find(sigind(1,:) == max(sigind(1,:))), max(sigind(1,:)), 'k*','linewidth',2);
text(find(sigind(1,:) == max(sigind(1,:))), max(sigind(1,:))+0.0004, ...
    ['(' num2str(find(sigind(1,:) == max(sigind(1,:)))) ' , ' num2str(max(sigind(1,:))) ')'])
plot(find(sigind(2,:) == max(sigind(2,:))), max(sigind(2,:)), 'k*','linewidth',2);
text(find(sigind(2,:) == max(sigind(2,:))), max(sigind(2,:))+0.0004, ...
    ['(' num2str(find(sigind(2,:) == max(sigind(2,:)))) ' , ' num2str(max(sigind(2,:))) ')'])
plot(find(sigind(3,:) == max(sigind(3,:))), max(sigind(3,:)), 'k*','linewidth',2);
text(find(sigind(3,:) == max(sigind(3,:))), max(sigind(3,:))+0.0004, ...
    ['(' num2str(find(sigind(3,:) == max(sigind(3,:)))) ' , ' num2str(max(sigind(3,:))) ')'])
plot(find(sigind(4,:) == max(sigind(4,:))), max(sigind(4,:)), 'k*','linewidth',2);
text(find(sigind(4,:) == max(sigind(4,:))), max(sigind(4,:))+0.0003, ...
    ['(' num2str(find(sigind(4,:) == max(sigind(4,:)))) ' , ' num2str(max(sigind(4,:))) ')'])

legend([S2m S4m S6m S8m],{'S(0)=20000','S(0)=40000', 'S(0)=60000', 'S(0)=80000'})
xlabel('Days from March 25, 2014')
ylabel('ratio')
title({'The ratio of $\frac{d}{dt}$(E/S)=0 in Guinea', '(Incubation period = 1/11.4)'}, 'Interpreter', 'latex');
hold off