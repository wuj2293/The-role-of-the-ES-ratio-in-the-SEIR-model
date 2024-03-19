clear
close all

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

Nt = linspace(20024,200024,200); % 1만부터 20만까지
St = linspace(20000,200000,200);

for i=1:length(Nt);
    N=Nt(i);
    S0 = N-(E0+I0+R0);
    S = S0 - cumtrapz(t,phi);
    sigind=E./S;
    test(i)=find(sigind==max(sigind));
end

%% Plot
figure(2)
hold on
plot(St,test,'LineWidth',2)
plot(St,237*ones(length(St),1),'LineWidth',2) % 하한선
hold off

title('Guinea')
ylim([236 243])
