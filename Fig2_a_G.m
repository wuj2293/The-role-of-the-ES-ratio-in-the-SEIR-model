close all
clear all
%%%%% Figure2 %%%%%
%% Load Data
load('Guinea_originaldata.mat')
load('Guinea_week.mat')

ind_used_data = ind(1:1:260);
Cases_used_data = Cases(1:1:260);
Deaths_used_data = Deaths(1:1:260);

%% Curve Fitting
% function [fitresult, gof] = createFit(ind_used_data, Cases_used_data)
% [xData, yData] = prepareCurveData(ind_used_data, Cases_used_data);
% % Set up fittype and options.
% ft = fittype( '1/(b+c*exp(-a*x))', 'independent', 'x', 'dependent', 'y' );
% opts = fitoptions( 'Method', 'NonlinearLeastSquares' );
% opts.Algorithm = 'Levenberg-Marquardt';
% opts.Display = 'Off';
% opts.StartPoint = [0.0461713906311539 0.0971317812358475 0.823457828327293];
% % Fit model to data.
% [fitresult, gof] = fit( xData, yData, ft, opts );
Curvefit_cases = createFit(ind_used_data, Cases_used_data);
Curvefit_Deaths = createFit(ind_used_data, Deaths_used_data);

t = (1:1:700);
x = Curvefit_cases(t);
y = Curvefit_Deaths(t);

%% Plot
figure(1)
hold on
scatter(ind(1:201),Cases(1:201),30,'o');
scatter(ind(1:201),Deaths(1:201),30,'*');
f1 = plot(t,x,'--','linewidth',2);
f2 = plot(t,y,'linewidth',2);
hold off

xlim([0 670])
xticks([1 167.3333 333.6667 500 667])
xticklabels({'25 Mar 2014','07 Sep 2014','20 Feb 2015','06 Aug 2015', '20 Jan 2016'})
legend('Cases: original data','Deaths: original data','Cases: fitted curve',...
    'Deaths: fitted curve','Location','northwest')
title('Guinea')
ylabel('Cumulative cases')
