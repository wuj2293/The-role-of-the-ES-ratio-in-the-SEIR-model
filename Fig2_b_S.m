close all
clear all
%%%%% Figure2 %%%%%
%% Load Data
load('SierraLeone_originaldata.mat')
load('SierraLeone_week.mat')

ind_used_data = ind(1:1:150);
Cases_used_data = Cases(1:1:150);
Deaths_used_data = Deaths(1:1:150);

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

t = (1:1:500);
x = Curvefit_cases(t);
y = Curvefit_Deaths(t);

%% Plot
figure(1)
hold on
scatter(ind(1:151),Cases(1:151),30,'o');
scatter(ind(1:151),Deaths(1:151),30,'*');
% scatter(ind(:),Cases(:),30,'o');
% scatter(ind(:),Deaths(:),30,'*');
f1 = plot(t,x,'--','linewidth',2);
f2 = plot(t,y,'linewidth',2);
hold off

xlim([0 520])
ylim([0 14200])
xticks([1 167.3333 333.6667 500])
xticklabels({'27 May 2014','09 Nov 2014','24 Apr 2015','08 Oct 2015'})
legend('Cases: original data','Deaths: original data','Cases: fitted curve',...
    'Deaths: fitted curve','Location','northwest')
title('Sierra Leone')
ylabel('Cumulative cases')
