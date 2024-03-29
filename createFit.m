function [fitresult, gof] = createFit11(ind_used_data, Cases_used_data)
%CREATEFIT3(IND_USED_DATA,CASES_USED_DATA)
%  Create a fit.
%
%  Data for 'untitled fit 1' fit:
%      X Input : ind_used_data
%      Y Output: Cases_used_data
%  Output:
%      fitresult : a fit object representing the fit.
%      gof : structure with goodness-of fit info.
%
%  See also FIT, CFIT, SFIT.

%  Auto-generated by MATLAB on 18-Dec-2019 10:23:52

%% Fit: 'untitled fit 1'.
[xData, yData] = prepareCurveData(ind_used_data, Cases_used_data);

% Set up fittype and options.
ft = fittype( '1/(b+c*exp(-a*x))', 'independent', 'x', 'dependent', 'y' );
opts = fitoptions( 'Method', 'NonlinearLeastSquares' );
opts.Algorithm = 'Levenberg-Marquardt';
opts.Display = 'Off';
opts.StartPoint = [0.0461713906311539 0.0971317812358475 0.823457828327293];

% Fit model to data.
[fitresult, gof] = fit( xData, yData, ft, opts );

% Plot fit with data.
% figure( 'Name', 'untitled fit 1' );
% h = plot( fitresult, xData, yData );
% legend( h, 'Cases_used_data vs. ind_used_data', 'untitled fit 1', 'Location', 'NorthEast' );
% % Label axes
% xlabel ind_used_data
% ylabel Cases_used_data
% grid on


