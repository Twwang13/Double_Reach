function [fitresult, gof] = createFit(a, b)
%CREATEFIT(A,B)
%  Create a fit.
%
%  Data for 'untitled fit 1' fit:
%      X Input : a
%      Y Output: b
%  Output:
%      fitresult : a fit object representing the fit.
%      gof : structure with goodness-of fit info.
%
%  See also FIT, CFIT, SFIT.

%  Auto-generated by MATLAB on 05-Mar-2021 17:02:05


%% Fit: 'untitled fit 1'.
[xData, yData] = prepareCurveData( a, b );

% Set up fittype and options.
ft = fittype( 'fourier1' );
opts = fitoptions( 'Method', 'NonlinearLeastSquares' );
opts.Display = 'Off';


% Fit model to data.
[fitresult, gof] = fit( xData, yData, ft, opts );




