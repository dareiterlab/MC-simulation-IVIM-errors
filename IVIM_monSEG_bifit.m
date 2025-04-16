function [ivim_simple_calculated] = IVIM_monSEG_bifit(b,y,bcut)

% normalize signal
y_e = y/max(y);

% define b-value cut point for segmented fit
i_blg = find(b>bcut);
b_lg = b(i_blg);

% fit constraints for stage 1
start_p = [0.8 2e-3];
lowr = [0.1 0.1e-3]; 
upr = [3 5e-3];

% STAGE 1 FIT: large bvalues for estimate of D
g_stg1 = fittype('a*exp(-x*b)');
opt_s1=fitoptions('method','NonlinearLeastSquares','Lower',lowr,...
    'Upper',upr,'Startpoint', start_p);

% stage 1 fit
Stg1FitModel = fit(b_lg, y_e(i_blg), g_stg1, opt_s1);

% STAGE 2 FIT: calculate fp based on monoexponential fit
% use small bvalues 
fp = y_e(1)-Stg1FitModel.a;

ivim_simple_calculated = [fp Stg1FitModel.b]; 




