function [bifit_calculated] = IVIM_SEG_bifit_3par_opt1(b,y,bcut)

% normalize signal
y_e = y/max(y);

% added a new variable "bcut" that can be passed in to define segmented fit
% cut poing
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
% use mono-exponential fit to estimate perfusion fraction
fp = y_e(1)-Stg1FitModel.a;

% STAGE 2 FIT: small bvalues for fp and D*
% use small bvalues 
[b_sm, i_bsm] = setdiff(b,b_lg);

% fit constraints for stage 2
start_p = [2e-2];
lowr = [0.5e-2];
upr = [2];

g_stg2 = fittype('a*exp(-x*b) + c*exp(-x*d)', 'problem', {'a','c','d'});
opt_s2=fitoptions('method','NonlinearLeastSquares','Lower',lowr,...
    'Upper',upr,'Startpoint', start_p); %, 'problem', Stg1FitModel.b);

% stage 2 fit
Stg2FitModel = fit(b_sm,y_e(i_bsm), g_stg2, opt_s2, 'problem', {fp, Stg1FitModel.a, Stg1FitModel.b});

bifit_calculated = [fp Stg2FitModel.b Stg1FitModel.b]; 




