%%%%%%%%%%%%%%%%%%%%%%%%%
% This scritp simulates multi-bvalue diffusion signals and adds noise
% this creates a magnitude signal for simulating IVIM fit results in imaging studies
% David Reiter 2017.03.06
% Revised for IVIM Whitepaper 2024.08.22
%%%%%%%%%%%%%%%%%%%%%%%%%

% bvalue range: 
bvalue_C = {[0 20 40 180 200 800]','Kidney';
            [0 10 120 200 1000]', 'Liver';
            [0 30 70 200 330 1000]', 'Muscle';
            [0 20 140 200 1000]', 'Breast malignant'; 
            [0 20 140 200 1000]', 'Breast benign'; 
            [0 40 60 200 290 1000]','Pancreas malignant'
            [0 40 60 200 290 1000]','Pancreas benign';
            [0 60 300 360 1000]','Brain'};

% mean parameters for kidney cortex, kidney medulla, liver, muscle, breast
% malignant, breast benign, pancrease malignant, pancreas benign, and
% brain (brain values averaged from Vieni paper)
% Currently, program only computes errors for mean tissue values
pf = [.189 .23 .103 .113 .07 .12 .2 .076];
ADC=[.00189 .00109 .00147 .00097 .00143 .0014 .00141 .00083]; 
APC=[.0405 .07 .0309 .0378 .0523 .0222 .0254 .0109]; 
SNR = logspace(log10(10), log10(1000), 15);
bcut=[190 190 190 190 190 190 190 190 300];

noise_std=1;
noise_N = 1000;

SEG_3parfit_result=zeros(length(pf),length(SNR),noise_N,3);

% loop structure steps through each organ first
for i=1:length(pf)
    bvalue = bvalue_C{i,1};
    for k=1:length(SNR)
        % MC loop
        for ll=1:noise_N
            e_r=noise_std*randn(length(bvalue),1);
            % assume signal and noise in real channel
            sig_r=SNR(k)*(pf(i)*exp(-bvalue.*APC(i))+(1-pf(i))*exp(-bvalue.*ADC(i)))+e_r;% forward model;

            % assume noise only in imaginary channel
            sig_i=noise_std*randn(length(bvalue),1);
            sig=sqrt(sig_r.^2+sig_i.^2);
            SEG_3parfit_result(i,k,ll,:) = IVIM_SEG_bifit_3par_opt1(bvalue,sig,bcut(i));

        end
        % display progress
        display([i k]);
    end
end
toc

%% create the solution matrix-true_par
true_par=zeros(length(pf),length(SNR),3);
for i=1:length(pf)
    for k=1:length(SNR)
       true_par(i,k,1) = pf(i);
       true_par(i,k,2) = APC(i);
       true_par(i,k,3) = ADC(i);
   end
end

save('whitepaper_sim_fin_fullb_20241125');



SEG_3parfit_ave_par = squeeze(mean(SEG_3parfit_result,3));
SEG_3parfit_acc_par = (SEG_3parfit_ave_par-true_par)./true_par;
flg=0;
SEG_3par_var_par = squeeze(std(SEG_3parfit_result,flg,3))./SEG_3parfit_ave_par;
save('whitepaper_acc_prec_sim_fin_fullb_20241125');


% Display Bias errors in f_p, D*, and D
figure; 
for i=1:length(pf)
    subplot(2,3,1); plot(SNR, SEG_3parfit_acc_par(i,:,1)*100);hold on;
end
xlabel('SNR'); ylabel('Rel Bias Error f_p (%)'); ylim([-10 10]);


for i=1:length(pf)
    subplot(2,3,2); plot(SNR, SEG_3parfit_acc_par(i,:,2)*100); hold on;
end
xlabel('SNR'); ylabel('Rel Bias Error D* (%)');ylim([-20 20]);

for i=1:length(pf)
    subplot(2,3,3); plot(SNR, SEG_3parfit_acc_par(i,:,3)*100); hold on;
end
xlabel('SNR'); ylabel('Rel Bias Error D (%)');ylim([-10 10]);


for i=1:length(pf)
    subplot(2,3,4); plot(SNR, SEG_3par_var_par(i,:,1)*100);hold on;
end
xlabel('SNR'); ylabel('Dispersion Error f_p (%)');ylim([-20 20]);

for i=1:length(pf)
    subplot(2,3,5); plot(SNR, SEG_3par_var_par(i,:,2)*100); hold on;
end
xlabel('SNR'); ylabel('Dispersion Error D* (%)');ylim([-20 20]);

for i=1:length(pf)
    subplot(2,3,6); plot(SNR, SEG_3par_var_par(i,:,3)*100); hold on;
end
xlabel('SNR'); ylabel('Dispersion Error D (%)');ylim([-20 20]);
legend(bvalue_C{:,2})


% find minimum SNR for error within 20% (i.e. [bias^2 + disp^2]^.5 )
SNRcut=.2;
min_SNR_fp_org=zeros(1,length(pf));
min_SNR_Dp_org=zeros(1,length(pf));
min_SNR_D_org=zeros(1,length(pf));
for i=1:length(pf) % for loop over each organ

    tmp_msnr_fp = sqrt(squeeze(SEG_3parfit_acc_par(i,:,1)).^2 + ...
                        squeeze(SEG_3par_var_par(i,:,1)).^2);
    minSNR_fp = SNR(find(tmp_msnr_fp<SNRcut));
    min_SNR_fp_org(i)=min(minSNR_fp);

    tmp_msnr_Dp = sqrt(squeeze(SEG_3parfit_acc_par(i,:,2)).^2 + ...
                        squeeze(SEG_3par_var_par(i,:,2)).^2);
    minSNR_Dp = SNR(find(tmp_msnr_Dp<SNRcut));
    min_SNR_Dp_org(i)=min(minSNR_Dp);

    tmp_msnr_D = sqrt(squeeze(SEG_3parfit_acc_par(i,:,3)).^2 + ...
                        squeeze(SEG_3par_var_par(i,:,3)).^2);
    minSNR_D = SNR(find(tmp_msnr_D<SNRcut));
    min_SNR_D_org(i)=min(minSNR_D);

end
[min_SNR_fp_org; min_SNR_Dp_org; min_SNR_D_org]


% find minimum SNR for Dispersion error within cutt-off % error
SNRcut=.2;
min_SNR_fp_org=zeros(1,length(pf));
min_SNR_Dp_org=zeros(1,length(pf));
min_SNR_D_org=zeros(1,length(pf));
for i=1:length(pf) % for loop over each organ

    tmp_msnr_fp = squeeze(SEG_3par_var_par(i,:,1));
    minSNR_fp = SNR(find(tmp_msnr_fp<SNRcut));
    min_SNR_fp_org(i)=min(minSNR_fp);

    tmp_msnr_Dp = squeeze(SEG_3par_var_par(i,:,2));
    minSNR_Dp = SNR(find(tmp_msnr_Dp<SNRcut));
    min_SNR_Dp_org(i)=min(minSNR_Dp);

    tmp_msnr_D = squeeze(SEG_3par_var_par(i,:,3));
    minSNR_D = SNR(find(tmp_msnr_D<SNRcut));
    min_SNR_D_org(i)=min(minSNR_D);

end
[min_SNR_fp_org; min_SNR_Dp_org; min_SNR_D_org]

% find minimum SNR for Bias error within cutt-off % error
SNRcut=.2;
min_SNR_fp_org=zeros(1,length(pf));
min_SNR_Dp_org=zeros(1,length(pf));
min_SNR_D_org=zeros(1,length(pf));
for i=1:length(pf) % for loop over each organ

    tmp_msnr_fp = squeeze(SEG_3parfit_acc_par(i,:,1));
    minSNR_fp = SNR(find(tmp_msnr_fp<SNRcut));
    min_SNR_fp_org(i)=min(minSNR_fp);

    tmp_msnr_Dp = squeeze(SEG_3parfit_acc_par(i,:,2));
    minSNR_Dp = SNR(find(tmp_msnr_Dp<SNRcut));
    min_SNR_Dp_org(i)=min(minSNR_Dp);

    tmp_msnr_D = squeeze(SEG_3parfit_acc_par(i,:,3));
    minSNR_D = SNR(find(tmp_msnr_D<SNRcut));
    min_SNR_D_org(i)=min(minSNR_D);

end
[min_SNR_fp_org; min_SNR_Dp_org; min_SNR_D_org]


