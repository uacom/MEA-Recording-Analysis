function [Results] = Cal_FreGC_Cui(data,ParaCau)
% Calculate the granger causality in the frequency domain based on the
% toolbox of MVGC_V1.0


X = data;
[nvars,nobs] = size(data);

% Model order estimation
[AIC,BIC] = tsdata_to_infocrit(X,ParaCau.momax,ParaCau.icregmode);
[~,bmo_AIC] = min(AIC);
[~,bmo_BIC] = min(BIC);
if strcmpi(ParaCau.morder,'AIC')
    morder = bmo_AIC;
    fprintf('\nusing AIC best model order = %d\n',morder);
elseif strcmpi(ParaCau.morder,'BIC')
    morder = bmo_BIC;
    fprintf('\nusing BIC best model order = %d\n',morder);
end

% Granger causality calculation
[A,SIG] = tsdata_to_var(X,morder,ParaCau.regmode); % Estimate VAR model of selected order from data.
[G,~] = var_to_autocov(A,SIG); % Autocovariance calculation

F = autocov_to_pwcgc(G); % Granger causality calculation: time domain
pval = mvgc_pval(F,morder,nobs,1,1,1,nvars-2); % Significance test using theoretical null distribution, adjusting for multiple hypotheses.
sig  = significance(pval,ParaCau.alpha,ParaCau.mhtc);
% figure(1); clf;
% subplot(1,3,1);
% plot_pw(F);
% title('Pairwise-conditional GC');
% subplot(1,3,2);
% plot_pw(pval);
% title('p-values');
% subplot(1,3,3);
% plot_pw(sig);
% title(['Significant at p = ' num2str(ParaCau.alpha)])

f = autocov_to_spwcgc(G,ParaCau.fres); % Granger causality calculation: frequency domain
f_f = sfreqs(ParaCau.fres,ParaCau.Fs)';
% figure(2); clf;
% plot_spw(f,ParaCau.Fs,[0 80]);
idex = f_f>= ParaCau.fpass(1) & f_f <= ParaCau.fpass(2); 
f_new = f(:,:,idex);
f_f_new = f_f(idex);

Results.TimeGC.F = F;
Results.TimeGC.SigF = sig;
Results.FreGC = f_new;
Results.FreGC_Fres = f_f_new;

