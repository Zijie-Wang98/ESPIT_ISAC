%%% Note: The scripts for figure plots are not included.
%%% Note: The scripts for figure plots are not included.
%%% Note: The scripts for figure plots are not included.

run("../../Parameter_setting.m");
load("../../Kmatrix_values.mat");

SNR_dB = 0:0.5:50;
SNR = 10.^(SNR_dB/10);

% sensinhg-opt

[eigenvec_r,eigenval_r] = eig(Krr);
[~,idx] = max(diag(eigenval_r));
k_r = eigenvec_r(:,idx);
mu_r = eigenval_r(idx,idx);       % \kappa_max{Kss}


lambda_s = 2*gain_2*abs(ro)^2*SNR*mu_r./(gain_2*sigma_ro2*SNR*mu_r+1); 
eta_original = 2/noise_var./(gain_2*sigma_ro2*SNR*mu_r+1);                     %% eta'/eta

PDs = marcumq(sqrt(real(lambda_s)),sqrt(real(eta_original*eta(3))),3);
Rs = real(log2(1+k_r'*Kur'*Kur*k_r/mu_r*SNR));

% communication-opt

[eigenvec_u,eigenval_u] = eig(Kuu);
[~,idx] = max(diag(eigenval_u));
k_u = eigenvec_u(:,idx);
mu_u = eigenval_u(idx,idx);     % \kappa_max{Kcc}


% noise_var normalized E-field power
mags = k_u'*Kur*Kur'*k_u/mu_u;
mags = SNR*abs(mags);
lambda_c = 2*gain_2*abs(ro)^2*mags./(gain_2*sigma_ro2*mags+1);
eta_c = 2/noise_var./(gain_2*sigma_ro2*mags+1);

PDc = marcumq(sqrt(real(lambda_c)),sqrt(real(eta_c*eta(3))),3);
Rc = log2(1+mu_u*SNR);

save('SNR_data.mat','-mat');
