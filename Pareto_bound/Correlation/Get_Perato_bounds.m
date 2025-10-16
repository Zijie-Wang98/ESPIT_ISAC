%%% Note: The scripts for figure plots are not included.
%%% Note: The scripts for figure plots are not included.
%%% Note: The scripts for figure plots are not included.

clc,clear all
run("../../Parameter_setting.m");

delta_list = 0:0.0001:1;
N_num = 3;
rho_list = zeros(N_num,1);
mags_eu = zeros(N_num,length(delta_list));
mags_er = mags_eu;
R = mags_eu;
PD = mags_eu;

for n = 1:N_num
    load(['Kmatrix_values',num2str(n),'.mat']);
    rho_list(n) = pho;
    for d = 1:length(delta_list)
        [eu,er,~] = K_delta(Krr,Kuu,Kur,delta_list(d),0);
        mags_eu(n,d) = vecnorm(eu)^2;
        mags_er(n,d) = vecnorm(er)^2;
    end
    R(n,:) = log2(1+mags_eu(n,:)/noise_var);
    lambda = 2*gain_2*abs(ro)^2*mags_er(n,:)./(gain_2*sigma_ro2*mags_er(n,:)+noise_var); 
    eta_original = 2./(gain_2*sigma_ro2*mags_er(n,:)+noise_var);                     %% eta'/eta
    PD(n,:) = marcumq(sqrt(lambda),sqrt(eta_original*eta(3)),3);
end

save('Pareto.mat','-mat');
