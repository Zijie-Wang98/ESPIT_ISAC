clc,clear all
run("../../Parameter_setting.m");
load("../../Kmatrix_values.mat");
delta_list = 0:0.0001:1;
%delta_list = [linspace(0.1,0.1276,10) linspace(0.1276,0.1279,100) linspace(0.1279,5,5)];
mags_eu = zeros(1,length(delta_list));
mags_er = mags_eu;

cd ../..
for d = 1:length(delta_list)
    [eu,er,~] = K_delta(delta_list(d),0,0);
    mags_eu(d) = vecnorm(eu)^2;
    mags_er(d) = vecnorm(er)^2;
end
cd Pareto_bound/Main
R = log2(1+mags_eu/noise_var);

Pfa_list = [1e-4, 10^(-4.5), 1e-5, 10^(-5.5)];
lambda = 2*gain_2*abs(ro)^2*mags_er./(gain_2*sigma_ro2*mags_er+noise_var); 
eta_original = 2./(gain_2*sigma_ro2*mags_er+noise_var);                     %% eta'/eta
for l = 1:length(Pfa_list)
    PD(l,:) = marcumq(sqrt(lambda),sqrt(eta_original*eta(l)),3);
end


save('Pareto.mat','-mat');
