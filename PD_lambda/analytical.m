clc,clear all,close all
run("../Parameter_setting.m");


SNR_dB = 35:1:55;        % SNR = |et(rs)|^2/noise_var
SNR = 10.^(SNR_dB/10);

%%% analytical PD

PD_analytical = zeros(length(b),length(SNR));

figure
lambda_p = 2*(gain_2)*abs(ro)^2*SNR./(gain_2*sigma_ro2*SNR+1);
for i = 1:length(b)
    eta_p = 2*eta(i)/noise_var./(gain_2*sigma_ro2*SNR+1);
    PD_analytical(i,:) = marcumq(sqrt(lambda_p),sqrt(eta_p),3);
    plot(SNR_dB,PD_analytical(i,:))
    hold on
end


%%% expected PD

et = [1-2j;2+3j;3+5j];
et = et/vecnorm(et);    % et can be any normalized 3-D vector 
iter_num = 10000;
PD_exp = zeros(length(b),length(SNR));
for i = 1:iter_num
    R = ro*eye(3) + sqrt(sigma_ro2)/sqrt(2) * (randn(3) + 1j*randn(3));

    lambda = 2*gain_2*vecnorm(R*et)^2*SNR;
    for j = 1:length(b)
        PD_exp(j,:) = PD_exp(j,:)+marcumq(sqrt(lambda),b(j),3);
    end
end
PD_exp = PD_exp/iter_num;

for j = 1:length(b)
    scatter(SNR_dB,PD_exp(j,:));
    hold on
end

save('analytical.mat','-mat');
