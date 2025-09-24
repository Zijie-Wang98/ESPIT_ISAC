%%% physical constant
c = 3e8;
fc = 12e9;
lambda = c/fc;
k = 2*pi/lambda;
noise_var = 5.6e-3;

Z0 = 376.7;

%%% antenna constant
Lx = 30*15e-3;
Ly = 30*15e-3;
N = 30*30;

%%% system setup

Pt = 100e-3;
ro = 0.2;
sigma_ro2 = 0.01;

u = [5;-5;15];
r = [10;5;30];

sx = -Lx/2+(0.5+0:sqrt(N-1))*Lx/sqrt(N);
sy = -Ly/2+(0.5+0:sqrt(N-1))*Ly/sqrt(N);

% g = zeros(N,1);
% for i = 1:sqrt(N)
%     for j = 1:sqrt(N)
%         p = vecnorm([sx(i);sy(j);0]-r);
%         g((i-1)*sqrt(N)+j) = exp(1j*k*p)/(4*pi*p);
%     end
% end
% 
% gain_2 = sum(abs(g).^2);
gain_2 = 0.0056;

b = [5.2779  5.5224, 5.7539, 5.9744]; % Pfa = -4 -4.5 -5 -5.5   sqrt{2\eta/noise_var^2}
eta = b.^2*noise_var/2;
