%% Propagation of time-harmonic waves in a random medium.
clc; clear all; close all; clf;

%%% Constants
h = 1; % Longitudinal (along z) step size
zc = 1;
xc = 4;
sigma = 1;
L = 10; % Position of mirror
k = 1; % Frequency ?
omega = 1; % Frequency
r0 = 2; % Gaussian beam radius
N = 2^10; % Number of points for discretization
xmax = 60;

n = (-N/2:N/2-1)';
x = n * 2 * xmax / N;

u0 = exp(-x.^2/r0^2); % Initial profile

rM = 2;
mirror = exp(-x.^2/rM^2); % Gaussian mirror profile

%%% split-step Fourier method
nb_GP = round(L / zc) + 1; % Number of Gaussian processes
nb_MC = 100; % Number of profiles to be averaged
U_L = zeros(nb_MC, N); % Transmitted wave profile
U_0_rand = zeros(nb_MC, N); % Refocused wave profile random medium
U_0_homo = zeros(nb_MC, N); % Refocused wave profile homogeneous medium
for i = 1:nb_MC
    GP_seq = sample_GP(x, sigma, xc, nb_GP);
    
    % Go forward in random medium
    U_L(i,:) = split_step_fourier_method(0, 1, round(L/h), u0, h, k, GP_seq)';
    
    % Go reverse in same random medium
    u = conj(U_L(i,:)') .* mirror;
    U_0_rand(i,:) = split_step_fourier_method(round(L/h), -1, 0, u, h, k, GP_seq)';
    
    % Go reverse in homogeneous medium
    u = conj(U_L(i,:)') .* mirror;
    U_0_homo(i,:) = split_step_fourier_method(round(L/h), -1, 0, u, h, k)';
end
U_L = mean(U_L, 1)';
U_0_rand = mean(U_0_rand, 1)';
U_0_homo = mean(U_0_homo, 1)';

%% Mean transmitted wave profile
rt = r0*sqrt(1 + 2*1i*L/k/r0^2);
gamma0 = sigma^2*zc;
mean_wave_L = r0/rt * exp(-x.^2/rt^2) * exp(-gamma0*omega^2*L/8);
figure(1); plot(x, abs(U_L), x,abs(mean_wave_L))
legend('empirical','theoretical')
title('Mean transmitted wave profile in random medium')

%% mean refocused wave profile in z=0 (random medium)
atr = sqrt(1 + 4*L^2/(k*r0*rM)^2 + 2*1i*L/k/rM^2);
gamma2 = 2*sigma^2*zc/xc^2;
ra_square = 48/L/gamma2/omega^2;
rtr_square =(1/rM^2+1/(r0^2-2*1i*L/k))^-1 + 2*1i*L/k;
mean_wave_0_rand = 1/atr * exp(-x.^2/rtr_square) .* exp(-x.^2/ra_square);

figure(2); plot(x, abs(U_0_rand), x, abs(mean_wave_0_rand));
legend('empirical','theoretical')
title('Mean refocused wave profile in random medium')

%% mean refocused wave profile in z=0 (homogeneous medium)
mean_wave_0_homo = 1/atr * exp(-x.^2/rtr_square) .* exp(-gamma0*omega^2*L/8);
figure(3); plot(x, abs(U_0_homo), x, abs(mean_wave_0_homo));
legend('empirical','theoretical')
title('Mean refocused wave profile in homegeneous medium')

%% Time reversal for time-dependent waves in a random medium.
omega0 = 1;
B = 0.75;
nb_discr = 20;
omega_discr = linspace(omega0-B,omega0+B,nb_discr);

U_L_w = zeros(nb_discr, N); % Transmitted wave profile
U_0_w = zeros(nb_discr, N); % Refocused wave profile random medium

for w = 1:nb_discr
    omega = omega_discr(w);
    k = omega_discr(w);

    GP_seq = sample_GP(x, sigma, xc, nb_GP);

    % Go forward in random medium
    U_L_w(w,:) = split_step_fourier_method(0, 1, round(L/h), u0, h, k, GP_seq)';

    % Go reverse in same random medium
    u = conj(U_L_w(w,:)') .* mirror;
    U_0_w(w,:) = split_step_fourier_method(round(L/h), -1, 0, u, h, k, GP_seq)';
end

U_0_w = sum(U_0_w,1);
figure(4); plot(x, abs(U_0_w));
%legend('empirical','theoretical')
title('Refocused time-dependent waves in random medium')