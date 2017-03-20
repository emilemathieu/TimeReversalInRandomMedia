%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%% Propagation of time-harmonic waves in a random medium %%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%% PART I - Homogeneous medium %%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all; close all; clc;

xmax = 60;
k = 1; 
w = 1; 
L = 10; 
r0 = 2; 
N = 2^10 + 1;
x = linspace(-xmax/2, xmax/2, N);

%% Part 1 : Propagation of time-harmonic waves in a homogeneous medium

% theoretical part 
Rt = r0  * (1 + 4 * L^2/k^2/r0^4)^0.5;
psi_t = r0/Rt * exp(-2 * x.^2./Rt^2); 

% simulation 
psi_0 = exp(-x.^2/r0^2);
dx = xmax /(N - 1);
psi_0_fft = fftshift(fft(psi_0));

fmax = 1/2/dx;
f = 2 * pi .* linspace(-fmax, fmax, N ); 

sol_L_fft = psi_0_fft .* exp(-1i * f.^2/2/k * L); 
sol = ifft(sol_L_fft);
sol_norm = abs(sol).^2 ;

h = figure(1);
plot(x, psi_t, '-r', x, sol_norm, '-b', x, psi_0, '-k');
legend('Theoretical solution on z = L', 'Simulation on z = L', 'Wave on z = 0')
xlabel('x'); ylabel('Square of the amplitude');grid('on')
title('Wave on z = L');

%% Part 1 : Time reversal for time-harmonic waves in a homogeneous medium
rt = r0* (1 + 2i * L/k/r0^2)^0.5;
psi_L = r0/rt * exp(-x.^2/rt^2);
psi_L_conj = conj(psi_L);

% chi_M compactly supported
rM = 2:4:22;
sol_0 = zeros(length(rM), length(x));

for i = 1 : length(rM)
    chi_M = (1 - (x./2/rM(i)).^2).^2;
    chi_M(x <-2*rM(i)) = 0 ;
    chi_M(x > 2*rM(i)) = 0 ;

    term1_fft = fftshift(fft(psi_L_conj .* chi_M));
    sol_0_fft = term1_fft.* exp(1i * f.^2/2/k * L);
    sol_0(i, :) = abs(ifft(sol_0_fft));
end
h = figure(2);
plot(x, sol_0, 'LineWidth', 2);
legend('rM = 2',  'rM = 6',  'rM = 10',  'rM = 14', 'rM = 18', 'rM = 22')
title('Wave on z = 0, time reversal')
xlabel('x'); ylim([min(min(sol_0)), max(max(sol_0))]); ylabel('Amplitude'); grid('off')

% chi_M gaussian
rM = 2;
chi_M = exp(- (x.^2)/rM^2);
term1_fft = fftshift(fft(psi_L_conj .* chi_M));
sol_0_fft = term1_fft.* exp(1i * f.^2/2/k * L);
sol_0 = abs(ifft(sol_0_fft));
%atr = (1 + 4 * L^2/k^2/r0^2/rM^2 + 2i*L/k/rM^2)^0.5;
atr = (1 -4i * L/k/r0^2 -  4 * L^2/k^2/r0^2/rM^2 - 2i*L/k/rM^2)^0.5;
rtr_square = (1/rM^2 + 1/(r0^2 - 2i * L/k))^(-1) - 2i * L/k;
psi_tr_0 = abs(1/atr * exp(- x.^2/rtr_square));

h = figure(3);
plot(x, sol_0, '-r', x, psi_tr_0, '-b');
legend('Simulation', 'Theoretical solution')
xlabel('x'); ylabel('Amplitude'); grid('on')
title('Wave on z = 0, Gaussian time-reversal mirror')