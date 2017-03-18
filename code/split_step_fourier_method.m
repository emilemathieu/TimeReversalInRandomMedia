function [ U ] = split_step_fourier_method(start, step, ending, u0, h, k, grid_f, GP_seq)
    u = u0;
    grid_f = fftshift(grid_f);
    c = fft(u);
    for m = start:step:ending % Start time loop
        c = exp(-h/2/k*1i*(grid_f).^2).*c; % Advance in Fourier Space
        u=ifft(c);
        if nargin > 7
            mu = GP_seq(floor(h*m)+1,:)';
            u = exp(h/2*1i*k.*mu).*u; % Solve non-constant part of LSE
        end
        c=fft(u);
    end
    U = ifft(c);

