function [ U ] = split_step_fourier_method(start, step, ending, u0, h, k, n, GP_seq)
    u = u0;
    for m = start:step:ending % Start time loop
        if nargin > 7
            mu = GP_seq(floor(h*m)+1,:)';
            u = exp(h/2*1i*k*mu).*u; % Solve non-constant part of LSE
        end
        c = fftshift(fft(u)); % Take Fourier transform
        c = exp(-h/2/k*1i*(2*pi*n).^2).*c; % Advance in Fourier Space
        u = ifft(fftshift(c)); % Return to Physical Space
    end
    U =u;
end

