function [ u ] = split_step_fourier_method(start, step, ending, u, h, k, GP_seq)

    for m = start:step:ending % Start time loop
        if nargin > 6
            mu = GP_seq(floor(h*m)+1,:)';
            u = exp(h/2*1i*k*mu).*u; % Solve non-constant part of LSE
        end
        c = fftshift(fft(u)); % Take Fourier transform
        c = exp(-h/2*1i*k).*c; % Advance in Fourier Space
        u = ifft(fftshift(c)); % Return to Physical Space
    end

end

