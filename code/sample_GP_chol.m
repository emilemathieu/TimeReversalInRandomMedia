function [ seq ] = sample_GP_chol( S, nb )
% Use fft to sample Gaussian processes

    N = size(S, 1);
    seq = zeros(nb, N);
    for i=1:nb
        W=randn(1,N);
        seq(i,:) = S * W';
    end 

end

