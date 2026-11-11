function rsnr = jun_rsnr(x, y)
% Reconstruction signal-to-noise ratio (dB) of y as an estimate of x.
%
% x: reference signal
% y: reconstructed signal (same size as x)
%
% rsnr = 10*log10( sum(x.^2) / sum((x-y).^2) )
%
% Julian Neri
% McGill University, Montreal, Canada

rsnr = 10*log10(sum(x(:).^2) / sum((x(:)-y(:)).^2));

end
