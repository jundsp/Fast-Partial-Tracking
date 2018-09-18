function out = jun_framph(alpha, p, p_prime)
% Calculates the frequency, amplitude, and phase from alpha.
%
% p = n^i for i=0:Q
% pD = dp/dn
%
% Julian Neri, 180914
% McGill University, Montreal, Canada

% eq. (7)
frequency = imag(p_prime*alpha(2:end,:));
% eq. (5)
amplitude = real(p*alpha);
% eq. (6)
phase =     imag(p*alpha);

out = [frequency' amplitude' phase'];
    
end