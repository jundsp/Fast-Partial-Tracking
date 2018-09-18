function [win, n] = jun_bhwindow(N,d)
% Blackman-Harris Window
%
% N: length of window
% d: order of derivative 
%
% Julian Neri, 180914
% McGill University, Montreal, Canada

if nargin < 2
    d = 1;
end

if mod(N,2)
    Nwinhf = (N-1)/2;
    n = (-Nwinhf:Nwinhf)';
else
    Nwinhf = N/2-1;
    n = (-Nwinhf-1:Nwinhf)';
end

a = [0.35875 0.48829 0.14128 0.01168];
in = [2 4 6]*pi/(N);

switch(d)
    case 1
        win = a(1) + a(2)*cos(in(1)*n)+a(3)*cos(in(2)*n) + a(4)*cos(in(3)*n);
    case 2
        win = -a(2)*in(1)*sin(in(1)*n) - a(3)*in(2)*sin(in(2)*n) - a(4)*in(3)*sin(in(3)*n);
    case 3
        win = -a(2)*in(1)^2*cos(in(1)*n) - a(3)*in(2)^2*cos(in(2)*n) - a(4)*in(3)^2*cos(in(3)*n);
end

end