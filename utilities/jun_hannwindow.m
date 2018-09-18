function win = jun_hannwindow(N,derivative)
% Hann Window
%
% N: length of window
% d: order of derivative 
%
% Julian Neri, 180914
% McGill University, Montreal, Canada

if nargin < 2
    derivative = 1;
end

if mod(N,2) == 1
    Nhf = (N-1)/2;

    n = (0:N-1)' - Nhf;
    in = 2*pi/(N-1);

    switch(derivative)
        case 1
            win = .5 + .5*cos(in*n);
        case 2
            win = -.5*in*sin(in*n);
        case 3
            win = -.5*in^2*cos(in*n);
    end

else
    n = (0:N-1)';
    in = 2*pi/(N-1);

    switch(derivative)
        case 1
            win = .5 - .5*cos(in*n);
        case 2
            win = .5*in*sin(in*n);
        case 3
            win = .5*in^2*cos(in*n);
    end
end
end

