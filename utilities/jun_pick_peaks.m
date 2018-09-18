function [peak_inds, n_peaks, left_inds, right_inds]  = jun_pick_peaks(X, G_g)
% Simple peak-picking routine
%
% A peak is a local maximum.
%
% X: spectrum
% G_g: a peak whose magnitude is less than G_g is disregarded.
%
% Julian Neri, 180914
% McGill University, Montreal, Canada

N = size(X,1);
M = size(X,2);

idx_peaks = zeros(N, M);
n_peaks = zeros(1, M);

LeftBin = zeros(N,M);
RightBin = zeros(N,M);

X = 20*log10(abs(X));

for m = 1:M
    num = 0;
    num_temp = 0;
    for n = 2:N-1
        if X(n, m) > G_g
            if X(n, m) > X(n-1, m) && X(n, m) > X(n+1, m)
                
                num_temp = num_temp + 1;

                left = n-1;
                while (left-1 > 0) && (X( left, m)>X( left-1, m))
                    left = left - 1;
                end

                right = n+1;
                while (right+1 < N+1) && (X( right, m) > X( right+1, m))
                    right = right + 1;
                end
                
                num = num + 1;
                idx_peaks(num, m) = n;
                LeftBin(num,m) = left;
                RightBin(num,m) = right;
            end
        end
    end
    n_peaks(m) = num;
end

maxpeaks = max(n_peaks);
peak_inds = idx_peaks(1:maxpeaks,:);
left_inds = LeftBin(1:maxpeaks,:);
right_inds = RightBin(1:maxpeaks,:);

end
