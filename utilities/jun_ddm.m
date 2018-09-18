function [Alpha, num_peaks, S] = jun_ddm(y, Q, R, G_g, win, winD, ...
                                    p, pD, centering, Ndft, ft_mat,omega)
% Distribution Derivative Method of estimating sinusoidal model parameters
%
% Reference:
%
% M. Betser, "Sinusoidal Polynomial Parameter Estimation Using the 
% Distribution Derivative", IEEE Transactions on Signal Processing, vol.
% 57, no. 12, pp. 4633-4645, Dec. 2009.
%
%
% y: input signal
% Q: polynomial order
% R: number of peaks to use for estimation
% G_g: peak amplitude threshold (dB)
% win: window
% winD: derivative of window
% p: time vector n^i, where i = 0:Q
% pD: time derivative of p
%
% Julian Neri, 180914
% McGill University, Montreal, Canada

% Input info
N = length(y);
max_y = max(abs(y))*2;

% Windowed signal
yst = win.*y;
% DFT (one frame of the signal's STFT)
S = fft(yst, Ndft);

% Pick Peaks
[peakBin, num_peaks, LeftBin, RightBin] = jun_pick_peaks(S(1:Ndft/2), G_g);

% Estimate parameters for each peak (alpha_ij in paper), DDM
Alpha = 0;
if num_peaks > 0
    % Stores the various DFTs involved with the DDM method
    Sp = zeros(Ndft, Q-1);
    % Polynomials * signal
    Sp(:,1) = S;
    Sp(:,1) = Sp(:,1).*centering;
    for i = 2:Q
        Sp(:,i) = fft(yst.*pD(:,i), Ndft);
        Sp(:,i) = Sp(:,i).*centering;
    end

    % Derivative Windowed Signal 
    yDst = winD.*y;
    SD = fft(yDst, Ndft);
    SD = SD.*centering; 
    SD = SD + Sp(:,1).*(-1j*omega);

    alpha_hat = zeros(Q+1,num_peaks);
    useful = 0;
    for jj = 1:num_peaks

        % Use the highest energy bins around the peak for DDM
        pbl = max(LeftBin(jj), peakBin(jj)-(R-1));
        pbr = min([RightBin(jj) peakBin(jj)+(R-1) Ndft/2]);
        pb_sides = sort(pbl:pbr, 'Descend');
        pbs = [peakBin(jj) pb_sides];

        % Define matrix A and vector b
        A = Sp(pbs, 1:Q);
        b = -SD(pbs);

        % Solve for alpha using least squares
        alpha_temp = [0; A\b];


        % Check edges
        if any(isnan(alpha_temp))
            fprintf('NAN estimate, skipping peak %0.0f of frame %0.0f \n', jj, m);
        elseif any(abs(p([1 N],2:end)*alpha_temp(2:end))>1e100)
            fprintf('Infinite edges, skipping peak %0.0f of frame %0.0f \n', jj, m)
        else
            % Solve for alpha0 (approximate)
            gam = win.*exp(p(:,2:end)*alpha_temp(2:end));
            T_gam = ft_mat(1:N,peakBin(jj))'*gam;

            alpha0 = log(Sp(peakBin(jj),1)) - log(T_gam);
            alpha_temp(1) = alpha0;

            % Final check and save
            if all(p([1 round(N/2) N],:)*real(alpha_temp) <= max_y)
                useful = useful + 1;
                alpha_hat(:,useful) = alpha_temp;
            end
        end 
    end

    % Save the estimates
    num_peaks = useful;
    Alpha = alpha_hat(:,1:num_peaks);
end
end