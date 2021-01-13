function [Partials, time, padding, L, S] = jun_track_partials(y, fs, varargin)
% ********************************************************************* %
% Fast partial tracking of audio.
%
% Involves the following two main procedures for each short-term frame
%
% 1. Short-term sinusoidal model parameters are estimated using the 
%    Distribution Derivative Method (DDM).
% 2. Peaks are connected over consecutive analysis frames by solving an
%    assignment problem, using the Hungarian algorithm (munkres.m).
%
%
% INPUTS
% y: input audio signal
% fs: sampling rate (Hz)
% Variable inputs:
% N: length of short-term analysis window
% HopFactor: Hop Size = N/HopFactor
% OverSample: Oversampling factor
% G_g: Peaks below this magnitude threshold (in dB) are not considered.
% Q: short-term model polynomial order (DDM)
% delta: controls the relative preference towards useful assignments.
%        larger values of delta promote useful assignments (0<delta<1)
% zetaF, zetaA: the approximate cross-over values from a useful
%        to spurious assignment. 
%
% OUTPUTS
% Partials: a cell array, one cell for each frame containing the
% instantaneous estimates of frequency, amplitude, and phase, and index of
% the partials active in that frame.
% Time: the locations in time of the short-term frames
% Padding: the amount of zeros padding the signal
% L: the signal length
% S: the STFT (for plotting later on).
%
% These outputs are also written to a binary file using the 
% jun_write_partials() function.  Read the data with jun_read_partials()
%
%
% Reference:
% 
% J. Neri and P. Depalle, "Fast Partial Tracking of Audio with Real-Time
% Capability through Linear Programming", In Proceedings of the 21st
% International Conference on Digital Audio Effects (DAFx-18), Aveiro,
% Portugal, pp. 325-333, Sep. 2018.
%
% Julian Neri, 210101
% McGill University, Montreal, Canada
% ********************************************************************* %


% Sinusoidal model data will be written into the following binary file,
outputFileName = 'demo_partials.bin';
fileID = fopen(outputFileName,'w');

% and also Pseudo-SDIF text file. This can be converted to SDIF using
% the 'tosdif' executable.
sdifFileName = 'demo_partials.txt';
sdifFileID = fopen(sdifFileName,'w');
% SDIF Header
fprintf(sdifFileID,'SDIF\n\n\nSDFC\n\n');


%%%%%%%%%%%%%%%%%%    Parse Inputs    %%%%%%%%%%%%%%%%%%
numvarargs = length(varargin);
if numvarargs > 9
    error('Requires at most 9 optional inputs');
end

% Defaults
N = 2^11-1;
HopFactor = 4;
OverSample = 1;
G_g = -40;
Q = 2;
delta = .2;
zetaF = 50;
zetaA = 15;

% Use given parameters if available
optargs = {N HopFactor OverSample G_g Q delta zetaF zetaA};
optargs(1:numvarargs) = varargin;
[N, HopFactor, OverSample, G_g, Q, delta, zetaF, zetaA] = optargs{:};

% Check/Set Values so they don't produce errors later on.
N = max(4,round(N));
HopFactor = 2^nextpow2(max(1,HopFactor));
OverSample = 2^nextpow2(max(1,OverSample));
G_g = abs(G_g);
Q = max(1,round(Q));
% These can't be <= 0
zetaF(zetaF<1e-1) = 1e-1;
zetaA(zetaA<1e-1) = 1e-1;
% Convert to log from dB
zetaA = zetaA/20*log(10);

%%%%%%%%%%%%%%%%%%%%%%%%%%%% SETUP %%%%%%%%%%%%%%%%%%%%%%%

% Convert to mono
y = y(:,1);
L = length(y);

% Calculate variance (sigma^2) according to given zeta
varF = -zetaF^2*log((delta-1)/(delta-2));
varA = -zetaA^2*log((delta-1)/(delta-2));

% Analysis Windows
if OverSample > 1
    win = jun_bhwindow(N,1); % Window
    winD = jun_bhwindow(N,2); % Window Derivative
    HopFactor = max(8,HopFactor);
else
    win = jun_hannwindow(N,1); % Window
    winD = jun_hannwindow(N,2); % Window Derivative
    HopFactor = max(4,HopFactor);
end

% Calculate Hop Size from Hop Factor and N
H = round((N-mod(N,2))/HopFactor);
% Calculate the size of the DFT from Oversampling Factor and N
Ndft = 2^(round(log2(OverSample*2^round(log2(N)))));
if Ndft < N
    Ndft = 2^round(log2(N)+1); 
end
ndft = (0:Ndft-1)';

% Calculate the number of short-term frames and zero padding.
padding = zeros(2,1);
padding(1) = (Ndft-H);
M = ceil((Ndft+L - 2*H-1)/H + 1);
padding(2) = ((M-1)*H + Ndft) - (padding(1) + L);

% Zero-pad the signal
ypad = [zeros(padding(1),1); y; zeros(padding(2),1)];

% time, frequency, and STFT arrays
time = zeros(M,1);
omega = 2*pi*ndft/Ndft;
S = zeros(Ndft, M);

% Center of analysis window
n_center = (N - mod(N,2))/2+1;

% Zero-centered time vector of DFT
ndft = ndft - n_center;
% For centering the spectra (DDM)
centering = exp(1j*omega*n_center);
% DFT Matrix for alpha0 estimation (DDM)
ft_mat = exp(1j*ndft*omega(1:Ndft/2)');
% Number of Bins for DDM esimtate
Rddm = Q*OverSample;

% DDM Polynomial Array
% Zero-centered time vector (n in paper)
n = (0:N-1)' - n_center;
i = (0:Q);
% Polynomial time vector
p = bsxfun(@power, n, i);
% Derivative of p time vector
pD = bsxfun(@times, bsxfun(@power, n, (0:Q-1)), 1:Q);
% Time index of mid-point between analysis frames (H/2 in paper)
n_midpoint = [-floor(H/2) ceil(H/2)] + n_center;

% Approximate max in the STFT, without computing it directly
% (for thresholding the peak selection)
[~, max_sample] = max(abs(ypad));
maxS = 20*log10(max(abs(fft(win.*ypad(max_sample + (0:N-1)-n_center),Ndft))));
% Set relative amplitude threshold for peak selection
G_g = maxS - G_g;

% For storing the short-term estimates (alpha_ij in paper)
Alpha = 0;
num_peaks = 0;

% For storing the parameters/index of the partials for each frame
Partials = cell(M,1);
pairs_ = 0;
tracks_ = 0;
num_tracks = 0;
num_active = 0;

% Write header information to file
fwrite(fileID,num_tracks,'uint');
fwrite(fileID,L,'uint');
fwrite(fileID,fs,'double');
fwrite(fileID,M,'uint');
fwrite(fileID,padding(1),'uint');
fwrite(fileID,padding(2),'uint');


%%%%%%%%%%%%%%%%%%%%% PARTIAL TRACKING STARTS HERE %%%%%%%%%%%%%%%%%%%%%%%

% For every short-term analysis frame 'm'
for m = 1:M
    
    %%%%%%%%%%%%%%%%%% SHORT-TERM PARAMETER ESTIMATION %%%%%%%%%%%%%%%%%%
    % Time at the center of the short-term analysis frame
    time(m) = 1 + (m-1)*H + n_center;
    % Short-term signal
    yst = ypad(time(m)-n_center:time(m)+N-1-n_center);

    % Save estimates from the previous frame
    Alpha_last = Alpha;
    num_peaks_last = num_peaks;
    
    % Estimate short-term model parameters of each peak using DDM
    [Alpha, num_peaks, S(:,m)] = jun_ddm(yst, Q, Rddm, G_g, win, winD, ...
                                    p, pD, centering, Ndft, ft_mat,omega);
                                
    %%%%%%%%%%%%%%%%%% TRACKING %%%%%%%%%%%%%%%%%%
    Partials{m}(:,1) = 0;
    if m > 1
        %%%%%%%%%%%%%%%%%% PEAK-to-PEAK ASSIGNMENTS %%%%%%%%%%%%%%%%%%
        num_assignments = 0;
        if num_peaks > 0 && num_peaks_last > 0
            
            % Amplitudes/Frequencies at midpoint of analysis frames
            % a(+H/2) [k-1]
            mA1 = real(p(n_midpoint(2),:)*Alpha_last);
            % a(-H/2) [k]
            mA2 = real(p(n_midpoint(1),:)*Alpha);
            % f(+H/2) [k]
            mF1 = fs/(2*pi)*imag(pD(n_midpoint(2),:)*Alpha_last(2:end,:));
            % f(-H/2) [k]
            mF2 = fs/(2*pi)*imag(pD(n_midpoint(1),:)*Alpha(2:end,:));

            % DELTA_a and DELTA_f
            % eq. (10)
            DeltaA = bsxfun(@minus, mA1', mA2);
            % eq. (11)
            DeltaF = bsxfun(@minus, mF1', mF2);
            
            % USEFUL COST
            % eq. (8)
            A_useful = 1-exp(-DeltaF.^2/varF - DeltaA.^2/varA);
            % SPURIOUS COST
            % eq. (9)
            B_spurious = 1-(1-delta)*A_useful;

            % COST MATRIX
            % eq. (14)
            [Costs, Type] = min([B_spurious(:) A_useful(:)],[],2);
            
            C_mtx = reshape(abs(Costs), num_peaks_last, num_peaks);

            % HUNGARIAN ALGORITHM SOLVES THE ASSIGNMENT PROBLEM
            jj = munkres(C_mtx);

            ii = find(jj);
            ind = ii + num_peaks_last*(jj(ii)-1);
            ii = ii(Type(ind)==2);

            % Stores the useful assignments
            num_assignments = length(ii);
            Assignments = [ii' jj(ii)'];
        end
            
        %%%%%%%%%%%%%%%%%% TRAJECTORY LABELING %%%%%%%%%%%%%%%%%%
        num_active_last = num_active;
        pairs_last = pairs_;
        tracks_last = tracks_;
        
        num_active = 0;
        pairs_ = zeros(num_assignments,1);
        tracks_ = zeros(num_assignments,1);
        
        if num_assignments > 0
            % An assignment either continues an existing partial or starts 
            % a new one (birth).
            [toContinue, ~, toBirth] = jun_match_sets(pairs_last, Assignments(:,1));
            
            num_toContinue = length(toContinue(:,1));
            num_toBirth = length(toBirth);
            num_active = num_toContinue + num_toBirth;
                        
            pairs_ = Assignments([toContinue(:,2); toBirth],2);
            tracks_ = [tracks_last(toContinue(:,1)); num_tracks + (1:num_toBirth)'];
            
            ii_active_last = num_active_last + (1:num_toBirth);
            
            % Save partials for frame m-1 and m
            Partials{m-1}(ii_active_last,1:3) = jun_framph(Alpha_last(:,Assignments(toBirth,1)), p(n_center+1,:), pD(n_center+1,:));
            Partials{m-1}(ii_active_last,4) = tracks_(num_toContinue+1:num_active);
            
            Partials{m}(1:num_active,1:3) = jun_framph(Alpha(:,pairs_), p(n_center+1,:), pD(n_center+1,:));
            Partials{m}(1:num_active,4) = tracks_;
            
            num_active_last = num_active_last + num_toBirth;
            num_tracks = num_tracks + num_toBirth;
        end
        % Write data into file for frame m-1 
        jun_write_partials(fileID, m-1, time(m-1), num_active_last, Partials{m-1});
        time_sdif = time(m-1) - (padding(1)+1);
        if time_sdif >= 0 && time_sdif <= (L + H)
            jun_write_sdif(sdifFileID, fs,time_sdif/fs, num_active_last, Partials{m-1});
        end
    end
end
% Write data into file for frame M
jun_write_partials(fileID, M, time(M), num_active, Partials{M});
% Write total number of partials detected
fseek(fileID,0,'bof');
fwrite(fileID,num_tracks,'uint');
fclose(fileID);

% SDIF Footer
fprintf(sdifFileID,'ENDC\nENDF');
% Close SDIF text File
fclose(sdifFileID);
end
% eof