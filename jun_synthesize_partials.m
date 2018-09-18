function out = jun_synthesize_partials(Partials, time, L)
% ********************************************************************* %
% Synthesizes the partials, using linear amplitude interpolation and the
% phase interpolation method described in:
%
% R.McAulay and T. Quatieri, "Speech analysis/synthesis based on a
% sinusoidal representation," IEEE Trans. Acoustics, Speech, Signal
% Process., vol. 34, no. 4, pp. 744-754, 1986.
%
% Julian Neri, 180914
% McGill University, Montreal, Canada
% ********************************************************************* %

out = zeros(L,1);

M = size(Partials);
num_active = 0;

for m = 1:M    
    num_active_last = num_active;
    num_active = nnz(Partials{m}(:,1));
    
    if num_active > 0 && num_active_last == 0
        % Fade in all
        if m == 1
            time_m = [1 time(m)];
        else
            time_m = [time(m-1) time(m)];
        end
        H = abs(diff(time_m));
        % Fade out all
        n_unmatched = num_active;
        for idx = 1:n_unmatched
            % Fade in
            Omega_m = [Partials{m}(idx,1) Partials{m}(idx,1)];
            Amp_m = [-inf Partials{m}(idx,2)];
            Phase_m = [(Partials{m}(idx,3) - Omega_m(2)*H) Partials{m}(idx,3)];
                        
            ytemp = jun_synthesize_one_partial(Omega_m, Amp_m, Phase_m, H);
            out(time_m(1):time_m(2)-1,1) = out(time_m(1):time_m(2)-1,1) + ytemp(1:end-1);
        end
    elseif num_active == 0 && num_active_last > 0
        time_m = [time(m-1) time(m)];
        H = abs(diff(time_m));
        % Fade out all
        n_unmatched_last = num_active_last;
        for idx = 1:n_unmatched_last
            Omega_m = [Partials{m-1}(idx,1) Partials{m-1}(idx,1)];
            Amp_m = [Partials{m-1}(idx,2) -inf];
            Phase_m = [Partials{m-1}(idx,3) (Partials{m-1}(idx,3)+Omega_m(1)*H)];
            
            ytemp = jun_synthesize_one_partial(Omega_m, Amp_m, Phase_m, H);
            out(time_m(1):time_m(2),1) = out(time_m(1):time_m(2),1) + ytemp;
        end
    elseif num_active > 0 && num_active_last > 0
        [matches, unmatched_last, unmatched] = jun_match_sets(Partials{m-1}(1:num_active_last,4), Partials{m}(1:num_active,4));
        n_matches = size(matches,1);
        n_unmatched_last = length(unmatched_last);
        n_unmatched = length(unmatched);
        
        time_m = [time(m-1) time(m)];
        H = abs(diff(time_m));
        
        for idx = 1:n_unmatched
            % Fade in
            Omega_m = [Partials{m}(unmatched(idx),1) Partials{m}(unmatched(idx),1)];
            Amp_m = [-inf Partials{m}(unmatched(idx),2)];
            Phase_m = [(Partials{m}(unmatched(idx),3) - Omega_m(2)*H) Partials{m}(unmatched(idx),3)];
            
            ytemp = jun_synthesize_one_partial(Omega_m, Amp_m, Phase_m, H);
            out(time_m(1):time_m(2)-1,1) = out(time_m(1):time_m(2)-1,1) + ytemp(1:end-1);
        end
        
        for idx = 1:n_unmatched_last
            % Fade out
            Omega_m = [Partials{m-1}(unmatched_last(idx),1) Partials{m-1}(unmatched_last(idx),1)];
            Amp_m = [Partials{m-1}(unmatched_last(idx),2) -inf];
            Phase_m = [Partials{m-1}(unmatched_last(idx),3) Partials{m-1}(unmatched_last(idx),3)+Omega_m(1)*H];
            
            
            ytemp = jun_synthesize_one_partial(Omega_m, Amp_m, Phase_m, H);
            out(time_m(1):time_m(2),1) = out(time_m(1):time_m(2),1) + ytemp;             
        end
            
        for idx = 1:n_matches
            % Continue
            Omega_m = [Partials{m-1}(matches(idx,1),1) Partials{m}(matches(idx,2),1)];
            Amp_m = [Partials{m-1}(matches(idx,1),2) Partials{m}(matches(idx,2),2)];
            Phase_m = [Partials{m-1}(matches(idx,1),3) Partials{m}(matches(idx,2),3)];
                        
            ytemp = jun_synthesize_one_partial(Omega_m, Amp_m, Phase_m, H);
            out(time_m(1):time_m(2)-1,1) = out(time_m(1):time_m(2)-1,1) + ytemp(1:end-1);
        end
    end
end

end


function [ out ] = jun_synthesize_one_partial(Omega, Amplitude, Phase, H)
% Synthesizes one partial, given two points for frequency, amplitude and
% phase, H samples apart in time.

t_span = [0 H];
t = (t_span(1):t_span(2))';

% Linear Amplitude 
Amp_ = interp1(t_span, exp(Amplitude), t);

% (for eq. 34)
T_mx = [ 3/H^2  -1/H
        -2/H^3   1/H^2];

% (eq. 36)
M_star = (Phase(1) + Omega(1)*H - Phase(2))+(Omega(2) - Omega(1))*H/2;
M_star = round(1/(2*pi)*M_star);

% (for eq. 34)
rowv = [Phase(2) - Phase(1) - Omega(1)*H + 2*pi*M_star
        Omega(2) - Omega(1)];

% (eq. 34)
ab = T_mx*rowv;
alpha = ab(1);
beta = ab(2);

% (eq. 37) 
Phase_ = Phase(1) + Omega(1)*t + alpha*t.^2 + beta*t.^3;


% Synthesize the Partial
out = 2*real(exp(log(Amp_) + 1j*Phase_));

end