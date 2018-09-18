function jun_plot_partials(Partials, time, fs, num_tracks)
% Plots the partials
% 
%
% Julian Neri, 180914
% McGill University, Montreal, Canada

M = length(time);

track_F = -inf+zeros(M,num_tracks);
track_A = -inf+zeros(M,num_tracks);


for m = 1:M    
    active = Partials{m}(:,4);

    track_F(m,active) = Partials{m}(:,1);
    track_A(m,active) = Partials{m}(:,2);
end

plot(time/fs, fs/(2*pi)*track_F, 'linewidth',1.5);
axis tight; grid on; box on;


end