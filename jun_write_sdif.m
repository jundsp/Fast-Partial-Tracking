function jun_write_sdif(fileID,fs,time,num_tracks,Partials)
% ********************************************************************* %
% Writes frame block to sdif text file.
%
%
% Julian Neri, 210101
% McGill University, Montreal, Canada
% ********************************************************************* %

fprintf(fileID,'1TRC\t1\t0\t%f\n',time);
fprintf(fileID,'  1TRC\t0x0004\t%d\t4\n',num_tracks);
for i = 1:num_tracks
    f = Partials(i,1)*fs/(2*pi);
    a = exp(Partials(i,2));
    p = Partials(i,3);
    track_id = Partials(i,4);
    fprintf(fileID,'\t%d\t%f\t%f\t%f\n',track_id,f,a,p);
end
fprintf(fileID,'\n');

end

