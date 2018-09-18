function jun_write_partials(fileID,frame,time,rows,Partials)
% ********************************************************************* %
% Writes frame block into file.
%
%
% Julian Neri, 180914
% McGill University, Montreal, Canada
% ********************************************************************* %

fwrite(fileID,frame,'uint');
fwrite(fileID,time,'uint');
fwrite(fileID,rows,'uint');
for r = 1:rows
    fwrite(fileID,Partials(r,1),'double');
    fwrite(fileID,Partials(r,2),'double');
    fwrite(fileID,Partials(r,3),'double');
    fwrite(fileID,Partials(r,4),'uint');
end

end