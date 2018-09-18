function [Partials, time, padding, L, fs, num_partials] = jun_read_partials(file_name)
% ********************************************************************* %
% Reads from file, accumulates frame data into Partials cell array.
%
% Julian Neri, 180914
% McGill University, Montreal, Canada
% ********************************************************************* %

fileID = fopen(file_name,'r');

% Read header information
num_partials = fread(fileID,1,'uint');
L = fread(fileID,1,'uint');
fs = fread(fileID,1,'double');
M = fread(fileID,1,'uint');
padding(1) = fread(fileID,1,'uint');
padding(2) = fread(fileID,1,'uint');

Partials = cell(M,1);
time = zeros(M,1);

frame = fread(fileID,1,'uint');
while frame
    time(frame) = fread(fileID,1,'uint');
    rows = fread(fileID,1,'uint');
    
    Partials{frame} = zeros(rows,4);
    for r = 1:rows
        Partials{frame}(r,1) = fread(fileID,1,'double');
        Partials{frame}(r,2) = fread(fileID,1,'double');
        Partials{frame}(r,3) = fread(fileID,1,'double');
        Partials{frame}(r,4) = fread(fileID,1,'uint');
    end
    frame = fread(fileID,1,'uint');
end
fclose(fileID);
end