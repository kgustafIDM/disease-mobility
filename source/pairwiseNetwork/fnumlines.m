function numLines = fnumlines(fid)
% finds number of lines in file and rewinds file to start state
% does not fclose
%
    numLines=0;
    while ~feof(fid)
        fgetl(fid);
        numLines=numLines+1;
    end
    frewind(fid);
end