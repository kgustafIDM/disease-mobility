function exportNexusSequence(obj,fileId,varargin)
% inputs: fileId (should be .xml)
%
%   varargin :: key-value pairs
%       'comments' :: 'comment string' possibly cell array for multiline
%       'CaseID' :: 'FastaName' (default) or 'EPID'
%       'dataType' :: 'dna' (default)
%       'gap' :: '-' (default) : gap character
%       'missing' :: 'N' (default) : unknown base character
%

% vararing handling
    options=struct('CaseID','FastaName','dataType','dna','gap','-','missing','N','comments',false);
    inputs=varargin;
    options = keyValuePairVararginHandler(options,inputs);


fid=fopen(fileId,'w');

fprintf(fid,'%s\r', '#nexus');
fprintf(fid,'%s\r','begin data;');
fprintf(fid,'%s\r',['dimensions ntax=',num2str(obj.numSeq),' nchar=',num2str(length(obj.sequence{1})),';']);
fprintf(fid,'%s\r',['format datatype=',options.dataType,' gap=',options.gap,' missing=',options.missing,';']);
fprintf(fid,'%s\r','matrix');

for k=1:obj.numSeq,
    fprintf(fid,'%s\t%s\r',regexprep(obj.(options.CaseID){k},'-','_'),upper(obj.sequence{k}));
end

fprintf(fid,'%s\r',';');
fprintf(fid,'%s\r','end;');

% provenance comment
    fprintf(fid,'%s\r','');
    fprintf(fid,'%s\r','provenance;');
    fprintf(fid,'%s\r',['Source files:: ',obj.sourceFile,';']);
    fprintf(fid,'%s\r',['This file was created on:: ',datestr(now,'dd-mmm-yyyy HH:MM:SS'),';']);
    fprintf(fid,'%s','end;');

% other comments
    if ~islogical(options.comments)
        fprintf(fid,'%s\r','');
        fprintf(fid,'%s\r','');
        fprintf(fid,'%s\r','comments;');
        switch class(options.comments)
            case 'char'
                fprintf(fid,'%s',[options.comments,';']);
            case 'cell'
                k=1;
                while k<length(options.comments)
                    fprintf(fid,'%s\r',[options.comments{k},';']);
                    k=k+1;
                end
                fprintf(fid,'%s',[options.comments{k},';']);
        end
    end

fclose(fid);
end