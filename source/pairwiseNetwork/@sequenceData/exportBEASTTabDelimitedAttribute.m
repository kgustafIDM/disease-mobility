function exportBEASTTabDelimitedAttribute(obj,fileId,attributeField,varargin)
%INPUT: sequenceData obj, export fileId, attributeField
%   varargin: leaf tip naming convention
%       'CaseID':: default 'FastaName'; allowed 'EPID'
%
%OUTPUT: tab delimited two column attribute file at fileId
% 

% vararing handling
    options=struct('CaseID','FastaName','isolateDateUnits','days');
    inputs=varargin;
    options = keyValuePairVararginHandler(options,inputs);

fid=fopen(fileId,'w');

fprintf(fid, ['traits\t',attributeField,'\n']);
for k=1:obj.numSeq,
    switch class(obj.(attributeField))
        case 'cell'
            tmp=obj.(attributeField){k};
        otherwise
            tmp=obj.(attributeField)(k);
            if strcmp(attributeField,'isolateDate') && strcmp(options.isolateDateUnits,'years')
                tmp=datenum2years(tmp);
            end
            tmp=num2str(tmp);
    end
    if k<obj.numSeq,
        fprintf(fid,'%s\t%s\n',obj.(options.CaseID){k},tmp);
    else
        fprintf(fid,'%s\t%s',obj.(options.CaseID){k},tmp);
    end
end
fclose(fid);


end

