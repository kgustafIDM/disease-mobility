function obj = addAttributesTwoColumnFile(obj,fileId,nameField,attributeField,varargin)
% WARNING: This function currently only works with two column tab-separated input files.

% varargin handling
    options=struct('overwrite',false,'isolateDateUnits','days');
    if nargin>4
        inputs=varargin;
        options = keyValuePairVararginHandler(options,inputs);
    end

dataIn=importdata(fileId,'\t');
obj.sourceFile = [obj.sourceFile,'; ',fileId];

% overwrite?
    if ~options.overwrite && ~isempty(obj.(attributeField))
        err=input('attribute field is not empty, do you want to overwrite (y/n)?  ','s');
        switch err
            case 'n'
                error('attribute field is not empty');
            case 'y'
                display('overwriting attribute field');
        end
    end

% pre-declare
    if any(strcmp(attributeField,obj.doubleFields)),
        obj.(attributeField)=nan(obj.numSeq,1);
    else
        obj.(attributeField)=cell(obj.numSeq,1);
    end
    
% data type and fixing wonky importdata problems
    switch class(dataIn)
        % if the whole thing is text, importdata seems to be unable to
        % split the columns, even with a '\t' tag
        case 'cell'
            if size(dataIn,2)==1,
                tmp=get_cell_tokens(dataIn,'\t');
                dataIn=struct('data',{''},'textdata',{''});
                for k=1:length(tmp)
                    dataIn.textdata{k}=tmp{k}{1};
                    dataIn.data{k}=tmp{k}{2};
                end
            end
        case 'struct'
            for k=1:length(dataIn.textdata),
                % remove second column of blank textdata
                    dataIn.textdata{k}=dataIn.textdata{k,1};
                % header filler
                    if length(dataIn.textdata)==(length(dataIn.data)+1)
                        dataIn.data = [nan; dataIn.data];
                    end
            end
    end
    
% datenums must go out!
    if strcmp(attributeField,'isolateDate')
        if any(isnan(dataIn.data))
            idx=~isnan(dataIn.data);
            dataIn.data(idx)=datenumCleaner(dataIn.data(idx),options.isolateDateUnits);
        else
            dataIn.data=datenumCleaner(dataIn.data,options.isolateDateUnits);
        end
    end
    
% find input indices by aligning fieldName to the input data
    NameIdx=[];
    for k=1:length(dataIn.textdata),
        %tmp=find(~cellfun(@isempty,strfind(obj.(nameField),dataIn.textdata{k})));
        tmp=find(strcmpi(dataIn.textdata{k},obj.(nameField)));
        if isempty(tmp)
            NameIdx=[NameIdx nan];
        else
            NameIdx=[NameIdx tmp];
        end
    end
    
% put data
    for k=1:length(dataIn.textdata),
        if ~isnan(NameIdx(k))
            switch class(dataIn.data(k))
                case 'cell'
                    if any(strcmp(attributeField,obj.doubleFields)),
                        obj.(attributeField)(NameIdx(k))=dataIn.data{k};
                    else
                        obj.(attributeField){NameIdx(k)}=dataIn.data{k};
                    end
                otherwise
                    if any(strcmp(attributeField,obj.doubleFields)),
                        obj.(attributeField)(NameIdx(k))=dataIn.data(k);
                    else
                        obj.(attributeField){NameIdx(k)}=dataIn.data(k);
                    end
            end
        end
    end
    
end