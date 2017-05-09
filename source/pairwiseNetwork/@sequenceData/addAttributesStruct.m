function obj = addAttributesStruct(obj,attributeStruct,varargin)
% adds attributes to sequenceData object using struct
% struct must be of form:
%       S.(EPID or FastaName or sequence) = {...}
%       S.(attributeName1) = {...} or [...] as appropriate
%       S.(attributeName2) ...
%

% varargin handling
    options=struct('overwrite',false,'isolateDateUnits','days');
    if nargin>4
        inputs=varargin;
        options = keyValuePairVararginHandler(options,inputs);
    end
    
    
% attributeStruct checker
    fields=fieldnames(attributeStruct);
    if ~ismember(fields{1},{'EPID','FastaName','sequence'})
        error('first element of attributeStruct must be EPID, FastaName, or sequence string')
    end
    
    for k=2:length(fields)
        if ~all(ismember(fields{k},fieldnames(obj)))
            error(['attribute ',fields{k},' is not allowed fieldname of sequenceData object'])
        end
    end
  
    lengths=zeros(length(fields),1);
    for k=1:length(fields)
        lengths(k)=length(attributeStruct.(fields{k}));
    end
    if ~all(lengths==lengths(1))
        error('attributes fields must all have same number of elements')
    end
    

% writing attributes
    for k=2:length(fields)
        
        % overwrite?
            if ~options.overwrite && ~isempty(obj.(fields{k}))
                err=input('attribute field is not empty, do you want to overwrite (y/n)?  ','s');
                switch err
                    case 'n'
                        error('attribute field is not empty');
                    case 'y'
                        display('overwriting attribute field');
                end
            end
            
        % pre-declare
            if any(strcmp(fields{k},obj.doubleFields)),
                obj.(fields{k})=nan(obj.numSeq,1);
            else
                obj.(fields{k})=cell(obj.numSeq,1);
            end
            
        % datenums must go out!
            if strcmp(fields{k},'isolateDate')
                if any(isnan(attributeStruct.(fields{k})))
                    idx=~isnan(attributeStruct.(fields{k}));
                    attributeStruct.(fields{k})(idx)=datenumCleaner(attributeStruct.(fields{k})(idx),options.isolateDateUnits);
                else
                    attributeStruct.(fields{k})=datenumCleaner(attributeStruct.(fields{k}),options.isolateDateUnits);
                end
            end

        % find input indices by aligning fieldName to the input data
            NameIdx=[];
            for m=1:length(attributeStruct.(fields{1})),
                %tmp=find(~cellfun(@isempty,strfind(obj.(fields{1}),attributeStruct.(fields{1}){m})));
                tmp=find(strcmpi(attributeStruct.(fields{1}){m},obj.(fields{1})));
                if isempty(tmp)
                    NameIdx=[NameIdx nan];
                else
                    NameIdx=[NameIdx tmp];
                end
            end

        % put data
            for m=1:length(attributeStruct.(fields{1})),
                if ~isnan(NameIdx(m))
                    switch class(attributeStruct.(fields{k})(m))
                        case 'cell'
                            if any(strcmp(fields{k},obj.doubleFields)),
                                obj.(fields{k})(NameIdx(m))=attributeStruct.(fields{k}){m};
                            else
                                obj.(fields{k}){NameIdx(m)}=attributeStruct.(fields{k}){m};
                            end
                        otherwise
                            if any(strcmp(fields{k},obj.doubleFields)),
                                obj.(fields{k})(NameIdx(m))=attributeStruct.(fields{k})(m);
                            else
                                obj.(fields{k}){NameIdx(m)}=attributeStruct.(fields{k})(m);
                            end
                    end
                end
            end
          
    end
    
end