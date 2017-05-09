classdef traitVisData < handle
    
    properties
        type = {} % location, other trait
        textAttribute = {}  
        numberAttribute = [];
        shapefileID = [];
        cmap = [];
        otherAttribute = {};
        otherAttributeType = {};
        parentType = {};
        parentTextAttribute = {};
        parentNumberAttribute = [];
        parentShapefileID= [];
        parentOtherAttribute = {};
        parentOtherAttributeType = {};        
    end
    
    methods % constructor
        function obj = traitVisData(varargin)
            if nargin>0
                options=struct('type','text','textAttribute',false,'numberAttribute',false ...
                            ,'shapefileID',false,'cmap',false,'otherAttribute',false ...
                            ,'otherAttributeType',false,'parentTextAttribute',false ...
                            ,'parentNumberAttribute',false,'parentShapefileID',false ...
                            ,'parentOtherAttribute',false,'parentOtherAttributeType',false ...
                            ,'parentType',false);
                inputs=varargin;
                options=keyValuePairVararginHandler(options,inputs);

                for F=reshape(fieldnames(options),1,[])
                    if ~islogical(options.(F{:}))
                        obj.(F{:})=options.(F{:});
                    end
                end
            end
        end
    end
        
    methods
        function obj = addAttribute(obj, varargin)
                inputs=varargin;
                options=struct('type',false,'textAttribute',false,'numberAttribute',false ...
                            ,'shapefileID',false,'cmap',false,'otherAttribute',false ...
                            ,'otherAttributeType',false,'parentTextAttribute',false ...
                            ,'parentNumberAttribute',false,'parentShapefileID',false ...
                            ,'parentOtherAttribute',false,'parentOtherAttributeType',false ...
                            ,'parentType',false);
                options=keyValuePairVararginHandler(options,inputs);
                for F=reshape(fieldnames(options),1,[])
                    if ~islogical(options.(F{:}))
                        if isempty(obj.(F{:}))
                            obj.(F{:})=options.(F{:});
                        else
                            val = input([F{:},' is not empty. Do you want to overwrite? y/n: '],'s');
                            if strcmpi(val,'y')
                                obj.(F{:})=options.(F{:});
                            else
                                display('ignoring field.')
                            end
                        end
                    end
                end
        end
        obj = createColormap(obj,center,varargin)
        plotMap(obj,varargin)
    end
    
    methods (Static)
        function objOut = getSubset(obj,input)
            objOut=traitVisData;
            if ~isa(input,'cell')
                input={input};
            end
            
            idx=[];
            for m=1:length(input)
                if isa(input{m},'char')
                    for k=1:length(obj.textAttribute), 
                        A=find(regexpi(obj.textAttribute{k},input{m}), 1); 
                        if ~isempty(A) && ~ismember(k,idx), 
                            idx = [idx k]; 
                        end
                    end
                elseif isa(input,'double')
                    for k=1:length(obj.numberAttribute), 
                        A=find(obj.numberAttribute{k}==input{m},1); 
                        if ~isempty(A) && ~ismember(k,idx), 
                            idx = [idx k]; 
                        end
                    end
                else
                    error('input type must be textAttribute or numberAttribute')
                end
                
                idx=sort(idx);
                objOut.type=obj.type;
                objOut.otherAttributeType=obj.otherAttributeType;
                if isstruct(obj.otherAttribute)
                    objOut.otherAttribute=obj.otherAttribute(idx);
                end
                for k=1:length(idx),
                    objOut.textAttribute{end+1}=obj.textAttribute{idx(k)};
                    objOut.numberAttribute(end+1)=obj.numberAttribute(idx(k));
                    objOut.cmap(end+1,:)=obj.cmap(idx(k),:);
                    if ~isstruct(obj.otherAttribute)
                        objOut.otherAttribute(end+1)=obj.otherAttribute(idx(k));
                    end
                end
            end            
        end
    end
    
    methods  %save and load struct
        function S = saveobj(obj,fileNameString)
            for F=reshape(fieldnames(obj),1,[])
                S.(F{:})=obj.(F{:});
            end
            save(fileNameString,'S');
        end
        
        function obj = reload(obj,fileNameString)
            load(fileNameString,'S');
            for F=reshape(fieldnames(S),1,[])
                obj.(F{:})=S.(F{:});
            end
        end
    end
    
    methods (Static) % load with class type
        function obj = loadobj(fileNameString)
            obj = traitVisData;
            obj = reload(obj,fileNameString);
        end
    end
end