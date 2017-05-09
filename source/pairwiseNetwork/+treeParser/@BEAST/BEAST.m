classdef BEAST < treeParser.baseClass

    properties (SetAccess = protected)
       leafNames = {}
       nodeNames = {}
       nodeAttributes = {}
    end
    
    properties (Hidden = true)
        leafNameTranslation
        fields
    end

    methods % BEAST class default constructor
        function obj = BEAST(varargin)
            if nargin>0
                Fields={};
                for pair = reshape(varargin,2,[])
                    obj.(pair{1})=pair{2}; 
                    Fields{end+1}=pair{1};
                end
                obj.fields=Fields;
                
                checkErrors(obj)
            end
        end           
    end
    
    methods 
        obj = loadTreeAttributes(obj,fidString,lastTipDate,lastTipDateUnits)   
        [ht, phyData] = plot(obj,varargin);
    end
   
    methods
        function obj = addNodeAttribute(obj,fieldName,data,varargin)
            options=struct('treeIdx',1);
            inputs=varargin;
            options = keyValuePairVararginHandler(options, inputs);
            
            if nargin==3,
                warning('no treeIdx specified: assuming data belongs to first tree only')
            end
                        
            if ~isfield(obj.nodeAttributes{options.treeIdx},fieldName)
                obj.nodeAttributes{options.treeIdx}.(fieldName)=data;
            else
                warning('attribute already exists: updating with new values and backing up old in _old')
                obj.nodeAttributes{options.treeIdx}.([fieldName,'_old'])=obj.nodeAttributes{options.treeIdx}.(fieldName);
                obj.nodeAttributes{options.treeIdx}.(fieldName)=data;
            end
        end
    end
    
    methods  %save and load struct
        function saveobj(obj,fileNameString)
            S = saveobj@treeParser.baseClass(obj,fileNameString);
            S.fields=obj.fields;
            for F=obj.fields,
                S.(F{:}) = obj.(F{:});
                save(fileNameString,'-append','S');
            end
        end
        
        function obj = reload(obj,fileNameString)
            obj = reload@treeParser.baseClass(obj,fileNameString);
            load(fileNameString,'S');
            obj.fields=S.fields;
            for F=obj.fields,
                obj.(F{:}) = S.(F{:});
            end
        end
    end
    
    methods (Static) % load with class type
        function obj = loadobj(fileNameString)
            obj = treeParser.BEAST;
            obj = reload(obj,fileNameString);
        end
    end
end