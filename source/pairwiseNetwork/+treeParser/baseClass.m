classdef baseClass < handle
    % this class stores tree data needed for coalescent likelihood calculations
    % is currently compatable only with simpleCoalecentLogLikelihood
    
    properties (SetAccess = protected)
        sourceFile
        timeUnits = ''
        nodeTimes = {}       % cell of ascending order vectors of inferred coalescent (interior nodes) and sample times (datenum)
        isCoalescentTime ={} % cell of logical vectors indicating if node time corresponds to a coalescent time and not a sample time
        numLineages          % cell of vectors of number of lineages between each time point
                                % was too slow if I kept it dependent on isCoalescentTime, which is what it should be.
        trees = {}           % cell of matlab phytree objects
    end
    
    properties (Dependent = true)
        leafDates        % leaf dates derived from tree properties
        numLeaves        % number of sequences in tree
        numTrees         % number of trees in ensemble
    end

    methods % baseClass default constructor
        function obj = baseClass(varargin)
            if nargin>0, 
                for pair = reshape(varargin,2,[])
                    inpName=pair{1};
                    if any(strcmp(inpName,fieldnames(obj)))
                        obj.(inpName)=pair{2}; 
                    else
                        error(['input name is not a valid class property: ',inpName])
                    end
                end
                getNumLineages(obj);
                treeParser.checkErrors(obj)
            end
        end
    end
    
    methods % dependent getters
        function NumTrees = get.numTrees(obj)
            NumTrees = length(obj.nodeTimes);
        end
        
        function NumLeaves = get.numLeaves(obj)
            if ~isempty(obj.nodeTimes)
                NumLeaves = (length(obj.nodeTimes{1})+1)/2;
            else
                NumLeaves=0;
            end
        end
        
        function LeafDates=get.leafDates(obj)
            if ~isempty(obj.nodeTimes)
                LeafDates=obj.nodeTimes{1}(~obj.isCoalescentTime{1});
            else
                LeafDates=[];
            end
        end
    end
       
    methods %plot
        [ht, phyData] = plot(obj,varargin);
    end
    
    methods
        function obj = getNumLineages(obj)
            obj.numLineages=cell(length(obj.nodeTimes),1);
            for k=1:length(obj.nodeTimes),
                obj.numLineages{k}=1+cumsum(2*obj.isCoalescentTime{k}-1);
            end
        end
    end
    
    methods  %save 
        function S = saveobj(obj,fileNameString)
            S.nodeTimes = obj.nodeTimes;
            S.isCoalescentTime = obj.isCoalescentTime;
            S.timeUnits = obj.timeUnits;
            S.sourceFile = obj.sourceFile;
            save(fileNameString,'S');
        end
        
        function obj = reload(obj,fileNameString)
            load(fileNameString,'S');
            obj.nodeTimes = S.nodeTimes;
            obj.isCoalescentTime = S.isCoalescentTime;
            obj.timeUnits = S.timeUnits;
            obj.sourceFile = S.sourceFile;
        end
    end
    
    methods (Static) % load with class type
        function obj = loadobj(fileNameString)
            obj = treeParser.baseClass;
            obj = reload(obj,fileNameString);
        end
    end
end