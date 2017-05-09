classdef sequenceData % < handle
    % example call of generic constructor (not all fields called)
    %   SD=sequenceData('sequence',{'gatc'},'isolateDate',[datenumCleaner(now,'days')],'FastaName',{'seq1'},'EPID',{'seqID'},'AdminL0',{'USA'},'AdminL1',{'Washington'},'AdminL2',{'Seattle'},'strain',{'LDPV4'})
    
    
    properties (SetAccess = protected)
        sequence = {}       % cell of sequence strings
        isolateDate = []    % vector of isolateDates as integer datenums
        EPID = {}           % cell of case EPID
        FastaName = {}      % cell of FastaNames
        AdminL0 = {}        % cell of location string
        AdminL1 = {}
        AdminL2 = {}
        AdminL3 = {}
        ClusterID = {}
        strain = {}           % cell of strings of WPV1 or WPV3 or VDPV2        
        phylogeneticDomain = {} % cell with 'VP1' or otherwise        
        sourceFile 
        constructorQuery
        attributes
    end

    properties (SetAccess = protected, Hidden = true)
        shapeFileIdAdminL0 = []  % these are hidden because the scientist shouldn't have to know about them.
        shapeFileIdAdminL1 = []  % all functions are written to know that if numbers are passed to adminLevel
        shapeFileIdAdminL2 = []  % commands, shapefiles are intended.
        cellFields={'sequence','EPID','FastaName','AdminL0','AdminL1','AdminL2','strain','phylogeneticDomain','ClusterID'}
        doubleFields={'isolateDate','shapeFileIdAdminL0','shapeFileIdAdminL1','shapeFileIdAdminL1'};        
    end
    
    properties (Dependent = true)
        numSeq
    end
    
    methods % default constructor
        function obj = sequenceData(varargin)
            if nargin>0, 
                for pair = reshape(varargin,2,[])
                    inpName=pair{1};
                    if any(strcmp(inpName,fieldnames(obj))) && ~any(strcmp(inpName,{'cellFields','doubleFields'}))
                        obj.(inpName)=pair{2}; 
                    else
                        error(['input name is not a valid class property: ',inpName])
                    end
                end
                %checkErrors(obj)
            end
        end
    end
    
    methods (Static) % additional constructors
        obj = constructFromDB(afp,Lineages,strain,strainIdx,varargin)
        % combine from mulitple sources
            obj = combine(varargin)
        % synthetic sequence data
            obj = generateSyntheticAlignmentOnTree(treeParserObj,varargin)
        
        function obj = bootstrap(objIn)
            obj=objIn;
            bootstrapIdx=randi(length(objIn.sequence{1}),length(objIn.sequence{1}),1);
            for m=1:length(obj.sequence)
                obj.sequence{m}=obj.sequence{m}(bootstrapIdx);
            end
        end
    end
    
    methods
        checkErrors(obj)
        obj = constructFromFasta(obj,fileId,varargin)
        obj = addAttributesTwoColumnFile(obj,fileId,varargin)
        obj = addAttributesStruct(obj,attributesStruct)
        exportNexusSequence(obj,fileId,varargin)
        % export to BEAST
            exportBEASTXML(obj,fileId,varargin)
            exportBEASTTabDelimitedAttribute(obj,fileId,attributeField,varargin)
    end
    
    methods % dependent getters
        function NumSeq = get.numSeq(obj)
                NumSeq = length(obj.sequence);
        end
    end
    
    methods %setters
        function obj=setStrain(obj,strainStr)
            obj.strain=cell(obj.numSeq,1);
            if obj.numSeq>0
                for k=1:obj.numSeq,
                    obj.strain{k}=strainStr;
                end
            end
        end
        function obj=setPhylogeneticDomain(obj,domainStr)
            obj.phylogeneticDomain=cell(obj.numSeq,1);
            if obj.numSeq>0
                for k=1:obj.numSeq,
                    obj.phylogeneticDomain{k}=domainStr;
                end
            end
        end
    end
    
    methods  %save 
        function S = saveobj(obj,fileNameString)
            for F=reshape(fieldnames(obj),1,[])
                S.(F{:})=obj.(F{:});
            end
            save(fileNameString,'S');
        end
        
        function obj = reload(obj,fileNameString)
            load(fileNameString,'S');
            for F=reshape(fieldnames(S),1,[])
                if ~strcmp(F{:},'numSeq') %don't try to load dependent properties
                    obj.(F{:})=S.(F{:});
                end
            end
        end
    end
    
    methods (Static) % load with class type
        function obj = loadobj(fileNameString)
            obj = sequenceData;
            obj = reload(obj,fileNameString);
        end
    end
end