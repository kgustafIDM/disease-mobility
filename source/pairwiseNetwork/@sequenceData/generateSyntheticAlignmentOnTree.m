function obj = generateSyntheticAlignmentOnTree(treeParserObj,varargin)
% works, but is delicate!!!
%
% fix name importation '''' multiple quote bug

% varargin handling
    options=struct('dateRange',false,'timeUnits','years','sequenceDataParent',false,'genomeLength',1,'mutationModel','jukesCantor','mutationModelParams',struct(),'clockModel','strict','clockModelParams',struct('meanRate',0.01/365.2431,'stdRate',0.0005/365.2431),'treeIdx',treeParserObj.numTrees);
    if nargin>1
        inputs=varargin;
        options = keyValuePairVararginHandler(options,inputs);
    end

    if isa(options.sequenceDataParent,'sequenceData') && options.genomeLength==1,
        options.genomeLength=length(options.sequenceDataParent.sequence{1});
    else
        options.genomeLength=906;
    end
    
% declare obj
    obj=sequenceData;

% source
    obj.sourceFile=treeParserObj.sourceFile;
    
% domain
    obj.phylogeneticDomain=cell(treeParserObj.numLeaves,1);
    obj.strain=cell(treeParserObj.numLeaves,1);
    for k=1:treeParserObj.numLeaves,
        obj.phylogeneticDomain{k}='synthetic data';
        obj.strain{k}=options.mutationModel;
    end
            
% constructor query
    obj.constructorQuery={'synethtic data from treeParserObj'};
    for k=1:length(varargin)
        obj.constructorQuery{end+1}=varargin{k};
    end

% names in ascending order on tree as parsed from the root
    numLeaves=sum(~treeParserObj.isCoalescentTime{options.treeIdx});
    obj.FastaName=cell(numLeaves,1);
    obj.EPID=cell(numLeaves,1);
    leafIdx=find(~treeParserObj.isCoalescentTime{options.treeIdx});
    for k=1:numLeaves
        obj.FastaName{k}=treeParserObj.nodeNames{options.treeIdx}{leafIdx(k)};
        obj.EPID{k}=treeParserObj.nodeNames{options.treeIdx}{leafIdx(k)};
    end
    
% dates in ascending order on tree as parsed from the root
    obj.isolateDate=treeParserObj.nodeTimes{options.treeIdx}(~treeParserObj.isCoalescentTime{options.treeIdx});

% inherit additonal info from sequenceDataParent if given
    if ~islogical(options.sequenceDataParent) && isa(options.sequenceDataParent,'sequenceData')
        obj.attributes=options.sequenceDataParent.attributes;
        % sort on names
            [isName, nameIdx] = ismember(obj.FastaName,options.sequenceDataParent.FastaName);
            if all(~isName)
                [isName, nameIdx] = ismember(obj.EPID,options.sequenceDataParent.EPID);
            end
            if all(isName)
                for k=1:length(isName)
                    fields={'AdminL0','AdminL1','AdminL2','shapeFileIdAdminL0','shapeFileIdAdminL1','shapeFileIdAdminL2','EPID','FastaName'};
                    for F=fields
                        if ~isempty(options.sequenceDataParent.(F{:}))
                            obj.(F{:})=options.sequenceDataParent.(F{:})(nameIdx);
                        end
                    end
                end
            else
                warning('sequenceDataParent names do not match all tree leaf names. Ignoring sequenceDataParent')
            end
    end
    
% pre-declare output sequence
    obj.sequence=cell(treeParserObj.numLeaves,1);
    
% generate fake data
    obj = syntheticAlignment(obj,treeParserObj,options);
    
% convert dates to years
    if strcmpi('years',options.timeUnits)
        obj.isolateDate=datenum2years(obj.isolateDate);
    end
    
% clean up    
    checkErrors(obj)
end


function obj = syntheticAlignment(obj,treeParserObj,options)

tree=treeParserObj.trees{options.treeIdx};

numNodes=get(tree,'numnodes');
nodeTimes=getPhyloDistancesFromRoot(tree)*365.2431;
nodeNames=get(tree,'nodeNames');
    
sequence=cell(numNodes,1);

D=zeros(numNodes,1);

sequence{numNodes}=randseq(options.genomeLength);

for k=(numNodes:-1:1);
    iChild=findProgenyNodes(tree,k,'recursive',false);
    if length(iChild)>1
        for c=2:3
            D(iChild(c))=nodeTimes(iChild(c))-nodeTimes(k);
            switch options.mutationModel
                case 'jukesCantor'
                    sequence{iChild(c)}=mutateJukesCantor(sequence{k},D(iChild(c))*options.clockModelParams.meanRate);
            end
        end
    end
end

% sort on names
    [isName, nameIdx] = ismember(obj.FastaName,nodeNames);
    if all(~isName)
        [isName, nameIdx] = ismember(obj.EPID,nodeNames);
    end
    nameIdx=nameIdx(isName);

    for k=1:treeParserObj.numLeaves,
        obj.sequence{k}=sequence{nameIdx(k)};
    end

end

function seqOut = mutateJukesCantor(seqIn,Dist)

num=round(length(seqIn)*3/4*(1-exp(-4/3*Dist)));
randIdx=randperm(length(seqIn),num);
seqOut=seqIn;
seqOut(randIdx)=randseq(num);

end
