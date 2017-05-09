classdef BEAST < network.baseClass
    % this class stores case directed network data 
    
    properties (SetAccess = protected)
        sourceData   % sequenceData object corresponding to case-directed Network
        location = struct('AdminL0',{''},'AdminL1',{''},'AdminL2',{''}) % struct with names and [N x 2] location array (long, lat)
        time = []                         % times for each node
        alpha = []  % upper alpha level for significance
    end
    
    methods % network.caseDirected class default constructor
        function obj = BEAST(varargin)
            if nargin>0,
                options=struct('treeIdx',1,'alpha',0.05,'FDR',true,'adminLevel','AdminL1' ...
                            ,'alternateLeafNames',[]);
                if nargin>1,
                    options=keyValuePairVararginHandler(options,varargin(2:end));
                end

                obj.sourceData=varargin{1};
                if isa(obj.sourceData,'treeParser.BEAST') 
                    % first pass through
                    obj.nodeNames = obj.sourceData.nodeNames{options.treeIdx};
                    obj.time = obj.sourceData.nodeTimes{options.treeIdx};
                    for P=reshape(fieldnames(obj.location),1,[])
                        if isfield(obj.sourceData.nodeAttributes{options.treeIdx},P{:})
                            obj.location.(P{:})=obj.sourceData.nodeAttributes{options.treeIdx}.(P{:});
                        end
                    end
                    
                    % significant nodes
                    tmpLeafNames=obj.sourceData.leafNames;
                    leafIdx = ismember(obj.nodeNames,tmpLeafNames);
                    if ~any(leafIdx)  % figtree does dumb shit to the names
                        obj.nodeNames=regexprep(obj.nodeNames,'''','');
                        tmpLeafNames=get_cell_tokens(tmpLeafNames,'\[');
                        tmpLeafNames=[tmpLeafNames{:}]';
                        leafIdx = ismember(obj.nodeNames,tmpLeafNames(:,1));
                    end
                    
                    switch options.adminLevel
                        case ''
                            nodeProb= obj.sourceData.nodeAttributes{options.treeIdx}.posterior;
                        case 'AdminL1'
                            nodeProb= obj.sourceData.nodeAttributes{options.treeIdx}.posterior.* ...
                                            obj.sourceData.nodeAttributes{options.treeIdx}.AdminL1_prob;
                        case 'AdminL2'
                            nodeProb= obj.sourceData.nodeAttributes{options.treeIdx}.AdminL2_prob;% .* ...
%                                             obj.sourceData.nodeAttributes{options.treeIdx}.posterior;
                    end
                    
                    % links
                    P=sort(nodeProb);
                    if options.FDR
                        Alpha=(1:length(nodeProb))'*options.alpha/length(nodeProb);
                    else
                        Alpha=options.alpha;
                    end
                    cutoff=find(P>=(1-Alpha),1,'first');
                    obj.p=1-P(cutoff);
                    obj.alpha=P(cutoff);

                    if ~isempty(P(cutoff))
                        nodeIdx = leafIdx | (nodeProb >=P(cutoff));% & cellfun(@isempty,regexp(obj.location.(options.adminLevel),'+')));
                    else
                        display('No internal nodes are significant, and so there is no network.');
                        return
                    end
                    
                    % cleanup nodes
                    obj.nodeNames=obj.nodeNames(nodeIdx);
                    obj.time=obj.time(nodeIdx);
                    for P=reshape(fieldnames(obj.location),1,[])
                        if isfield(obj.sourceData.nodeAttributes{options.treeIdx},P{:})
                            obj.location.(P{:})=obj.location.(P{:})(nodeIdx);
                        end
                    end
                    
                    % find links
                    leafIdx = ismember(obj.nodeNames,tmpLeafNames);
                    tree=obj.sourceData.trees{options.treeIdx};
                    nodeNamesTreeOrder=regexprep(get(tree,'nodenames'),'''','');
					
                    [ix,nodeObjOrderToTreeOrder] = ismember(nodeNamesTreeOrder,obj.nodeNames);
                    nodeObjOrderToTreeOrder=nodeObjOrderToTreeOrder(ix);
                                        
 					root = findNodeCommonToLeaves(tree,find(leafIdx(nodeObjOrderToTreeOrder)));
                    nodesUnderRoot = findProgenyNodes(tree, root);
                    nodeNamesFromRoot=nodeNamesTreeOrder(nodesUnderRoot);
                    [ix, loc] = ismember(nodeNamesFromRoot,obj.nodeNames);
                    eligibleNodesFromRoot=find(ix);
                    loc=loc(ix);
                    
                    network=nan(2*obj.sourceData.numLeaves-1,2);
                    
                    count=0;
                    for k=1:length(eligibleNodesFromRoot)
                        parent=find(strcmp(nodeNamesFromRoot(eligibleNodesFromRoot(k)),nodeNamesTreeOrder));
                        children = findProgenyNodes(tree,parent);
                        childrenNames = nodeNamesTreeOrder(children);
                        [ixP, locP] = ismember(childrenNames,obj.nodeNames);
                        eligibleNodesUnderParent=find(ixP);
                        locP=locP(ixP);
                        for m=2:length(eligibleNodesUnderParent),
                            count=count+1;
                            network(count,:)=[loc(k),locP(m)];
                        end
                    end
                    network=network(~isnan(network(:,1)),:);
                    
                    % children can only have 1 parent in a tree
                    [~, IA]=unique(network(:,2),'last');
                    obj.network=network(IA,:);
                    
					% relabel names if necessary
                    if ~isempty(options.alternateLeafNames)
                        [ix, loc]=ismember(obj.nodeNames,options.alternateLeafNames(:,1));
                        obj.nodeNames(ix)=options.alternateLeafNames(loc(ix),2);
                    end
                  
                    % re-packaging
                    [~,sortIdx]=sort(obj.network(:,1));
                    obj.network=obj.network(sortIdx,:);

                    [nodes,~,nodeLoc]=unique(obj.network);
                    obj.time=obj.time(nodes);
                    obj.nodeNames=obj.nodeNames(nodes);
                    for P=reshape(fieldnames(obj.location),1,[])
                        if numel(obj.location.(P{:}))>1
                            obj.location.(P{:})=obj.location.(P{:})(nodes);
                        end
                    end

                    tmp=(1:length(nodes))';
                    obj.network=reshape(tmp(nodeLoc),[],2);
                    
                else
                    error('network:BEAST:wrongInputClass','Input must be treeParser.BEAST object.')
                end


                %checkErrors(obj)
            end
        end
    end
    
    methods %plot
        function [he, hv]=plotMap(obj,varargin)
            [he, hv]=network.plotMap(obj,varargin{:});
        end
        function [he, hv]=plotTimeseries(obj,varargin)
            [he, hv]=network.plotTimeseries(obj,varargin{:});
        end
        function [options, selected] =plotOptionHandler(obj,varargin)
            [options, selected] =network.plotOptionHandler(obj,varargin{:});
        end
    end

    
    
    methods  %save and load struct
        function saveobj(obj,fileNameString)
            S = saveobj@network.baseClass(obj,fileNameString);
            for F=fieldnames(obj)',
                if ~ismember(F{:},{'sourceData','numNodes','numLinks'})
                    S.(F{:}) = obj.(F{:});
                end
                save(fileNameString,'-append','S');
            end
        end
        
        function obj = reload(obj,fileNameString)
            obj = reload@network.baseClass(obj,fileNameString);
            load(fileNameString,'S');
            for F=fieldnames(S)',
               % if ~ismember(F{:},{'numNodes','numLinks'})
                    obj.(F{:}) = S.(F{:});
               % end
            end
        end
    end
    
    methods (Static) % load with class type
        function obj = loadobj(fileNameString)
            obj = network.caseDirected;
            obj = reload(obj,fileNameString);
        end
    end
end