function [options, selected] = plotOptionHandler(obj,varargin)

options=struct('includeOnlyNodes',[],'includeOnlyLinks',[],'timeInterval',[],'timeUnits','years' ...
        ,'duration',[],'slice',[0,1],'randSeed',100,'displayNodes',true ...
        ,'sourceAdminL0',[],'sourceAdminL1',[],'sourceAdminL2',[],'sourceAdminL3',[],'sourceClusterID',[], 'receiveAdminL1',[],'receiveAdminL2',[],'receiveClusterID',[],'receiveAdminL0',[],'receiveAdminL3',[] ...
        ,'forceAdminLevel',false,'upAdminLevel',40 ...
        ,'dateGradient',false,'durationGradient',false,'linewidth',2,'linecolormap',colormap(winter(128)),'markersize',8,'markercolormap',[],'marker','o' ...
        ,'mapColor',false,'newFigure',true,'isolateDateUnits','days','uniqueNames',false,'map',true);
options=keyValuePairVararginHandler(options,varargin);

% options
linkIdx=true(obj.numLinks,1);
nodeIdx=true(obj.numNodes,1);
selected.title='';

if options.dateGradient && options.durationGradient
    error('must specify either dateGradient or durationGradient only')
end

% includeOnly logical or index array
if ~isempty(options.includeOnlyNodes)
    if islogical(options.includeOnlyNodes)
        nodeIdx=nodeIdx & options.includeOnlyNodes;
    elseif isa(options.includeOnlyNodes,'double')
        tmp=false(obj.numNodes,1);
        tmp(options.includeOnlyNodes)=true;
        nodeIdx=nodeIdx & tmp;
        clear tmp;
    else
        error('includeOnlyNodes must be logical or index array')
    end
    linkIdx= linkIdx & ismember(obj.network(:,1),find(nodeIdx)) & ismember(obj.network(:,2),find(nodeIdx));
end
if ~isempty(options.includeOnlyLinks)
    if islogical(options.includeOnlyLinks)
        linkIdx=linkIdx & options.includeOnlyLinks;
    elseif isa(options.includeOnlyLinks,'double')
        tmp=false(obj.numLinks,1);
        tmp(options.includeOnlyLinks)=true;
        linkIdx=linkIdx & tmp;
        clear tmp;
    else
        error('includeOnlyLinks must be logical or index array')
    end
    tmp=false(size(nodeIdx));
    tmp(unique(obj.network(linkIdx,:)))=true;
    nodeIdx = nodeIdx & tmp;
    clear tmp;
end


% select by case time intervals
if ~isempty(options.timeInterval)
    titleTimes=options.timeInterval;
    if strcmp(options.timeUnits,'years')
        options.timeInterval=datenumCleaner(options.timeInterval,'years');
    end
    if numel(options.timeInterval)==2  % all cases in range
        nodeIdx= nodeIdx & obj.time>options.timeInterval(1) & obj.time<=options.timeInterval(2);
        linkIdx= linkIdx & ismember(obj.network(:,1),find(nodeIdx)) & ismember(obj.network(:,2),find(nodeIdx));
        selected.title=['Cases from ',num2str(titleTimes(1)),' to ',num2str(titleTimes(2))];
    elseif numel(options.timeInterval)==4 % starting interval and ending interval
        nodeIdx1=obj.time>options.timeInterval(1) & obj.time<=options.timeInterval(2);
        nodeIdx2=obj.time>options.timeInterval(3) & obj.time<=options.timeInterval(4);
        nodeIdx = nodeIdx & (nodeIdx1 | nodeIdx2);
        linkIdx= linkIdx & ismember(obj.network(:,1),find(nodeIdx1)) & ismember(obj.network(:,2),find(nodeIdx2));
        selected.title={['Source cases from ',num2str(titleTimes(1)),' to ',num2str(titleTimes(2))] ...
                        ['Recipient cases from ',num2str(titleTimes(3)),' to ',num2str(titleTimes(4))],};
    end
end

% exclude links with durations above cutoff
if ~isempty(options.duration)
    if strcmp(options.timeUnits,'years')
        options.duration=datenumCleaner(options.duration,'years');
    end
    duration=obj.time(obj.network(:,2))-obj.time(obj.network(:,1));
    if length(options.duration)==1
        linkIdx = linkIdx & (duration <= options.duration); 
    elseif length(options.duration)==2
        linkIdx = linkIdx & (duration > options.duration(1) & duration <= options.duration(2) );
    else
        error('duration option takes either 1 value (upper cutoff) or 2 (interval)');
    end
    tmp=false(size(nodeIdx));
    tmp(unique(obj.network(linkIdx,:)))=true;
    nodeIdx = nodeIdx & tmp;
    clear tmp;
end

% locations
sourceGeoNodeIdx=true(size(nodeIdx));
receiveGeoNodeIdx=true(size(nodeIdx));
if ~isempty(options.sourceAdminL0)
    if ischar(options.sourceAdminL0)
        options.sourceAdminL0={options.sourceAdminL0};
    end
    sourceGeoNodeIdx = sourceGeoNodeIdx & ismember(lower(obj.location.AdminL0),lower(options.sourceAdminL0));    
end
if ~isempty(options.sourceAdminL1)
    if ischar(options.sourceAdminL1)
        options.sourceAdminL1={options.sourceAdminL1};
    end
    sourceGeoNodeIdx = sourceGeoNodeIdx & ismember(lower(obj.location.AdminL1),lower(options.sourceAdminL1));    
end
if ~isempty(options.sourceAdminL2)
    if ischar(options.sourceAdminL2)
        options.sourceAdminL2={options.sourceAdminL2};
    end
    sourceGeoNodeIdx = sourceGeoNodeIdx & ismember(lower(obj.location.AdminL2),lower(options.sourceAdminL2));
end
if ~isempty(options.sourceAdminL3)
    if ischar(options.sourceAdminL3)
        options.sourceAdminL3={options.sourceAdminL3};
    end
    sourceGeoNodeIdx = sourceGeoNodeIdx & ismember(lower(obj.location.AdminL2),lower(options.sourceAdminL3));
end
if ~isempty(options.sourceClusterID)
    if ischar(options.sourceClusterID)
        options.sourceClusterID={options.sourceClusterID};
    end
    sourceGeoNodeIdx = sourceGeoNodeIdx & ismember(lower(obj.location.ClusterID),lower(options.sourceClusterID));
end
if ~isempty(options.receiveAdminL0)
    if ischar(options.receiveAdminL0)
        options.receiveAdminL0={options.receiveAdminL0};
    end
    receiveGeoNodeIdx = receiveGeoNodeIdx & ismember(lower(obj.location.AdminL0),lower(options.receiveAdminL0));    
end
if ~isempty(options.receiveAdminL1)
    if ischar(options.receiveAdminL1)
        options.receiveAdminL1={options.receiveAdminL1};
    end
    receiveGeoNodeIdx = receiveGeoNodeIdx & ismember(lower(obj.location.AdminL1),lower(options.receiveAdminL1));    
end
if ~isempty(options.receiveAdminL2)
    if ischar(options.receiveAdminL2)
        options.receiveAdminL2={options.receiveAdminL2};
    end
    receiveGeoNodeIdx = receiveGeoNodeIdx & ismember(lower(obj.location.AdminL2),lower(options.receiveAdminL2));    
end
if ~isempty(options.receiveAdminL3)
    if ischar(options.receiveAdminL3)
        options.receiveAdminL3={options.receiveAdminL3};
    end
    receiveGeoNodeIdx = receiveGeoNodeIdx & ismember(lower(obj.location.AdminL2),lower(options.receiveAdminL3));    
end
if ~isempty(options.receiveClusterID)
    if ischar(options.receiveClusterID)
        options.receiveClusterID={options.receiveClusterID};
    end
    receiveGeoNodeIdx = receiveGeoNodeIdx & ismember(lower(obj.location.ClusterID),lower(options.receiveClusterID));    
end


linkIdx= linkIdx & ismember(obj.network(:,1),find(nodeIdx & sourceGeoNodeIdx)) & ismember(obj.network(:,2),find(nodeIdx & receiveGeoNodeIdx));

% clean up output
selected.network=obj.network(linkIdx,:);
selected.numLinks=sum(linkIdx);
% nodeIdx=unique(selected.network);
nodeIdx=find(nodeIdx);
selected.numNodes=length(nodeIdx);
selected.nodes=nodeIdx;

if length(obj.location.AdminL0)>1,
    selected.location.AdminL0=obj.location.AdminL0(nodeIdx);
end
if length(obj.location.AdminL1)>1,
    selected.location.AdminL1=obj.location.AdminL1(nodeIdx);
end
if length(obj.location.AdminL2)>1,
    selected.location.AdminL2=obj.location.AdminL2(nodeIdx);
end
if length(obj.location.AdminL3)>1,
    selected.location.AdminL3=obj.location.AdminL3(nodeIdx);
end
if isfield(obj.location,'ClusterID') && length(obj.location.ClusterID)>1,
    selected.location.ClusterID=obj.location.ClusterID(nodeIdx);
end
selected.time=obj.time(nodeIdx);
selected.nodeNames=obj.nodeNames(nodeIdx);


end