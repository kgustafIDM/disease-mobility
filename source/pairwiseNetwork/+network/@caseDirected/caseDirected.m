classdef caseDirected < network.baseClass
    % this class stores case directed network data 
    
    properties (SetAccess = protected)
        sourceData   % sequenceData object corresponding to case-directed Network
        location = struct('AdminL0',{''},'AdminL1',{''},'AdminL2',{''},'AdminL3',{''},'ClusterID',{''}) % struct with names and [N x 2] location array (long, lat)
        time = []                         % times for each node
        alpha = []  % upper alpha level for significance
        clockViolation=[]
    end
    
    methods % network.caseDirected class default constructor
        function obj = caseDirected(varargin)
            if nargin>0,
                obj.sourceData=varargin{1};
                if isa(obj.sourceData,'sequenceData')
                    options=struct('clockRate',0.01/365.2431,'alpha',0.05,'distanceMethod','Hasegawa','FDR',true);
                    if nargin>1,
                        options=keyValuePairVararginHandler(options,varargin(2:end));
                    end
                % first pass through
                    obj.nodeNames = obj.sourceData.EPID;
                    if isempty(obj.sourceData.EPID)
                        obj.nodeNames = obj.sourceData.FastaName;
                    end
                    obj.time = obj.sourceData.isolateDate;
                    for P=reshape(fieldnames(obj.location),1,[])
                        obj.location.(P{:})=obj.sourceData.(P{:});
                    end    

                    obj.findMLNetwork('clockRate',options.clockRate,'distanceMethod',options.distanceMethod,'alpha',options.alpha,'FDR',options.FDR);

                elseif isa(obj.sourceData,'treeParser.BEAST') || isa(obj.sourceData,'treeParser')
                    options=struct('treeIdx',1,'clockRate',0.01/365.2431,'treeGenomeLength',903,'alpha',0.05,'distanceMethod','Hasegawa','FDR',true);
                    if nargin>1,
                        options=keyValuePairVararginHandler(options,varargin(2:end));
                    end
                % first pass through
                    obj.nodeNames = obj.sourceData.leafNames;
                    obj.time = obj.sourceData.leafDates;
                    for P=reshape(fieldnames(obj.location),1,[])
                        if isfield(obj.sourceData.nodeAttributes{options.treeIdx},P{:})
                            obj.location.(P{:})=obj.sourceData.nodeAttributes{options.treeIdx}.(P{:});
                        end
                    end

                    obj.findMLNetwork('clockRate',options.clockRate,'treeIdx',options.treeIdx,'treeGenomeLength',options.treeGenomeLength,'alpha',options.alpha,'distanceMethod',options.distanceMethod,'FDR',options.FDR);

                else
                    for pair = reshape(varargin,2,[])
                    inpName=pair{1};
                        if any(strcmp(inpName,fieldnames(obj)))
                            obj.(inpName)=pair{2}; 
                        else
                            error(['input name is not a valid class property: ',inpName])
                        end
                    end
%                     error('network:caseDirected:wrongInputClass','Input must be sequenceData or treeParser object.')
                end


                %checkErrors(obj)
            end
        end
    end

    methods %plot
        function [he, hv]=plotMap(obj,traitVisDataObj,varargin)
            [he, hv]=network.plotMap(obj,traitVisDataObj,varargin{:});
        end
        function [he, hv]=plotTimeseries(obj,varargin)
            [he, hv]=network.plotTimeseries(obj,varargin{:});
        end
        function [options, selected] =plotOptionHandler(obj,varargin)
            [options, selected] =network.plotOptionHandler(obj,varargin{:});
        end
    end
    
    methods 
        function obj = randomize(objIn)
            obj=network.caseDirected;
            for F=fieldnames(obj)',
                if ~ismember(F{:},{'sourceData','numNodes','numLinks'})
                    obj.(F{:}) = objIn.(F{:});
                end
            end
            obj.sourceData=[objIn.sourceData,'.randomize'];
            obj.network(:,1)=objIn.network(randperm(length(objIn.network(:,1))),1);
            obj.network(:,2)=objIn.network(randperm(length(objIn.network(:,2))),2);
        end
    end

    
    methods
        function obj = findMLNetwork(obj,varargin)
            options=struct('clockRate',0.01/365.2431,'treeGenomeLength',903,'treeIdx',1,'alpha',0.05,'distanceMethod','Nei-Tamura','FDR',true);
            if nargin>1
                options=keyValuePairVararginHandler(options,varargin);
            end
            
            priorTimes=((2*min(obj.time)-max(obj.time)):1:max(obj.time));
            
            if isa(obj.sourceData,'sequenceData')
                % setup workspace
                    genomeLength=length(obj.sourceData.sequence{1});
                % calculate distances
                    display('calculating distance matrix')
                    v=ver;
                    if any(strcmpi({v.Name},'Bioinformatics Toolbox'))
                        switch options.distanceMethod
                            case {'jukes-cantor','p-distance'}
                                distanceMatrix=seqpdist(obj.sourceData.sequence,'squareform',true,'Alphabet','NT','indels','pairwise-del','method',options.distanceMethod);
                            otherwise
                                seq=regexprep(obj.sourceData.sequence,'[RYSWKMBDHVN]','-');
                                if iscell(seq)
                                    seq=cell2mat(seq);
                                end

                                distanceMatrix=cell(3,1);
                                distanceMatrix{1}=seqpdist(seq(:,1:3:end),'squareform',true,'alphabet','NT','indels','pairwise-del','method',options.distanceMethod); 
                                distanceMatrix{2}=seqpdist(seq(:,2:3:end),'squareform',true,'alphabet','NT','indels','pairwise-del','method',options.distanceMethod); 
                                distanceMatrix{3}=seqpdist(seq(:,3:3:end),'squareform',true,'alphabet','NT','indels','pairwise-del','method',options.distanceMethod); 
                                distanceMatrix=1/3*(distanceMatrix{1}+distanceMatrix{2}+distanceMatrix{3});
                                clear seq;
                        end
                    else
                        display('Bioinformatics not available. Calculating p-distance.')
                        distanceMatrix=zeros(obj.sourceData.numSeq);
                        for ii=1:length(obj.sourceData.sequence)
                            for jj=ii+1:length(obj.sourceData.sequence)
                                tmp1=regexprep(obj.sourceData.sequence{ii},'[nN]','-');
                                tmp2=regexprep(obj.sourceData.sequence{jj},'[nN]','-');
                                distanceMatrix(ii,jj)=sum(tmp1~=tmp2 & tmp1~='-' & tmp2~='-')/sum(tmp1~='-' & tmp2~='-');
                                distanceMatrix(jj,ii)=distanceMatrix(ii,jj);
                            end
                        end
                    end
            elseif isa(obj.sourceData,'treeParser.BEAST') || isa(obj.sourceData,'treeParser')
                genomeLength=options.treeGenomeLength;
                display('calculating distance matrix')
                if strcmp(obj.sourceData.timeUnits{1},'years')
                    distanceMatrix=365.2431*pdist(obj.sourceData.trees{options.treeIdx},'squareform',true)*options.clockRate;
                else
                    distanceMatrix=pdist(obj.sourceData.trees{options.treeIdx},'squareform',true)*options.clockRate;
                end
            end 
                    
            % pre-declare
                ClockProbMatrix=zeros(size(distanceMatrix,1)*(size(distanceMatrix,1)-1)/2,9);  %ii, jj, deltaT1, p(dij|deltaT1), deltaTPrior, p(dij|deltaTPrior), logRatio, p, isDirected

            % Using a fixed clock
            count=0;
            for ii = 1:size(distanceMatrix,1)
                if mod(ii,ceil(size(distanceMatrix,1)/10))==1,
                    display([num2str(ii),' of ',num2str(size(distanceMatrix,1))]);
                end
                for jj = (ii+1):size(distanceMatrix,1)
                    count=count+1;
                    
                    distance = distanceMatrix(ii,jj) * genomeLength; % convert to total number of mutations
                    
                    deltaT1 = obj.time(jj)-obj.time(ii);
                                        
                    expectedDistance=abs(deltaT1)*options.clockRate*genomeLength;

                    dt=max(0,1/(2*options.clockRate*genomeLength)*(distance-options.clockRate*genomeLength*abs(deltaT1)));
                                        
                    cutoff=10;
                    if expectedDistance<=cutoff
                        tmp1 = (expectedDistance.^distance) .* exp(-expectedDistance) / gamma(distance + 1);
                        tmp2 = ((options.clockRate * genomeLength * (abs(deltaT1)+2*dt)).^distance) .* exp(-options.clockRate * genomeLength * (abs(deltaT1)+2*dt)) / gamma(distance + 1);
                    else
                        tmp1 = exp(-(expectedDistance - distance).^2 / (2*distance)) / sqrt(2*pi*distance);
                        tmp2 = exp(-(options.clockRate * genomeLength * (abs(deltaT1)+2*dt) - distance).^2 / (2*distance)) / sqrt(2*pi*distance);
                    end
                    
                    ClockProbMatrix(count,1)=ii;
                    ClockProbMatrix(count,2)=jj;
                    ClockProbMatrix(count,3)=deltaT1;
                    ClockProbMatrix(count,4)=tmp1;
                    
                    ClockProbMatrix(count,5)=dt;
                    ClockProbMatrix(count,6)=tmp2;
                    ClockProbMatrix(count,7)=log(tmp1/tmp2);
                    if isnan(ClockProbMatrix(count,7))
                        ClockProbMatrix(count,7)=-inf;
                    end
                    pvalue=gammainc(max(0,-ClockProbMatrix(count,7)),1/2);
                    ClockProbMatrix(count,8)=pvalue;
                end
            end 
            
            % significance testing
            
                [~, Pidx]=sort(ClockProbMatrix(:,8)); 
                ClockProbMatrix=ClockProbMatrix(Pidx,:); %sorting by p value
                
                if options.FDR
                    Alpha=(1:length(ClockProbMatrix))'*options.alpha/length(ClockProbMatrix);
                else
                    Alpha=options.alpha;
                end
                cutoff=find(ClockProbMatrix(:,8)<Alpha,1,'last');
                if options.FDR  %BKY 2004
                    cutoff=find(ClockProbMatrix(:,8)<Alpha./(1+options.alpha).*(length(ClockProbMatrix))/(length(ClockProbMatrix)-cutoff),1,'last');
                end
                obj.p=ClockProbMatrix(1:cutoff,8);
                obj.alpha=ClockProbMatrix(cutoff+1,8);

                ClockProbMatrix(1:cutoff,9)=1;
                
                obj.network=ClockProbMatrix(1:cutoff,1:2);
                for k=1:length(obj.network)
                    if ClockProbMatrix(k,3)<0
                        obj.network(k,:)=obj.network(k,2:-1:1);
                    end
                end
                obj.clockViolation=ClockProbMatrix(1:cutoff,5);
            
            % re-packaging
                [~,sortIdx]=sort(obj.network(:,1));
                obj.network=obj.network(sortIdx,:);

%                 [nodes,~,nodeLoc]=unique(obj.network);
%                 obj.time=obj.time(nodes);
%                 obj.nodeNames=obj.nodeNames(nodes);
%                 for P=reshape(fieldnames(obj.location),1,[])
%                     if numel(obj.location.(P{:}))>1
%                         obj.location.(P{:})=obj.location.(P{:})(nodes);
%                     end
%                 end
% 
%                 tmp=(1:length(nodes))';
%                 obj.network=reshape(tmp(nodeLoc),[],2);
            
        end
    end
    
    methods
        function obj = pruneTree(obj,varargin)  
            % keep only the parent-child links with the lowest p-value when
            % more than 1 parent per child
            options=struct('includeOnlyNodes',[]);
            inputs=varargin;
            options=keyValuePairVararginHandler(options,inputs);
            childNodes=unique(obj.network(:,2));
            if ~isempty(options.includeOnlyNodes)
                childNodes=unique(obj.network(:,2));
                if islogical(options.includeOnlyNodes)
                    options.includeOnlyNodes=find(options.includeOnlyNodes);
                end
                childNodes=childNodes(ismember(childNodes,(options.includeOnlyNodes)));
            end
            for k=1:length(childNodes)
                idx=find(obj.network(:,2)==childNodes(k));
                if numel(idx)>1
                    parentDates=obj.time(obj.network(idx,1));
                    tmp=(parentDates==max(parentDates));
                    while sum(tmp)>1
                        p=obj.p(idx);
                        if iscell(p);
                            %p=cellfun(@mean,p);
                            P=zeros(length(p),1);
                            for n=1:length(p)
                                P(n)=p{n}(1);
                            end
                            p=P;
                        end
                        tmp=tmp & (p==min(p));
                        tmp(find(tmp,1,'first'))=false;
                    end
                    obj.network(idx(~tmp),:)=nan(sum(~tmp),2);
                end
            end
            idx=~isnan(obj.network(:,1));
            if ~isempty(options.includeOnlyNodes)
                idx=idx & ismember(obj.network(:,1), options.includeOnlyNodes) & ismember(obj.network(:,2), options.includeOnlyNodes);
            end
            obj.network=obj.network(idx,:);
            obj.p=obj.p(idx);
            obj.clockViolation=obj.clockViolation(idx);
        end

        function obj = pruneToMostRecentParentAtEachLocation(obj,varargin)  
            % keep only the parent-child links with the lowest p-value for each parental location when
            % more than 1 parent per child
            options=struct('adminLevel','AdminL2','includeOnlyNodes',[]);
            inputs=varargin;
            options=keyValuePairVararginHandler(options,inputs);
            childNodes=unique(obj.network(:,2));
            if ~isempty(options.includeOnlyNodes)
                childNodes=unique(obj.network(:,2));
                if islogical(options.includeOnlyNodes)
                    options.includeOnlyNodes=find(options.includeOnlyNodes);
                end
                childNodes=childNodes(ismember(childNodes,(options.includeOnlyNodes)));
            end
            for k=1:length(childNodes)
                idx=find(obj.network(:,2)==childNodes(k));
                if numel(idx)>1
                    parentDates=obj.time(obj.network(idx,1));
                    parentLoc=obj.location.(options.adminLevel)(obj.network(idx,1));
                    [~,~,parentLocIdx]=unique(parentLoc);
                    
                    for m=1:max(parentLocIdx)
                        tmpParentLocIdx=parentLocIdx==m;
                        tmp=(parentDates==max(parentDates(tmpParentLocIdx)) & tmpParentLocIdx);

                        while sum(tmp)>1
                            p=obj.p(idx);
                            if iscell(p);
                                %p=cellfun(@mean,p);
                                P=zeros(length(p),1);
                                for n=1:length(p)
                                    P(n)=p{n}(1);
                                end
                                p=P;
                            end
                            tmp=tmp & (p==min(p));
                            tmp(find(tmp,1,'first'))=false;
                        end
                        obj.network(idx(~tmp & tmpParentLocIdx),:)=nan(sum(~tmp & tmpParentLocIdx),2);
                    end
                end
            end
            idx=~isnan(obj.network(:,1));
            if ~isempty(options.includeOnlyNodes)
                idx=idx & ismember(obj.network(:,1), options.includeOnlyNodes) & ismember(obj.network(:,2), options.includeOnlyNodes);
            end
            obj.network=obj.network(idx,:);
            obj.p=obj.p(idx);
            obj.clockViolation=obj.clockViolation(idx);
        end

        function pruneTriangles(obj)
            parentNodes=unique(obj.network(:,1));
            for k=1:length(parentNodes)
                if mod(k,round(length(parentNodes)/100))==0
                    display([num2str(k/length(parentNodes)*100),' percent complete.']);
                end
                for m=1:size(obj.network,1)
                    if all(ismember(obj.network(m,:),obj.network(obj.network(:,1)==parentNodes(k),2)))
                        if obj.time(parentNodes(k)) <= obj.time(obj.network(m,1));
                            obj.network(obj.network(:,1)==parentNodes(k) & obj.network(:,2)==obj.network(m,2),:)=[nan, nan];
                        else
                            obj.network(m,:)=[nan, nan];
                        end
                    end
                end
            end
            idx=~isnan(obj.network(:,1));
            obj.network=obj.network(idx,:);
            obj.p=obj.p(idx);
            obj.clockViolation=obj.clockViolation(idx);
        end
        
        function obj = pruneDuration(obj,durationRange)
            idx=obj.duration<durationRange(1) | obj.duration>durationRange(2);
            obj.network=obj.network(idx,:);
            obj.p=obj.p(idx);
            obj.clockViolation=obj.clockViolation(idx);
        end
    end
    
    methods  %save and load struct
        function saveobj(obj,fileNameString)
            S = saveobj@network.baseClass(obj,fileNameString);
            for F=fieldnames(obj)',
                if ~ismember(F{:},{'sourceData','numNodes','numLinks','duration'})
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