classdef baseClass < handle
    % this class stores tree data needed for coalescent likelihood calculations
    % is currently compatable only with simpleCoalecentLogLikelihood
    
    properties (SetAccess = protected)
        nodeNames = {}      % cellarray of names/IDs for each node
        network = []        % [N x 2] set of indices (younger, older) for each link in the network
        p =[];
    end
    
    properties (Dependent = true)
        numNodes
        numLinks
        duration
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
                %checkErrors(obj)
            end
        end
    end
    
    methods % dependent getters
        function NumLinks = get.numLinks(obj)
            NumLinks = size(obj.network,1);
        end
        
        function NumNodes = get.numNodes(obj)
            NumNodes = length(obj.nodeNames);
        end
        
        function Duration = get.duration(obj)
            Duration = obj.time(obj.network(:,2))-obj.time(obj.network(:,1));
        end
    end
    
    methods
        function obj = pruneTriangles(obj)
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
            if isfield(obj.clockViolation)
                obj.clockViolation=obj.clockViolation(idx);
            end
        end
        
        function obj = pruneLatentPeriod(obj,latentPeriod)
            idx=obj.duration>latentPeriod;
            obj.network=obj.network(idx,:);            
            obj.p=obj.p(idx);
        end
    end
    
    methods %plot
        [he, hv] = plotMap(obj,traitVisDataObj,varargin);
        [he, hv] = plotTimeseries(obj,traitVisDataObj,varargin);
        %h = plotDispersion(obj,traitVisDataObj);
        [options, selected] = plotOptionHandler(obj,varargin);
    end
    
    methods
        function [coeff, h] = plotDispersion(obj,traitVisDataObj,varargin)

            options=struct('meanWindowSize',15,'timeCutoff',4*365.244);
            options=keyValuePairVararginHandler(options,varargin);

            duration = obj.time(obj.network(:,2))-obj.time(obj.network(:,1));
            dispersion = nan(size(duration));
            figure
            hist(duration,20);
            mean(duration)/15*40/15,

            for k=1:obj.numLinks

                if strcmp(traitVisDataObj.type,'AdminL2'),
                    traitIdx1 = strcmpi(obj.location.AdminL2{obj.network(k,1)},traitVisDataObj.textAttribute) & ...
                                strcmpi(obj.location.AdminL1{obj.network(k,1)},traitVisDataObj.parentTextAttribute);
                    traitIdx2 = strcmpi(obj.location.AdminL2{obj.network(k,2)},traitVisDataObj.textAttribute) & ...
                                strcmpi(obj.location.AdminL1{obj.network(k,2)},traitVisDataObj.parentTextAttribute);
                elseif strcmp(traitVisDataObj.type,'ClusterID'),
                    traitIdx1 = strcmpi(obj.location.ClusterID{obj.network(k,1)},traitVisDataObj.textAttribute) & ...
                                strcmpi(obj.location.AdminL2{obj.network(k,1)},traitVisDataObj.parentTextAttribute);
                    traitIdx2 = strcmpi(obj.location.ClusterID{obj.network(k,2)},traitVisDataObj.textAttribute) & ...
                                strcmpi(obj.location.AdminL2{obj.network(k,2)},traitVisDataObj.parentTextAttribute);
                elseif strcmp(traitVisDataObj.type,'AdminL1'),
                    traitIdx1 = strcmpi(obj.location.AdminL1{obj.network(k,1)},traitVisDataObj.textAttribute);
                    traitIdx2 = strcmpi(obj.location.AdminL1{obj.network(k,2)},traitVisDataObj.textAttribute);
                end
                
                if any(traitIdx1) && any(traitIdx2)
                    dispersion(k)=(traitVisDataObj.numberAttribute(traitIdx2,1)-traitVisDataObj.numberAttribute(traitIdx1,1)).^2 + ...
                                (traitVisDataObj.numberAttribute(traitIdx2,2)-traitVisDataObj.numberAttribute(traitIdx1,2)).^2;
                end
            end
            timeBins=options.meanWindowSize*(0:1:ceil(max(duration/options.meanWindowSize)));
            meanDisp=zeros(size(timeBins));
            for k=2:length(timeBins)
                idx=duration>timeBins(k-1) & duration <=timeBins(k);
                meanDisp(k)=nanmean(dispersion(idx));
            end

            % fit mean
            dispersionFun = @(a,x) a(1).*(x.^exp(a(2))); 
            %idx=duration<=options.timeCutoff;
            %alpha = nlinfit(duration(idx),distance(idx),dispersionFun,0.5);
            idx=timeBins<=options.timeCutoff;
            coeff = nlinfit(timeBins(idx),meanDisp(idx),dispersionFun,[1,log(1)]);

            h=figure;
            idx=duration<options.timeCutoff;
            plot(duration(idx),dispersion(idx),'rx');
            hold on;
            plot(duration(~idx),dispersion(~idx),'bx');
            plot(timeBins,meanDisp,'m');
            plot(timeBins,dispersionFun(coeff,timeBins),'k');
            coeff(2)=exp(coeff(2));
            xlabel('duration [days]')
            ylabel('dispersion (degree^2)')
            ylim([0 1.1*max(dispersion)]);
            text(options.meanWindowSize, max(dispersion)*0.8,['alpha = ',num2str(coeff(2))]);
            
        end
    end
    
    methods  %save 
        function S = saveobj(obj,fileNameString)
            S.nodeNames = obj.nodeNames;
            S.network = obj.network;
            save(fileNameString,'S');
        end
        
        function obj = reload(obj,fileNameString)
            load(fileNameString,'S');
            obj.nodeNames = S.nodeNames;
            obj.network = S.network;
        end
    end
    
    methods (Static) % load with class type
        function obj = loadobj(fileNameString)
            obj = network.baseClass;
            obj = reload(obj,fileNameString);
        end
    end
end