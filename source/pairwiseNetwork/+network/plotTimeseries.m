function [he, hv] = plotTimeseries(obj,traitVisDataObj,varargin)
% must pass traitVisDataObject

[options, selected] = plotOptionHandler(obj,varargin{:});
rng(options.randSeed);
he=[]; hv=[];

if isfield(selected.location, traitVisDataObj.type) && numel(selected.location.(traitVisDataObj.type))>1

    % time
        coords=zeros(selected.numNodes,2);
        if ~strcmp('years',options.isolateDateUnits)
            coords(:,1)=datenum2years(selected.time);
        else
            coords(:,1)=selected.time;
        end
        
        if islogical(options.uniqueNames)
            [uniqueNames]=regexprep(unique(selected.location.(traitVisDataObj.type)),' ','_');
            [~,sortByX]=sort(traitVisDataObj.numberAttribute(:,1));
            tmp=regexprep(lower(traitVisDataObj.textAttribute(sortByX)),' ','_');

            [ix]=ismember(tmp,lower(uniqueNames));
            uniqueNames=tmp(ix);
            
            [~,sortByTime]=sort(selected.time);
            locations=lower(regexprep(selected.location.(traitVisDataObj.type)(sortByTime),' ','_'));
            [~,loc]=ismember(uniqueNames,locations);
            [~,sortByTime]=sort(loc);
            uniqueNames=uniqueNames(sortByTime);
        else
            uniqueNames=options.uniqueNames;
        end
        
        if ~options.forceAdminLevel && length(uniqueNames)>options.upAdminLevel && ~isempty(traitVisDataObj.parentTextAttribute) && ~all(strcmp('',unique(selected.location.(traitVisDataObj.parentType))))
            uniqueNames=unique(selected.location.(traitVisDataObj.parentType));
%             [~,sortByX]=sort(traitVisDataObj.numberAttribute(:,1));
%             tmp=regexprep(lower(traitVisDataObj.textAttribute(sortByX)),' ','_');
% 
%             [ix]=ismember(tmp,lower(uniqueNames));
%             uniqueNames=tmp(ix);
%             
            [~,sortByTime]=sort(selected.time);
            locations=lower(regexprep(selected.location.(traitVisDataObj.parentType)(sortByTime),' ','_'));
            [~,loc]=ismember(uniqueNames,locations);
            [~,sortByTime]=sort(loc);
            uniqueNames=uniqueNames(sortByTime);
            for k=1:selected.numNodes,
                coords(k,2)=find(strcmp(selected.location.(traitVisDataObj.parentType){k},uniqueNames))+0.05*randn;
            end
        else
            for k=1:selected.numNodes,
                coords(k,2)=find(strcmpi(regexprep(selected.location.(traitVisDataObj.type){k},' ','_'),uniqueNames))+0.05*randn;
            end
        end
            
        if options.newFigure
            figure
        end
        [he, hv] = plotGraph(selected.network, selected.nodes, coords,'linewidth',options.linewidth,'markersize',options.markersize);
        grid on
        set(gca,'YTick',1:length(uniqueNames),'YTickLabel',uniqueNames)
        for k=1:length(hv)
            pos=round(get(hv(k),'ydata'));
            col=traitVisDataObj.cmap(ismember(regexprep(lower(traitVisDataObj.textAttribute),'[ -]','_'),lower(uniqueNames(pos))),:);
            set(hv(k),'markerfaceColor',col,'markerEdgecolor',col);
        end
        ylim([0.8 length(uniqueNames)+0.2]);
        if ceil(max(coords(:,1)))-floor(min(coords(:,1)))>=2
            xlim([floor(min(coords(:,1))), ceil(max(coords(:,1)))]);
            set(gca,'XTick',floor(min(coords(:,1))):1:ceil(max(coords(:,1))));
        else
            xlim([floor(min(coords(:,1)*10)), ceil(max(coords(:,1)*10))]/10);
            set(gca,'XTick',[(floor(10*min(coords(:,1))):1:ceil(10*max(coords(:,1))))]/10);
        end
        
else
    % 1:N on y-axis
    figure

end


end