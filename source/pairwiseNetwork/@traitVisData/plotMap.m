function plotMap(obj,varargin)
% plots Map from traitVisData object
%
% options:
%   centers :: 'on'  : centers (mean of polygons) shown
%   'data' :: vector of data with as many elements as included sites
%   'includeOnly' :: index or names vector with map elements to include
%
options=struct('centers','off','data',[],'includeonly',[],'title','on','colorbar','off','colormap',[],'clim',[] ...
    ,'parentTraitVisDataObj',[],'mapColor',true);
inputs=varargin;
options=keyValuePairVararginHandler(options,inputs);

if isa(obj,'traitVisData')
   
    if isempty(options.includeonly)
        S.polygons=obj.otherAttribute;
        S.names=obj.textAttribute;
        cmap=obj.cmap;
        S.centers=obj.numberAttribute;
    elseif isa(options.includeonly,'double')
        S.polygons=obj.otherAttribute(options.includeonly);
        S.names=obj.textAttribute(options.includeonly);
        cmap=obj.cmap(options.includeonly,:);
        S.centers=obj.numberAttribute(options.includeonly,:);
    elseif isa(options.includeonly,'cell')
        if isa(options.includeonly{1},'char')
            [ix, loc]=ismember(lower(options.includeonly),lower(obj.textAttribute));
            S.polygons=obj.otherAttribute(loc(ix));
            S.names=obj.textAttribute(loc(ix));
            cmap=obj.cmap(loc(ix),:);
            S.centers=obj.numberAttribute(loc(ix),:);
        end
    end
    numTrait=length(S.names);  	      
 
    if ~isempty(options.colormap)
        cmap=options.colormap;
    end
    if ~options.mapColor
        cmap=colormap(nan(size(cmap)));
    end
    
    colormap(cmap);
         
    if ~isempty(obj.otherAttributeType) && strcmp(obj.otherAttributeType,'polygons')
        if numel(options.data)==length(S.names)
            if isempty(options.colormap)
                if isempty(options.clim)
                   options.clim=[min(options.data),max(options.data)];
                end
                if diff(options.clim)==0
                    options.clim=[];
                end
                plotmap(S,options.data,'clim',options.clim)
            else
                if diff(options.clim)==0
                    options.clim=[];
                end
                plotmap(S,options.data,'colormap',cmap,'clim',options.clim)
            end
        else
            plotmap(S,[1:numTrait],'colormap',cmap,'clim',options.clim);
        end
    
        switch options.centers,
            case {'on' , true}
                hold on
                for k=1:length(S.names)
                   plot(S.centers(k,1),S.centers(k,2),'o','MarkerFaceColor','k','MarkerEdgeColor','k','MarkerSize',4)
                end
        end

        if ~isempty(obj.parentOtherAttribute) && strcmpi(obj.parentOtherAttributeType,'polygons')
            [S2.names, loc]=unique(obj.parentTextAttribute);
            S2.polygons=obj.parentOtherAttribute(loc);
            data=nan(1,length(S2.polygons));
            plotmap(S2,data,'linewidth',2);
        end
        axis tight
        axis equal
    else
        if isempty(options.parentTraitVisDataObj)
            switch options.centers,
                case {'on' , true}
                    bbox=plotWorldMap(S.centers(:,1),S.centers(:,2),regexprep(S.names,'_',' '));
                    xlim(bbox(1:2));
                    ylim(bbox(3:4));
                otherwise
                    plotWorldMap;
            end
        else
            options.parentTraitVisDataObj.plotMap;%('data',nan(size(options.parentTraitVisDataObj.textAttribute)));
            switch options.centers,
                case {'on' , true}
                    hold on
                    for k=1:length(S.names)
                       plot(S.centers(k,1),S.centers(k,2),'wo','MarkerFaceColor','w','MarkerEdgeColor','k','MarkerSize',10)
                       text(S.centers(k,1),S.centers(k,2),S.names{k},'color',[0 0.8 0],'fontWeight','bold', 'VerticalAlign','bottom', 'HorizontalAlign','right')            
                    end
                otherwise
%                     for k=1:length(S.names)
%                         text(S.centers(k,1),S.centers(k,2),S.names{k},'color',[0 0.8 0],'fontWeight','bold')            
%                     end
            end            
        end
    end
 
    axis off
    
    switch options.colorbar
        case 'off'
            colorbar off
    end
    switch options.title
        case 'char'
            if strcmp(options.title,'on')
                title(obj.type)
            else
                title(options.title)
            end
    end
 
else
    error('invalide traitVisDataObj')
end

end