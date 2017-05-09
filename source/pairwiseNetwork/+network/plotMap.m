function [he, hv] = plotMap(obj,traitVisDataObj,varargin)
% must pass traitVisDataObject

he=[];
hv=[];

[options, selected] = plotOptionHandler(obj,varargin{:});
if options.newFigure
    figure;
end
rng(options.randSeed);


if isfield(selected.location,traitVisDataObj.type) && numel(selected.location.(traitVisDataObj.type))>1
    % map
        coords=nan(selected.numNodes,2);
        for k=1:selected.numNodes,
            idx=find(strcmpi(regexprep(selected.location.(traitVisDataObj.type){k},'[\._ ]',''),regexprep(traitVisDataObj.textAttribute,'[\._ ]','')),1);
            if ~isempty(traitVisDataObj.otherAttribute)
                dcenter(1)=max(traitVisDataObj.otherAttribute(idx).X)-min(traitVisDataObj.otherAttribute(idx).X);
                dcenter(2)=max(traitVisDataObj.otherAttribute(idx).Y)-min(traitVisDataObj.otherAttribute(idx).Y);
                coords(k,:)=traitVisDataObj.numberAttribute(idx,:)+dcenter/20.*randn(1,2);
            else
                dcenter=[2 2];
                coords(k,:)=traitVisDataObj.numberAttribute(idx,:)+dcenter/20.*randn(1,2);
            end
        end
        % world map in Mercator projection
        if isempty(traitVisDataObj.otherAttribute)
            fname = '..\..\pairwiseNetwork\World_Blank_Map_(Mercator_projection).png';

            img = imread(fname);
            [imgH,imgW,~] = size(img);
            [coords(:,1),coords(:,2)]=mercatorProjection(coords(:,1),coords(:,2),imgW,imgH);
        end
        
        if options.map
            if ~options.mapColor
                traitVisDataObj.plotMap('colormap',nan(size(traitVisDataObj.cmap)));
            else
                traitVisDataObj.plotMap;
            end
        end
        hold on
        if ~options.dateGradient && ~options.durationGradient
            [he, hv] = plotGraph(selected.network, selected.nodes, coords,'linewidth',options.linewidth,'displayNodes',options.displayNodes,'linecolormap',options.linecolormap,'markersize',options.markersize,'markercolormap',options.markercolormap);
        elseif options.dateGradient && ~options.durationGradient
            [he, hv] = plotGraph(selected.network, selected.nodes, coords,'linewidth',options.linewidth ...
                        ,'linegradientdata',selected.time,'gradientBy','absolute','slice',options.slice,'displayNodes',options.displayNodes,'linecolormap',options.linecolormap,'markersize',options.markersize);
        elseif ~options.dateGradient && options.durationGradient
            [he, hv] = plotGraph(selected.network, selected.nodes, coords,'linewidth',options.linewidth ...
                        ,'linegradientdata',selected.time,'gradientBy','difference','slice',options.slice,'displayNodes',options.displayNodes,'linecolormap',options.linecolormap,'markersize',options.markersize);
        end
else
    % circle graph
    
end

title(selected.title)

if options.dateGradient && ~options.durationGradient
    hm = colorbar;    
    YT=get(hm,'YTick');
    %set(hm, 'YTick',[YT(1) YT(end)],'YTickLabel',{num2str(floor(10*datenum2years(min(selected.time)))/10),num2str(ceil(10*datenum2years(max(selected.time)))/10)})
    if strcmpi(options.isolateDateUnits,'years')
        set(hm, 'YTick',[YT(1) YT(end)],'YTickLabel',{datestr(datenumCleaner(min(selected.time),'years')),datestr(datenumCleaner(max(selected.time),'years'))})
    else
        set(hm, 'YTick',[YT(1) YT(end)],'YTickLabel',{datestr(min(selected.time)),datestr(max(selected.time))})
    end
elseif ~options.dateGradient && options.durationGradient
    hm = colorbar;
    minDur=1e300;
    maxDur=0;
    for k=1:selected.numLinks,
        idx1=selected.nodes==selected.network(k,1);
        idx2=selected.nodes==selected.network(k,2);
        diff=selected.time(idx2)-selected.time(idx1);
        maxDur=max(maxDur,diff);
        minDur=min(minDur,diff);
    end
    YT=get(hm,'YTick');
    if strcmp(options.isolateDateUnits,'years')
        set(hm, 'YTick',[YT(1) YT(end)],'YTickLabel',{num2str(floor(10*(minDur))/10),num2str(ceil(10*(maxDur))/10)})
    else
        set(hm, 'YTick',[YT(1) YT(end)],'YTickLabel',{num2str(floor(10*datenum2years(minDur))/10),num2str(ceil(10*datenum2years(maxDur))/10)})
    end
end

if isfield(selected.location,traitVisDataObj.type) && numel(selected.location.(traitVisDataObj.type))>1
    dcm_obj=datacursormode(gcf);
    set(dcm_obj,'UpdateFcn',{@InfoOnTips,coords,selected})
end
end

function output_txt=InfoOnTips(obj,event_obj,coords,selected)
    pos=get(event_obj,'Position');
    x=pos(1); y=pos(2);
    idxX=x==coords(:,1);
    idxY=y==coords(:,2);
    idx=idxX&idxY;
    
    if ~isempty(idx)
        output_txt=selected.nodeNames(idx);
        output_txt{end+1}=datestr(selected.time(idx),'yyyy-mmm-dd');
        if isfield(selected.location,'AdminL1')
            output_txt{end+1}=['AdminL1: ',selected.location.AdminL1{idx}];
        end
        if isfield(selected.location,'AdminL2')
            output_txt{end+1}=['AdminL2: ',selected.location.AdminL2{idx}];
        end
        if isfield(selected.location,'ClusterID')
            output_txt{end+1}=['ClusterID: ',selected.location.ClusterID{idx}];
        end
    else
        output_txt={};
    end
end