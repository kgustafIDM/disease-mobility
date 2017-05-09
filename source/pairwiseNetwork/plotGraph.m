function [he, hv] = plotGraph(linkedPair,nodeID,nodeCoord,varargin)
% INPUT
%   linkedPair :: [ numLink x 2 ] array where each pair indicates linked nodes
%   coords     :: [ numNode x 2 ] array with positions of each node
%
% OPTIONS
%
% OUTPUT
%   [he, hv] = handles to edges and vertices
%

options=struct( 'gradientBy','absolute'... %'difference'
                ,'marker','o','markersize',3,'markercolor','b','displayNodes',true ...
                ,'linecolor',[.5 .5 .5],'linewidth',1,'linegradientdata',[],'linecolormap',colormap(autumn(128)) ...
                ,'slice',[0,1] ... % fraction of total line gradient data to slice network at 
                ,'jitter','off');
options=keyValuePairVararginHandler(options,varargin);


% orientation
    if size(linkedPair,1)==2 && size(linkedPair,2)~=2
        linkedPair = linkedPair';
    end
    if size(nodeCoord,1)==2 && size(nodeCoord,2)~=2
        nodeCoord = nodeCoord';
    end
    if ~isempty(options.linegradientdata)
        if size(options.linegradientdata,1)==2 && size(options.linegradientdata,2)~=2
            options.linegradientdata = options.linegradientdata';
        end 
    end

hold on;
colormap(options.linecolormap);
  
% edges
    he=cell(size(linkedPair,1),1);
    [~,loc1]=ismember(linkedPair(:,1),nodeID);
    [~,loc2]=ismember(linkedPair(:,2),nodeID);
    if isempty(options.linegradientdata)
        he = plot([nodeCoord(loc1,1),nodeCoord(loc2,1)]',[nodeCoord(loc1,2),nodeCoord(loc2,2)]','-' ...
                    ,'color',options.linecolor,'LineWidth',options.linewidth);
    elseif size(options.linegradientdata,1)==size(nodeID,1)
        cmax=size(options.linecolormap,1);
        switch options.gradientBy
            case 'difference'
                duration=options.linegradientdata(loc2)-options.linegradientdata(loc1);
                colorSpace=linspace(0,max(duration),cmax);
            case 'absolute'
                colorSpace=linspace(min(options.linegradientdata),max(options.linegradientdata),cmax);
        end
        % graph slicing
        sliceMax=options.slice(2)*cmax;
        sliceMin=options.slice(1)*cmax+1;
        
        for k=1:size(linkedPair,1)
            switch options.gradientBy
                case 'absolute'
                    idx1=find(options.linegradientdata(loc1(k))>=colorSpace,1,'last');
                    idx2=find(options.linegradientdata(loc2(k))<=colorSpace,1,'first');
                    if idx1~=1
                    end
                case 'difference'
                    idx1=1;
                    idx2=find(duration(k)<=colorSpace,1,'first');
                case 'nodeorder'
                    idx1=1;
                    idx2=(ceil(sliceMax));
            end
            dx=(nodeCoord(loc2(k),1)-nodeCoord(loc1(k),1))/abs(idx1-idx2);
            dy=(nodeCoord(loc2(k),2)-nodeCoord(loc1(k),2))/abs(idx1-idx2);
            x=nodeCoord(loc1(k),1):dx:nodeCoord(loc2(k),1);
            y=nodeCoord(loc1(k),2):dy:nodeCoord(loc2(k),2);
            if idx1<idx2
                c=(idx1+1):idx2;
            elseif idx1>idx2
                c=(idx2+1):idx1;
            else 
                c=idx1;
                x=[nodeCoord(loc1(k),1),nodeCoord(loc2(k),1)];
                y=[nodeCoord(loc1(k),2),nodeCoord(loc2(k),2)];
            end
            switch options.jitter
                case 'on'
                    x=x+randn*(x(end)-x(1))/20;
                    y=y+randn*(y(end)-y(1))/20;
            end
            
            he{k}=nan(min(length(c),ceil(sliceMax))-max(1,sliceMin)+1,1);
            for m=max(1,sliceMin):min(length(c),ceil(sliceMax)),
                he{k}(m)=plot(x(m:m+1),y(m:m+1),'-','color',options.linecolormap(c(m),:),'linewidth',options.linewidth);
            end
        end
        colorbar;
    else
        error('size(linegradientdata) must equal size(nodeID)')
    end

% nodes
    hv=[];
    if options.displayNodes
        for k=1:length(nodeCoord(:,1))
            hv(end+1) = plot(nodeCoord(k,1),nodeCoord(k,2),'.','Marker',options.marker,'MarkerSize',options.markersize ...
                    ,'MarkerFaceColor',options.markercolor,'MarkerEdgeColor',options.markercolor);
        end
    end

end