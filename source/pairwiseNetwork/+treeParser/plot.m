function [ht, phyData]=plot(treeParserObj,varargin)
% Generate tree plots with observed cases colored by attribute
%
% function generateAndPlotPhyloClusters(Lineages,LocationData,strainRoot)
%
%   INPUT
%       treeParserObj :: tree data object
%
%       display options ::
%
%           'type', plottype :: Options are: 'square' (default), 'angular', 'radial',
%               'equalangle', and 'equaldaylight'
%           'colorbar', 'off' :: turn off colorbar
%              'showalllocations', 'on' :: by default, colorbar only shows locations 
%                               that are used in the plot. setting this to 'on'
%                               will show all locations on colorbar if
%                               colorbar is on. ignored if colorbar off
%           'titles' , 'off' :: turn off titles
%
%   OUTPUT
%       [ht phyData] :: plot handles and plot data
%

if isempty(treeParserObj.trees)
    ht=[];
    phyData=[];
    return
end

% default options    
    options = struct('treeIdx',1,'traitVisData',false,'leafVisData',false,'nodeVisData',false,'trait',false ...
                ,'leaftrait',false,'nodetrait',false,'type', 'square','colorbar','on','showalllocations','off','titles','on','leaflabels','off' ...
                ,'showTraitsFor',true,'markersize',6);
    inputs = varargin;
    options = keyValuePairVararginHandler(options,inputs);

    % show all traits as a logical vector in order of nodeAttributes
    if all(options.showTraitsFor)
        options.showTraitsFor=true(2*treeParserObj.numLeaves-1,1);
    end
    % is show traits passed as number list, convert to logical
    if ~islogical(options.showTraitsFor)
        tmp=false(2*treeParserObj.numLeaves-1,1);
        options.showTraitsFor=tmp(options.showTraitsFor);
        clear tmp;
    end
  
% generate Layout 
    figure
    phyHF=gcf();
    phyHA=gca();

    [ht, phyData] = plotphylo(treeParserObj.trees{options.treeIdx}, phyHF, phyHA, @phyPlotFcn,  'Type', options.type, 'Orientation', 'left', 'Rotation', 0, 'BranchLabels', false, 'LeafLabels', false, 'TerminalLabels', false, 'LLRotation',true);

% align to case times
    switch options.type
        case 'square'
            % position leaves by explicit isolate date in years
            if strcmp(treeParserObj.timeUnits,'years')
                DT=max(datenum2years(treeParserObj.leafDates))-max(max(phyData.BranchLinesX));
                phyData.BranchLinesX=phyData.BranchLinesX +DT;
                phyData.BranchDotsX=phyData.BranchDotsX +DT;
                phyData.LeafDotsX=phyData.LeafDotsX +DT;
                phyData.BranchNodeLabelX=phyData.BranchNodeLabelX +DT;
                phyData.LeafNodeLabelX=phyData.LeafNodeLabelX +DT;
                phyData.TerminalNodeLabelsX=phyData.TerminalNodeLabelsX+DT;
            elseif strcmp(treeParserObj.timeUnits,'days')
                DT=max(treeParserObj.leafDates)-max(max(phyData.BranchLinesX));
                phyData.BranchLinesX=datenum2years(phyData.BranchLinesX +DT);
                phyData.BranchDotsX=datenum2years(phyData.BranchDotsX +DT);
                phyData.LeafDotsX=datenum2years(phyData.LeafDotsX +DT);
                phyData.BranchNodeLabelX=datenum2years(phyData.BranchNodeLabelX +DT);
                phyData.LeafNodeLabelX=datenum2years(phyData.LeafNodeLabelX +DT);
                phyData.TerminalNodeLabelsX=datenum2years(phyData.TerminalNodeLabelsX+DT);                
            else
                DT=0;
                warning('unknown time units')
            end
    end

% plotting with desired labels and colors
    clf

    % plot branchlines
        hold on
        ht=plot(phyData.BranchLinesX,max(max(phyData.BranchLinesY))-phyData.BranchLinesY,'color',[.5, .5, .5]);

                
        
% post-layout relabeling if given trait data
    leafTraitName={};
    textLeafTraitName={};
    nodeTraitName={};
    textNodeTraitName={};
    if islogical(options.traitVisData) && islogical(options.leafVisData) && islogical(options.nodeVisData) && islogical(options.trait)% no trait given
        plot(phyData.LeafDotsX,max(phyData.LeafDotsY)-phyData.LeafDotsY,'o','MarkerEdgeColor','b','MarkerFaceColor','b','MarkerSize',options.markersize)
        %plot(phyData.BranchDotsX,max(phyData.LeafDotsY)-phyData.BranchDotsY,'s','MarkerEdgeColor','r','MarkerFaceColor','b','MarkerSize',5)
    else
        % display trait that does not carry color and set information
        if (~islogical(options.trait)&& any(ismember(options.trait,fieldnames(treeParserObj.nodeAttributes{options.treeIdx}))))
            % leaves
                for k=1:treeParserObj.numLeaves,
                    [ix, leafLocIdx]=ismember(phyData.LeafNodeLabelName(k),treeParserObj.nodeNames{options.treeIdx});
                    leafLocIdx=leafLocIdx(ix);
                    if options.showTraitsFor(leafLocIdx)
                        switch class(treeParserObj.nodeAttributes{options.treeIdx}.(options.trait))
                            case 'double'
                                switch options.trait
                                    case 'height_95_HPD'
                                    otherwise
                                        if ~all(isnan(treeParserObj.nodeAttributes{options.treeIdx}.(options.trait)(leafLocIdx,:)))
                                            plot([phyData.LeafDotsX(k)+mean(diff(phyData.LeafDotsX)),max(phyData.LeafDotsX)],(max(phyData.LeafDotsY)-phyData.LeafDotsY(k))*ones(1,2),':','color',[0.7,0.7,0.7]);                    
                                            text(max(phyData.LeafDotsX)+mean(diff(phyData.LeafDotsX)),max(phyData.LeafDotsY)-phyData.LeafDotsY(k),regexprep([num2str(treeParserObj.nodeAttributes{options.treeIdx}.(options.trait)(leafLocIdx,:),4)],'_',' '));
                                            textLeafTraitName={options.trait};
                                        end
                                end
                            case 'cell'
                                    plot([phyData.LeafDotsX(k)+mean(diff(phyData.LeafDotsX)),max(phyData.LeafDotsX)],(max(phyData.LeafDotsY)-phyData.LeafDotsY(k))*ones(1,2),':','color',[0.7,0.7,0.7]);                    
                                    text(max(phyData.LeafDotsX)+mean(diff(phyData.LeafDotsX)),max(phyData.LeafDotsY)-phyData.LeafDotsY(k),regexprep([treeParserObj.nodeAttributes{options.treeIdx}.(options.trait){leafLocIdx,:}],'_',' '));
                                    textLeafTraitName={options.trait};
                        end
                    end
                end
            % nodes
                textNodeTraitName={options.trait};
                for k=1:(treeParserObj.numLeaves-1),
                    [ix, nodeLocIdx]=ismember(phyData.BranchDotsLabelName(k),treeParserObj.nodeNames{options.treeIdx});
                    nodeLocIdx=nodeLocIdx(ix);
                    if options.showTraitsFor(nodeLocIdx)
                        switch class(treeParserObj.nodeAttributes{options.treeIdx}.(options.trait))
                            case 'double'
                                switch options.trait
                                    case 'height_95_HPD'
                                        switch options.type
                                            case 'square'
                                                    plot(max(phyData.LeafDotsX)-treeParserObj.nodeAttributes{options.treeIdx}.(options.trait)(nodeLocIdx,:),max(phyData.LeafDotsY)-phyData.BranchDotsY(k)*ones(1,2),'color',[0, 0, 1],'LineWidth',1.5);
        %                                         case 'radial'  % doesn't really make sense! there's no absolute time on a radial tree---the root jumps around
        %                                             radius=sqrt(max(phyData.LeafDotsX.^2)+max(phyData.LeafDotsY).^2);
        %                                             theta=atan((phyData.LeafDotsY(end)-phyData.BranchDotsY(k))/(phyData.BranchDotsX(k)-phyData.BranchDotsX(end)));
        %                                             plot(phyData.BranchDotsX(end)+(radius-treeParserObj.nodeAttributes{options.treeIdx}.(options.trait)(nodeLocIdx,:))*cos(theta),phyData.BranchDotsY(end)+(radius-treeParserObj.nodeAttributes{options.treeIdx}.(options.trait)(nodeLocIdx,:))*sin(theta),'color',[0, 0, 1],'LineWidth',1.5);
                                            otherwise
                                        end
                                    otherwise
                                        text(phyData.BranchDotsX(k),max(phyData.LeafDotsY)-phyData.BranchDotsY(k),regexprep([num2str(treeParserObj.nodeAttributes{options.treeIdx}.(options.trait)(nodeLocIdx,:),4)],'_',' '));
                                end
                            case 'cell'
                                text(phyData.BranchDotsX(k),max(phyData.LeafDotsY)-phyData.BranchDotsY(k),regexprep([treeParserObj.nodeAttributes{options.treeIdx}.(options.trait){nodeLocIdx,:}],'_',' '));
                        end
                    end                
                end
         end
        
        % plot case leaves with trait info from a traitVisData object
        if (~islogical(options.traitVisData) && any(ismember(options.traitVisData.type,fieldnames(treeParserObj.nodeAttributes{options.treeIdx}))))...
                || (~islogical(options.nodeVisData) && any(ismember(options.nodeVisData.type,fieldnames(treeParserObj.nodeAttributes{options.treeIdx}))))...
                || (~islogical(options.leafVisData) && any(ismember(options.leafVisData.type,fieldnames(treeParserObj.nodeAttributes{options.treeIdx}))))
            
            % leaves
                leafTraitName={};
                for k=1:treeParserObj.numLeaves,
                    [ix, leafLocIdx]=ismember(phyData.LeafNodeLabelName(k),treeParserObj.nodeNames{options.treeIdx});
                    leafLocIdx=leafLocIdx(ix);         
                    if options.showTraitsFor(leafLocIdx)
                        if islogical(options.leafVisData) && islogical(options.traitVisData)
                            plot(phyData.LeafDotsX(k),max(phyData.LeafDotsY)-phyData.LeafDotsY(k),'o','MarkerEdgeColor','b','MarkerFaceColor','b','MarkerSize',options.markersize)
                            leafTraitName={};
                        elseif islogical(options.leafVisData) && ~islogical(options.traitVisData)
                            [~, leafColorIdx]=ismember(lower(treeParserObj.nodeAttributes{options.treeIdx}.(options.traitVisData.type)(leafLocIdx)),lower(options.traitVisData.textAttribute));
                            leafTraitName{end+1}=treeParserObj.nodeAttributes{options.treeIdx}.(options.traitVisData.type){leafLocIdx};
                            if ~leafColorIdx,
                                plot(phyData.LeafDotsX(k),max(phyData.LeafDotsY)-phyData.LeafDotsY(k),'o','MarkerEdgeColor','k','MarkerFaceColor','k','MarkerSize',options.markersize)
                            else
                                plot(phyData.LeafDotsX(k),max(phyData.LeafDotsY)-phyData.LeafDotsY(k),'o','MarkerEdgeColor',options.traitVisData.cmap(leafColorIdx,:),'MarkerFaceColor',options.traitVisData.cmap(leafColorIdx,:),'MarkerSize',options.markersize)
                            end
                        else
                            [~, leafColorIdx]=ismember(lower(treeParserObj.nodeAttributes{options.treeIdx}.(options.leafVisData.type)(leafLocIdx)),lower(options.leafVisData.textAttribute));
                            if iscell(treeParserObj.nodeAttributes{options.treeIdx}.(options.leafVisData.type))
                                leafTraitName{end+1}=treeParserObj.nodeAttributes{options.treeIdx}.(options.leafVisData.type){leafLocIdx};
                            else
                                leafTraitName{end+1}=treeParserObj.nodeAttributes{options.treeIdx}.(options.leafVisData.type)(leafLocIdx);
                            end
                            if ~leafColorIdx,
                                plot(phyData.LeafDotsX(k),max(phyData.LeafDotsY)-phyData.LeafDotsY(k),'o','MarkerEdgeColor','k','MarkerFaceColor','k','MarkerSize',options.markersize)
                            else
                                plot(phyData.LeafDotsX(k),max(phyData.LeafDotsY)-phyData.LeafDotsY(k),'o','MarkerEdgeColor',options.leafVisData.cmap(leafColorIdx,:),'MarkerFaceColor',options.leafVisData.cmap(leafColorIdx,:),'MarkerSize',options.markersize)
                            end
                        end
                    end
                end
                
            % internal nodes
                nodeTraitName={};
                for k=1:(treeParserObj.numLeaves-1)
                    [ix, nodeLocIdx]=ismember(phyData.BranchDotsLabelName(k),treeParserObj.nodeNames{options.treeIdx});
                    nodeLocIdx=nodeLocIdx(ix);
                    if options.showTraitsFor(nodeLocIdx)
                        if islogical(options.nodeVisData) && islogical(options.traitVisData)
                            %plot(phyData.BranchDotsX(k),max(phyData.LeafDotsY)-phyData.BranchDotsY(k),'s','MarkerEdgeColor','r','MarkerFaceColor','b','MarkerSize',options.markersize)
                            nodeTraitName={};
                        elseif islogical(options.nodeVisData) && ~islogical(options.traitVisData) && ~isempty(treeParserObj.nodeAttributes{options.treeIdx}.(options.traitVisData.type){nodeLocIdx})
                            [~, nodeColorIdx]=ismember(lower(treeParserObj.nodeAttributes{options.treeIdx}.(options.traitVisData.type)(nodeLocIdx)),lower(options.traitVisData.textAttribute));
                            nodeTraitName{end+1}=treeParserObj.nodeAttributes{options.treeIdx}.(options.traitVisData.type){nodeLocIdx};
                            if ~nodeColorIdx
                                plot(phyData.BranchDotsX(k),max(phyData.LeafDotsY)-phyData.BranchDotsY(k),'o','MarkerEdgeColor','k','MarkerFaceColor','k','MarkerSize',options.markersize)
                            else
                                plot(phyData.BranchDotsX(k),max(phyData.LeafDotsY)-phyData.BranchDotsY(k),'o','MarkerEdgeColor',options.traitVisData.cmap(nodeColorIdx,:),'MarkerFaceColor',options.traitVisData.cmap(nodeColorIdx,:),'MarkerSize',options.markersize)
                            end
                        elseif ~islogical(options.nodeVisData) && ~isempty(treeParserObj.nodeAttributes{options.treeIdx}.(options.traitVisData.type){nodeLocIdx})
                            [~, nodeColorIdx]=ismember(lower(treeParserObj.nodeAttributes{options.treeIdx}.(options.nodeVisData.type)(nodeLocIdx)),lower(options.nodeVisData.textAttribute));
                            nodeTraitName{end+1}=treeParserObj.nodeAttributes{options.treeIdx}.(options.nodeVisData.type){nodeLocIdx};
                            if ~nodeColorIdx
                                plot(phyData.BranchDotsX(k),max(phyData.LeafDotsY)-phyData.BranchDotsY(k),'o','MarkerEdgeColor','k','MarkerFaceColor','k','MarkerSize',options.markersize)
                            else
                                plot(phyData.BranchDotsX(k),max(phyData.LeafDotsY)-phyData.BranchDotsY(k),'o','MarkerEdgeColor',options.nodeVisData.cmap(nodeColorIdx,:),'MarkerFaceColor',options.nodeVisData.cmap(nodeColorIdx,:),'MarkerSize',options.markersize)
                            end
                        end
                    end
                end
            
        else
            warning('trait not found')
            plot(phyData.LeafDotsX,max(phyData.LeafDotsY)-phyData.LeafDotsY,'o','MarkerEdgeColor','b','MarkerFaceColor','b','MarkerSize',options.markersize)
            leafTraitName={};
            nodeTraitName={};
        end    
        
    % display options
        switch options.colorbar
            case 'on'
                switch options.showalllocations
                    case 'on'
                        leafColorIdx=1:size(options.treeVisData.textAttribute);
                end
                if ~islogical(options.traitVisData)
                    [~, leafColorIdx]=ismember(lower(treeParserObj.nodeAttributes{options.treeIdx}.(options.traitVisData.type)),lower(options.traitVisData.textAttribute));
                    leafColorIdx=leafColorIdx(leafColorIdx>0); % BEAST combines names if exactly equiprobable max aposteriori
                    colormap(options.traitVisData.cmap(unique(leafColorIdx),:))
                    h=colorbar;
                    NumLocUsed=length(options.traitVisData.textAttribute(unique(leafColorIdx)));
                    set(h,'YTick',(.5+(1:NumLocUsed)))
                    set(h,'YTickLabel',options.traitVisData.textAttribute(unique(leafColorIdx)));
                end
        end

        switch options.titles
            case 'on'
                if ~islogical(options.traitVisData)
                   title(options.traitVisData.type);
                end
            case 'off'
                switch options.type 
                    case 'square'
                    otherwise
                        set(textHandle,'visible','off')
                end
        end        

    end
    
% axis display properties

    switch options.leaflabels
        case 'on'
            for k=1:treeParserObj.numLeaves
                plot([phyData.LeafDotsX(k)+mean(diff(phyData.LeafDotsX)),max(phyData.LeafDotsX)],(max(phyData.LeafDotsY)-phyData.LeafDotsY(k))*ones(1,2),':','color',[0.7,0.7,0.7]);
                text(max(phyData.LeafDotsX)+mean(diff(phyData.LeafDotsX)),max(phyData.LeafDotsY)-phyData.LeafDotsY(k),phyData.LeafNodeLabelName{k});
            end
    end
    
    axis tight
    box off
    dcm_obj=datacursormode(gcf);

    switch options.type
        case 'square'
            xlabel(['Isolate Date (',treeParserObj.timeUnits{1},')'])
            set(gca,'YTick',[])
            % data browsing on plot
                set(dcm_obj,'UpdateFcn',{@InfoOnTipsSquare,treeParserObj.timeUnits{1},phyData.LeafDotsX,phyData.LeafDotsY,leafTraitName,phyData.LeafNodeLabelName,phyData.BranchDotsX,phyData.BranchDotsY,nodeTraitName,phyData.BranchDotsLabelName,textLeafTraitName,textNodeTraitName})
        otherwise 
            xdiff=round((max(phyData.LeafDotsX)-min(phyData.LeafDotsX))/2);
            if xdiff>4,
                xdiff=4;
            end
            if xdiff==0
                xdiff=(max(phyData.LeafDotsX)-min(phyData.LeafDotsX))/4;
            end
            
            set(gca,'YTick',[],'XTick',[])
            axis off
            axis square
            % data browsing on plot
            switch options.type
                case 'radial'
                    set(dcm_obj,'UpdateFcn',{@InfoOnTipsRadial,treeParserObj.timeUnits{1},phyData.LeafDotsX,phyData.LeafDotsY,leafTraitName,phyData.LeafNodeLabelName,phyData.BranchDotsX,phyData.BranchDotsY,nodeTraitName,phyData.BranchDotsLabelName,textLeafTraitName,textNodeTraitName,treeParserObj.leafDates})
                    plot([min(phyData.LeafDotsX) min(phyData.LeafDotsX)+xdiff],zeros(1,2),'k','LineWidth',2)
                    text(min(phyData.LeafDotsX)+xdiff/2,0,[num2str(xdiff,2),' ',treeParserObj.timeUnits{1}],'HorizontalAlignment','center','VerticalAlignment','top');
                otherwise
                    set(dcm_obj,'UpdateFcn',{@InfoOnTipsOther,treeParserObj.timeUnits{1},phyData.LeafDotsX,phyData.LeafDotsY,leafTraitName,phyData.LeafNodeLabelName,phyData.BranchDotsX,phyData.BranchDotsY,nodeTraitName,phyData.BranchDotsLabelName,textLeafTraitName,textNodeTraitName,treeParserObj.leafDates})
                    plot([min(phyData.LeafDotsX) min(phyData.LeafDotsX)+xdiff],zeros(1,2),'k','LineWidth',2)
                    text(min(phyData.LeafDotsX)+xdiff/2,0,['~ ',num2str(xdiff,2),' ',treeParserObj.timeUnits{1}],'HorizontalAlignment','center','VerticalAlignment','top');
            end
    end
    
end


function output_txt=InfoOnTipsSquare(obj, event_obj,timeUnits,leafX,leafY,leafTraitName,leafNames,nodeX,nodeY,nodeTraitName,nodeNames,textLeafTraitName,textNodeTraitName)
    pos=get(event_obj,'Position');
    x=pos(1); y=pos(2);
    output_txt={datestr(datenumCleaner(x,timeUnits),1)};
    % leaf data
        Idx=(x==leafX);
        if ismember(y,max(leafY)-leafY(Idx)),
            name=leafNames{max(leafY)-y};
            if ~isempty(name)
                output_txt{end+1}=['CaseID: ',name];
            end
            if ~isempty(leafTraitName)
                att=leafTraitName{max(leafY)-y};
                if ~isempty(att)
                    output_txt{end+1}=['color: ',att];
                end
            end
            if ~isempty(textLeafTraitName)
                output_txt{end+1}=['text: ',textLeafTraitName{:}];
            end
        end
    % node data
        Idx=(x==nodeX);
        if ismember(y,max(leafY)-nodeY(Idx)),
            name=nodeNames{Idx};
            if ~isempty(name)
                output_txt{end+1}=name;
            end
            if ~isempty(nodeTraitName)
                att=nodeTraitName{Idx};
                if ~isempty(att)
                    output_txt{end+1}=['color: ',att];
                end
            end
            if ~isempty(textNodeTraitName)
                if strcmp(textNodeTraitName,'height_95_HPD')
                    output_txt{end+1}=['interval: ',textNodeTraitName{:}];
                else
                    output_txt{end+1}=['text: ',textNodeTraitName{:}];
                end
            end
        end
end

function output_txt=InfoOnTipsRadial(obj, event_obj,timeUnits,leafX,leafY,leafTraitName,leafNames,nodeX,nodeY,nodeTraitName,nodeNames,textLeafTraitName,textNodeTraitName,leafDates)
    pos=get(event_obj,'Position');
    x=pos(1); y=pos(2);
    % leaf data
        IdxX=(x==leafX);
        IdxY=(y==(max(leafY)-leafY));
        IdxX2=(x==nodeX);
        IdxY2=(y==(max(leafY)-nodeY));

        output_txt={datestr(max(leafDates)-datenumCleaner(max(sqrt(leafX.^2+leafY.^2))-sqrt(x^2+(y-max(leafY))^2),timeUnits),1)};

        if any(IdxX) && all(IdxX==IdxY),
            name=leafNames{IdxX};
            if ~isempty(name)
                output_txt{end+1}=['CaseID: ',name];
            end
            if ~isempty(leafTraitName)
                att=leafTraitName{IdxX};
                if ~isempty(att)
                    output_txt{end+1}=['color: ',att];
                end
            end
            if ~isempty(textLeafTraitName)
                output_txt{end+1}=['text: ',textLeafTraitName{:}];
            end

        % node data
        elseif any(IdxX2) && all(IdxX2==IdxY2)
            name=nodeNames{IdxX2};
            if ~isempty(name)
                output_txt{end+1}=name;
            end
            if ~isempty(nodeTraitName)
                att=nodeTraitName{IdxX2};
                if ~isempty(att)
                    output_txt{end+1}=['color: ',att];
                end
            end
            if ~isempty(textNodeTraitName)
                if strcmp(textNodeTraitName,'height_95_HPD')
                    output_txt{end+1}=['interval: ',textNodeTraitName{:}];
                else
                    output_txt{end+1}=['text: ',textNodeTraitName{:}];
                end
            end
        end
end

function output_txt=InfoOnTipsOther(obj, event_obj,timeUnits,leafX,leafY,leafTraitName,leafNames,nodeX,nodeY,nodeTraitName,nodeNames,textLeafTraitName,textNodeTraitName,leafDates)
    pos=get(event_obj,'Position');
    x=pos(1); y=pos(2);
    % leaf data
        IdxX=(x==leafX);
        IdxY=(y==(max(leafY)-leafY));
        IdxX2=(x==nodeX);
        IdxY2=(y==(max(leafY)-nodeY));

        % date currently busted!  needs sorted date info as input since
        % dates can't be reconstructed from layout.
        % output_txt={datestr(max(leafDates)-datenumCleaner(max(sqrt(leafX.^2+leafY.^2))-sqrt(x^2+(y-max(leafY))^2),timeUnits),1)};
        output_txt={};
        
        if any(IdxX) && all(IdxX==IdxY),
            name=leafNames{IdxX};
            if ~isempty(name)
                output_txt{end+1}=['CaseID: ',name];
            end
            if ~isempty(leafTraitName)
                att=leafTraitName{IdxX};
                if ~isempty(att)
                    output_txt{end+1}=['color: ',att];
                end
            end
            if ~isempty(textLeafTraitName)
                output_txt{end+1}=['text: ',textLeafTraitName{:}];
            end

        % node data
        elseif any(IdxX2) && all(IdxX2==IdxY2)
            name=nodeNames{IdxX2};
            if ~isempty(name)
                output_txt{end+1}=name;
            end
            if ~isempty(nodeTraitName)
                att=nodeTraitName{IdxX2};
                if ~isempty(att)
                    output_txt{end+1}=['color: ',att];
                end
            end
            if ~isempty(textNodeTraitName)
                if strcmp(textNodeTraitName,'height_95_HPD')
                    output_txt{end+1}=['interval: ',textNodeTraitName{:}];
                else
                    output_txt{end+1}=['text: ',textNodeTraitName{:}];
                end
            end
        end
end
