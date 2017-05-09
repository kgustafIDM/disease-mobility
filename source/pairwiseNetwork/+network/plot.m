function [ht] = plot(obj,traitVisDataObj,varargin)
% must pass traitVisDataObject

ht=[];

isDirectedMatrix=sparse(obj.network(:,1),obj.network(:,2),true(obj.numLinks,1));

if numel(obj.location.(traitVisDataObj.type))>1
    % map
        centers=zeros(length(traitVisDataObj.otherAttribute),2);
        for k=1:length(traitVisDataObj.otherAttribute)
            centers(k,1)=mean(traitVisDataObj.otherAttribute(k).X);
            centers(k,2)=mean(traitVisDataObj.otherAttribute(k).Y);
        end
        dcenter=mean(abs(diff(centers,1)));

        coords=zeros(obj.numNodes,2);
        for k=1:obj.numNodes,
            idx=find(strcmpi(obj.location.(traitVisDataObj.type){k},regexprep(traitVisDataObj.textAttribute,'\.','')),1);
            coords(k,:)=centers(idx,:)+dcenter/20.*randn(1,2);
        end

        figure;
        traitVisDataObj.plotMap('colormap',nan(size(traitVisDataObj.cmap)));
        hold on
        wgPlot(isDirectedMatrix,coords);

    % time
        coords=zeros(obj.numNodes,2);
        coords(:,1)=datenum2years(obj.time);

        uniqueNames=unique(obj.location.(traitVisDataObj.type));
        for k=1:obj.numNodes,
            coords(k,2)=find(strcmp(obj.location.(traitVisDataObj.type){k},uniqueNames))+0.05*randn;
        end
            
        figure
        wgPlot(isDirectedMatrix,coords);
        axis on
        set(gca,'YTick',1:length(uniqueNames),'YTickLabel',uniqueNames)
        ylim([0.8 length(uniqueNames)+0.2]);
        if ceil(max(coords(:,1)))-floor(min(coords(:,1)))>=2
            xlim([floor(min(coords(:,1))), ceil(max(coords(:,1)))]);
            set(gca,'XTick',floor(min(coords(:,1))):1:ceil(max(coords(:,1))));
        else
            xlim([floor(min(coords(:,1)*10)), ceil(max(coords(:,1)*10))]/10);
            set(gca,'XTick',[(floor(10*min(coords(:,1))):1:ceil(10*max(coords(:,1))))]/10);
        end
        
else
    % circle graph
    
    % time
        coords=zeros(obj.numNodes,2);
        coords(:,1)=datenum2years(obj.time);
        coords(:,2)=1:obj.numNodes;
        

        figure
        wgPlot(isDirectedMatrix,coords);
        axis on
        set(gca,'YTick',(1:obj.numNodes),'YTickLabel',obj.nodeNames)
        ylim([0.8 obj.numNodes+0.2]);
        if ceil(max(coords(:,1)))-floor(min(coords(:,1)))>=2
            xlim([floor(min(coords(:,1))), ceil(max(coords(:,1)))]);
            set(gca,'XTick',floor(min(coords(:,1))):1:ceil(max(coords(:,1))));
        else
            xlim([floor(min(coords(:,1)*10)), ceil(max(coords(:,1)*10))]/10);
            set(gca,'XTick',[floor(10*min(coords(:,1))):1:ceil(10*max(coords(:,1)))]/10);
        end
end


end