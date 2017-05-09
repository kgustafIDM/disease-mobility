function obj = createColormap(obj,center,varargin)
% Creates Colormap where each entry gives color of trait determined by
% center and radial distance
%
%
%   INPUT
%       traitVisDataObj :: 
%       center :: string specifying center value
%       OPTIONAL : 'TJK', 'Nigeria','Pakistan',cmap
%
%   OUTPUT
%       traitVisDataObj.cmap :: (NumTraits x 3) colormap
%

if isa(center,'char')
    [ix, loc]=ismember(lower(center),lower(obj.textAttribute));
   
elseif isa(center,'double') && size(center,2)==3 && size(center,1)==length(obj.textAttribute)
    obj.cmap=center;
    return;
elseif isa(center,'double') && size(center,2)~=3
    
    dim=length(center);
    loc=zeros(dim,1);
    for k=1:dim,
        loc(k)=find(center(k),obj.numberAttribute(:,k));
    end
    if all(loc==loc(1))
        ix=1;
    else
        ix=0;
    end
else
    error('invalid inputs.')
end

if ~ix
   error('center is not in Attribute set')
else
    if ~isempty(obj.numberAttribute),
        dim=size(obj.numberAttribute,2);
        nameMap=0;
    else
        dim=1;
        nameMap=1;
    end
end

if dim>3
    error('can only handle up to 3d Attributes')
end

if nameMap,
    if dim>1
        error('name data must be one-dimensional')
    end
    
    obj.cmap=colormap(lines(length(obj.textAttribute)));
    
else
    if dim==1,
        cmap=colormap(lines(length(obj.numberAttribute)));
        [~, sortIdx] = sort(obj.numberAttribute);
        obj.cmap=cmap(sortIdx,:);
    elseif dim==2,
        centerDistance=sqrt((obj.numberAttribute(:,1)-obj.numberAttribute(loc,1)).^2+(obj.numberAttribute(:,2)-obj.numberAttribute(loc,2)).^2);
        if ~isempty(varargin)
            if ischar(varargin{1})
                switch varargin{1}
                    case 'Nigeria'
                        R=1-1*centerDistance/max(centerDistance);
                        G=1-1*(obj.numberAttribute(:,1)-min(obj.numberAttribute(:,1)))/(max(obj.numberAttribute(:,1)-min(obj.numberAttribute(:,1))));
                        B=1-0.8*(obj.numberAttribute(:,2)-min(obj.numberAttribute(:,2)))/(max(obj.numberAttribute(:,2)-min(obj.numberAttribute(:,2))));
                    case 'TJK'
                        R=.2+.08./(.1+(centerDistance-min(centerDistance))./(max(centerDistance)-min(centerDistance)));
                        G=.2+.08./(.1+(obj.numberAttribute(:,1)-min(obj.numberAttribute(:,1)))/(max(obj.numberAttribute(:,1)-min(obj.numberAttribute(:,1)))));
                        B=.2+.08./(.1+(obj.numberAttribute(:,2)-min(obj.numberAttribute(:,2)))/(max(obj.numberAttribute(:,2)-min(obj.numberAttribute(:,2)))));
                    case 'Pakistan'
                        R=1-1*(centerDistance-min(centerDistance))/(max(centerDistance)-min(centerDistance));
                        G=1-1*(obj.numberAttribute(:,1)-min(obj.numberAttribute(:,1)))/(max(obj.numberAttribute(:,1)-min(obj.numberAttribute(:,1))));
                        B=1-0.8*(obj.numberAttribute(:,2)-min(obj.numberAttribute(:,2)))/(max(obj.numberAttribute(:,2)-min(obj.numberAttribute(:,2))));
                end
            else
                    R=varargin{1}(:,1);
                    G=varargin{1}(:,2);
                    B=varargin{1}(:,3);
            end
        else
            R=1-1*(centerDistance-min(centerDistance))/(max(centerDistance)-min(centerDistance));
            G=1-1*(obj.numberAttribute(:,1)-min(obj.numberAttribute(:,1)))/(max(obj.numberAttribute(:,1)-min(obj.numberAttribute(:,1))));
            B=1-1*(obj.numberAttribute(:,2)-min(obj.numberAttribute(:,2)))/(max(obj.numberAttribute(:,2)-min(obj.numberAttribute(:,2))));
        end
        obj.cmap=[R,G,B];
    elseif dim==3,
        R=1-1*(obj.numberAttribute(:,1)-min(obj.numberAttribute(:,1)))/(max(obj.numberAttribute(:,1)-min(obj.numberAttribute(:,1))));
        G=1-.8*(obj.numberAttribute(:,2)-min(obj.numberAttribute(:,2)))/(max(obj.numberAttribute(:,2)-min(obj.numberAttribute(:,2))));
        B=1-1*(obj.numberAttribute(:,3)-min(obj.numberAttribute(:,3)))/(max(obj.numberAttribute(:,3)-min(obj.numberAttribute(:,3))));
        obj.cmap=[R,G,B];
    end
end

end