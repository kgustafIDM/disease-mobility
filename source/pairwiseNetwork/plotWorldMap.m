function bbox = plotWorldMap(long,lat,labels,varargin)
% plot points on mercator map of world
% thanks to user Amro at http://stackoverflow.com/questions/11655532/plot-geo-locations-on-worldmap-with-matlab
%
% INPUT: 
%   long :: longitude vector
%   lat  :: latitude vector
%   labels :: labels (optional)
%
%   empty input just displays background map. 
%


if nargin<3
    options=struct('markersize',6,'cmap',[]);
else
    options=struct('markersize',6,'cmap',colormap(autumn(length(unique(labels)))));
end
options=keyValuePairVararginHandler(options,varargin);



% world map in Mercator projection
% fname = 'http://upload.wikimedia.org/wikipedia/commons/thumb/7/74/Mercator-projection.jpg/773px-Mercator-projection.jpg';
% fname = 'http://upload.wikimedia.org/wikipedia/commons/thumb/7/76/World_V2.0.svg/2000px-World_V2.0.svg.png';
fname = 'World_Blank_Map_(Mercator_projection).png';
img = imread(fname);
[imgH,imgW,~] = size(img);

% plot markers on map
image(img), 
hold on

% Mercator projection
if nargin==2
    [x,y] = mercatorProjection(long, lat, imgW, imgH);
    plot(x,y, 'b.', 'MarkerSize',options.markersize)
elseif nargin==3
    [x,y] = mercatorProjection(long, lat, imgW, imgH);
    [uniqueLabels,~,labelPos]=unique(labels);
    for k=1:length(uniqueLabels)
        plot(x(labelPos==k),y(labelPos==k), '.', 'MarkerSize',options.markersize, 'color',options.cmap(k,:));
    end
    
    if isa(labels,'cell') || isa(labels,'char')
        text(x, y, labels, 'Color','w', 'VerticalAlign','bottom', 'HorizontalAlign','right')
    end
end

if nargin>1
    dx=max(x)-min(x);
    dy=max(y)-min(y);
    bbox=[min(x)-dx/5, max(x)+dx/5,min(y)-dy/5, max(y)+dy/5 ];
end

axis off
hold off

end


