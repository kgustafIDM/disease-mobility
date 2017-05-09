function [x,y] = mercatorProjection(lon, lat, width, height)
    borderX=0.995;%85;
    borderY=0.985;
    width=width*borderX;
    height=height*borderY;
    x = -8.35*width/360+ mod((lon+180)*width/360, width) ;
    y = (.99)*height/2 - log(tan((lat+90)*pi/360))*width/(2*pi);
    if y>.99*height/2
        y = (1.01)*height/2 - log(tan((lat+90)*pi/360))*width/borderX/(2*pi);
    end
end