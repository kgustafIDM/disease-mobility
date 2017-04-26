% written by Kyle B. Gustafson at the Institute for Disease Modeling, Bellevue, WA 
% between June 2016 and April 2017
%% Useful tools for a subset of the analysis

%% look at fraction of stationary linkages cumulative in time

[sortuniqorigins,isort,iorig] = unique(originL3.median);

rests = zeros(size(sortuniqorigins));
moves = zeros(size(sortuniqorigins));
for ii = 1:size(sortuniqorigins,1)
    % recover the chiefdom id for the origin
    oid = sortuniqorigins(ii);
    destys = destinationL3.median(find(oid==originL3.median));
    for jj = 1:size(destys)
        if oid == destys(jj)
            rests(ii) = rests(ii) + 1;
        else
            moves(ii) = moves(ii) + 1;
        end
    end
end

[loma moma] = sort(pop_counts);

touched = find(mean(allrestprob(moma,:),2)>0);
figure; plot(pop_counts(moma(touched)),mean(allrestprob(moma(touched),:),2),'r.')
xlabel('population');
ylabel('mean probability of self-infection');
title('chiefdoms of maximum population');

%% fraction of stationary linkages at each time window

figure; hold on;
for jj = 1:13
    plot(windomain, allrestprob(moma(touched(jj)),:),'DisplayName', ebolaAdminL3.textAttribute{moma(touched(jj))})
end


%% color-coded maps of probability for spatial models: gravity and power law

scalecol = 14;

cmap=jet(180);
cmap = flipud(cmap);

for chid=153
clear gfcolor gmyColor

gfcolor = -log(gravprobnorm(chid,:));
gmyColor = zeros(size(gfcolor,2),3);
gmyColor = cmap(floor(gfcolor*scalecol),:);

ebolaAdminL3=traitVisData('type','AdminL3','textAttribute',regexprep(SierraLeone.AdminL3.NAME_3,' ','') ...
    ,'numberAttribute',newcenters ...
    ,'shapefileID',SierraLeone.AdminL3.ID_3,'cmap',gmyColor,'otherAttribute',newpolygons ...
    ,'otherAttributeType','polygons','parentTextAttribute',SierraLeone.AdminL3.NAME_2 ...
    ,'parentNumberAttribute',false,'parentShapefileID',SierraLeone.AdminL3.ID_2 ...
    ,'parentOtherAttribute',false,'parentOtherAttributeType','AdminL2' ...
    ,'parentType','AdminL2');

figure; title({ebolaAdminL3.textAttribute{chid},'gravity model '});
ebolaAdminL3.plotMap;

clear lfcolor lmyColor

lfcolor = -log(levyprobnorm(chid,:));
lmyColor = zeros(size(lfcolor,2),3);
lmyColor = cmap(floor(lfcolor*scalecol),:);

ebolaAdminL3=traitVisData('type','AdminL3','textAttribute',regexprep(SierraLeone.AdminL3.NAME_3,' ','') ...
    ,'numberAttribute',newcenters ...
    ,'shapefileID',SierraLeone.AdminL3.ID_3,'cmap',lmyColor,'otherAttribute',newpolygons ...
    ,'otherAttributeType','polygons','parentTextAttribute',SierraLeone.AdminL3.NAME_2 ...
    ,'parentNumberAttribute',false,'parentShapefileID',SierraLeone.AdminL3.ID_2 ...
    ,'parentOtherAttribute',false,'parentOtherAttributeType','AdminL2' ...
    ,'parentType','AdminL2');

figure; title({ebolaAdminL3.textAttribute{chid},'Levy flights'});
%ebolaAdminL2.plotMap('colormap',nan(size(ebolaAdminL2.cmap)))
ebolaAdminL3.plotMap;

chid 

end
