function [aX3cdflessfit, aX3cdfmorefit, likliplcdfless, likliplcdfmore] = fit_gravlevy(EVDpop,EVDroaddist)
% written by Kyle B. Gustafson at the Institute for Disease Modeling, Bellevue, WA 
% between June 2016 and April 2017
% ensure that the fixed parameters are in the local space
LR_params;

%%
% acquire the nonzero distances and populations between those locations
% y should be count of occurrences of unique distances

    if strcmp(chiefpop,'median')
        nzidEVDrdm = find(EVDroaddist.median>0);
        pop1 = EVDpop.median(nzidEVDrdm(:),1);
        pop2 = EVDpop.median(nzidEVDrdm(:),2);
        dist = EVDroaddist.median(nzidEVDrdm);
    elseif strcmp(chiefpop,'max')
        nzidEVDrdm = find(EVDroaddist.max>0);
        pop1 = EVDpop.max(nzidEVDrdm(:),1);
        pop2 = EVDpop.max(nzidEVDrdm(:),2);
        dist = EVDroaddist.max(nzidEVDrdm);
    elseif strcmp(chiefpop,'mean')
        nzidEVDrdm = find(EVDroaddist.mean>0);
        pop1 = EVDpop.mean(nzidEVDrdm(:),1);
        pop2 = EVDpop.mean(nzidEVDrdm(:),2);
        dist = EVDroaddist.mean(nzidEVDrdm);
    elseif strcmp(chiefpop,'min')
        nzidEVDrdm = find(EVDroaddist.min>0);
        pop1 = EVDpop.min(nzidEVDrdm(:),1);
        pop2 = EVDpop.min(nzidEVDrdm(:),2);
        dist = EVDroaddist.min(nzidEVDrdm);
    end
    XEVDrdd = [pop1,pop2,dist];

% simple linear fit to the histogram of distances
linkdist_km = XEVDrdd(:,3)./1000;
%
nbinsEVD = 20;

[countX3,edgesX3,binX3] = histcounts(linkdist_km,nbinsEVD,'Normalization','count');
[probX3,edgesX3,binX3] = histcounts(linkdist_km,nbinsEVD,'Normalization','probability');

%% binning to find generalized gravity law fit


binmean_pop1 = zeros(size(probX3));
binmean_pop2 = zeros(size(probX3));
binmean_dist = zeros(size(probX3));

for jj = 1:size(probX3,2)
    gig = find(binX3==jj);
    if isempty(gig)
        0;
    else
        binmean_pop1(jj) = mean(XEVDrdd(gig,1));
        binmean_pop2(jj) = mean(XEVDrdd(gig,2));
        binmean_dist(jj) = mean(XEVDrdd(gig,3));
    end
end

lesslink = find(linkdist_km<xcutvalue);
morelink = find(linkdist_km>xcutvalue);
midlink = find(linkdist_km>xminvalue & linkdist_km<xmaxvalue);

linkdist_less = linkdist_km(lesslink);
linkdist_more = linkdist_km(morelink);

if size(unique(linkdist_less),1)<2
    aX3cdfless = 0;
    likliplcdfless = 0;
else
    [aX3cdfless,bmX3cdfless,likliplcdfless] = plfit(linkdist_less,'xmin',xminvalue);
end
aX3cdflessfit = aX3cdfless;
% aX3cdfless = 0.5;
if size(unique(linkdist_more),1)<2
    aX3cdfmore = 0;
    likliplcdfmore = 0;
else
    [aX3cdfmore,bmX3cdfmore,likliplcdfmore] = plfit(linkdist_more,'xmin',xminvalue);
%     plplot(linkdist_more,1,aX3cdfmore); drawnow; pause(1);
end
aX3cdfmorefit = aX3cdfmore;

end