% written by Kyle B. Gustafson at the Institute for Disease Modeling, Bellevue, WA 
% between June 2016 and April 2017
%%
LR_params;

%%

Or_id = EVDlinkagesSLE(:,1);
De_id = EVDlinkagesSLE(:,2);

[uOrigin, uOr_id, uOr_id2] = unique(Or_id);
EVDlinkagesmintime = [];

for jj = 1:size(uOrigin,1)
    uO_dest = find(jj==uOr_id2);
    [lmt_id_row, lmt_id_col] = find(EVDlink_dur(uO_dest)==min(EVDlink_dur(uO_dest)));
    
    for kk = 1:size(lmt_id_row)
        destmintime = uO_dest(lmt_id_row(kk));
        EVDlinkagesmintime = [EVDlinkagesmintime; uOrigin(jj) De_id(destmintime)];
    end
end

[~,mintime_id,~] = intersect(EVDlinkagesSLE,EVDlinkagesmintime,'rows');

EVDpvaluemintime = EVDpvalueSLE(mintime_id,:);
EVDlink_durmintime = EVDlink_dur(mintime_id,:);

if prunetime
    EVDpvalue = EVDpvaluemintime;
    EVDlinkage = EVDlinkagesmintime;
elseif filtertime > 0
    shortdur = find(EVDlink_dur<filtertime);
    EVDpvalue = EVDpvalueSLE(shortdur,:);
    EVDlinkage = EVDlinkagesSLE(shortdur,:);
else
    EVDpvalue = EVDpvalueSLE;
    EVDlinkage = EVDlinkagesSLE;
end

%%

if strcmp(timertype,'window')
    numtpoints = size(firstday,2);
    windowcases = cell(1,numtpoints);
    %     datecutset = 0:twinsize:dateboundhigh-dateboundlow;
elseif strcmp(timertype,'cumulative')
    daterange = dateboundhigh-dateboundlow;
    numtpoints = cumultpts;
    beforecases = cell(1,numtpoints);
    timestep = floor(daterange/numtpoints);
    datecutset = timestep:timestep:timestep*numtpoints;
elseif strcmp(timertype,'MLEscan')
    numtpoints = 1;
    windowcases = cell(1,numtpoints);
else
    error('timer type bad');
end

% non-normalized likelihood ratio
RL1L2test = zeros(1,numtpoints);
Ptest = zeros(1,numtpoints);
casesat = zeros(1,numtpoints);
linkagestest = cell(1,numtpoints);
% normalized likelihood ratio
normRtest = zeros(1,numtpoints);
keepidtest = cell(1,numtpoints);
alldrvtest = cell(1,numtpoints);
llGravtest = cell(1,numtpoints);
llLevytest = cell(1,numtpoints);

originL3_maxtest = cell(1,numtpoints);
originL3_mintest = cell(1,numtpoints);
originL3_meantest = cell(1,numtpoints);
originL3_mediantest = cell(1,numtpoints);
destinationL3_maxtest = cell(1,numtpoints);
destinationL3_mintest = cell(1,numtpoints);
destinationL3_meantest = cell(1,numtpoints);
destinationL3_mediantest = cell(1,numtpoints);

RL1L2val = zeros(1,numtpoints);
normRval = zeros(1,numtpoints);
Pval = zeros(1,numtpoints);
keepidval = cell(1,numtpoints);
alldrvval = cell(1,numtpoints);
llGravval = cell(1,numtpoints);
llLevyval = cell(1,numtpoints);
linkagesval = cell(1,numtpoints);

originL3_maxval = cell(1,numtpoints);
originL3_minval = cell(1,numtpoints);
originL3_meanval = cell(1,numtpoints);
originL3_medianval = cell(1,numtpoints);
destinationL3_maxval = cell(1,numtpoints);
destinationL3_minval = cell(1,numtpoints);
destinationL3_meanval = cell(1,numtpoints);
destinationL3_medianval = cell(1,numtpoints);

alphaCDF_lesstest = zeros(1,numtpoints);
alphaCDF_moretest = zeros(1,numtpoints);
likelipl_lesstest = zeros(1,numtpoints);
likelipl_moretest = zeros(1,numtpoints);

alphaCDF_lessval = zeros(1,numtpoints);
alphaCDF_moreval = zeros(1,numtpoints);
likelipl_lessval = zeros(1,numtpoints);
likelipl_moreval = zeros(1,numtpoints);

aX3cdflf=0;aX3cdfmf=0;likliplcdfl=0;likliplcdfm=0;

% setting the p-value for significant links in the POTN
siglinksnetid = find(EVDpvalue<0.05);

clear datelow datehigh windomain

testrestprob = [];
testmoveprob = [];

clear windowmain
windomain = (lastday+firstday)./2;
%%
for k = 1:numtpoints
    
    if strcmp(timertype,'MLEscan')
        datelow = dateboundlow;
        datehigh = dateboundhigh;
    else
        datelow = date0+firstday(k);
        datehigh = date0+lastday(k);
    end
    
    if strcmp(timertype,'cumulative')
        datebound = dateboundlow + (k)*timestep;
    end
    
    if strcmp(timertype,'window') | strcmp(timertype,'MLEscan')
        windowcases{k} = find(ebolaTime<=datehigh & ebolaTime>=datelow);
        casesat(k) = size(windowcases{k},1);
        postcasecheck = ismember(EVDlinkage(:,2),windowcases{k});
    elseif strcmp(timertype,'cumulative')
        beforecases{k} = find(ebolaTime<datebound); % option for cumulative counting
        casesat(k) = size(beforecases{k},1);
        postcasecheck = ismember(EVDlinkage(:,2),beforecases{k});
    else
        error('bad timertype');
    end
    
    postcaseid = find(postcasecheck>0);
    keepid = intersect(postcaseid,siglinksnetid);
    
    % take a random sample, sized by insampfrac, to create a test set
    testsetid = randperm(size(keepid,1),floor(size(keepid,1)*insampfrac));
    % samples not in the test set are in the validation set
    valsetid = setdiff(1:size(keepid,1),testsetid);
    % linkages in the test set
    EVDlinkagesSLEtest = EVDlinkage(keepid(sort(testsetid)),:);
    linkagestest{k} = EVDlinkagesSLEtest;
    % linkages in the validation set
    EVDlinkagesSLEval = EVDlinkage(keepid(valsetid),:);
    linkagesval{k} = EVDlinkagesSLEval;
    
    if MLE_use
        rho = rhoMLE(k);
        tau1 = 1;
        tau2 = tau2MLE(k);
        alevy=rho;
    end
    %
    % ran gmap_dist.m and saved: save drive_durdist_fullto_kk152.mat drive_dur_seconds_full drive_dist_meters_full
    % load C:\Users\kgustafson\Cell_ParkBedford2015Supp2\seqTrack\drive_durdist_fullto_kk152.mat drive_dur_seconds_full drive_dist_meters_full
    
    % this function gets the origins and destinations for input genetic linkages
    % SEPARATE call later using validation set of data
    [originL3,destinationL3,EVDp,EVDrd] = ...
        get_orig_dest(EVDlinkagesSLEtest,drive_dist_meters_fullsym,drive_dur_seconds_fullsym,...
        nodeL3id,SLPpop,nameL3SLE,sortednameL3SLE);
    
    % if there are cases in the window or cumulatively
    if casesat(k)>0
        % this function compute the fits for gravity model generalization and from Clauset
        [aX3cdflf, aX3cdfmf, likliplcdfl, likliplcdfm] = ...
            fit_gravlevy(EVDp,EVDrd);
        alphaCDF_lesstest(k) = aX3cdflf;
        alphaCDF_moretest(k) = aX3cdfmf;
        likelipl_lesstest(k) = likliplcdfl;
        likelipl_moretest(k) = likliplcdfm;
    end
    
    % this function computes the likelihood ratio
    [allp2,alldrv,llGrav,llLevy,RL1L2,normRgravlevy,pgravlevy,restprob,moveprob,gravgamma,gravprobnorm,...
            levyprobnorm] = ...
        LP_gravlevy(rho,tau1,tau2,alevy,...
        drive_dist_meters_fullsym,casesat(k),originL3,destinationL3,SLPpop,nameL3SLE,...
        sortednameL3SLE,aX3cdfmf,aX3cdflf);
    
    testrestprob = [testrestprob restprob];
    testmoveprob = [testmoveprob moveprob];
    
    RL1L2test(k) = RL1L2;
    normRtest(k) = normRgravlevy;
    Ptest(k) = -log10(pgravlevy);
    
    keepidtest{k} = testsetid;
    alldrvtest{k} = alldrv;
    llGravtest{k} = llGrav;
    llLevytest{k} = llLevy;
    
    originL3_maxtest{k} = originL3.max;
    originL3_mintest{k} = originL3.min;
    originL3_meantest{k} = originL3.mean;
    originL3_mediantest{k} = originL3.median;
    destinationL3_maxtest{k} = destinationL3.max;
    destinationL3_mintest{k} = destinationL3.min;
    destinationL3_meantest{k} = destinationL3.mean;
    destinationL3_mediantest{k} = destinationL3.median;
    
    % this function gets the origins and destinations for input genetic linkages
    % EARLIER evaluation for test dataset
    if size(valsetid,2)>0
        [originL3,destinationL3,EVDp,EVDrd] = ...
            get_orig_dest(EVDlinkagesSLEval,drive_dist_meters_fullsym,drive_dur_seconds_fullsym,...
            nodeL3id,SLPpop,nameL3SLE,sortednameL3SLE);
        
        % if there are cases in the window or cumulatively
        if casesat(k)>0
            % this function compute the fits for gravity model generalization and from Clauset
            [aX3cdflf, aX3cdfmf, likliplcdfl, likliplcdfm] = ...
                fit_gravlevy(EVDp,EVDrd);
            alphaCDF_lessval(k) = aX3cdflf;
            alphaCDF_moreval(k) = aX3cdfmf;
            likelipl_lessval(k) = likliplcdfl;
            likelipl_moreval(k) = likliplcdfm;
        end
        
        % this function computes the likelihood ratio
        [allp2,alldrv,llGrav,llLevy,RL1L2,normRgravlevy,pgravlevy,restprob,moveprob,gravgamma,gravprobnorm,...
            levyprobnorm] = ...
            LP_gravlevy(rho,tau1,tau2,alevy,...
            drive_dist_meters_fullsym,casesat(k),originL3,destinationL3,SLPpop,nameL3SLE,...
            sortednameL3SLE,aX3cdfmf,aX3cdflf);
      
        RL1L2val(k) = RL1L2;
        normRval(k) = normRgravlevy;
        Pval(k) = -log10(pgravlevy);
        keepidval{k} = valsetid;
        alldrvval{k} = alldrv;
        llGravval{k} = llGrav;
        llLevyval{k} = llLevy;

        originL3_maxval{k} = originL3.max;
        originL3_minval{k} = originL3.min;
        originL3_meanval{k} = originL3.mean;
        originL3_medianval{k} = originL3.median;
        destinationL3_maxval{k} = destinationL3.max;
        destinationL3_minval{k} = destinationL3.min;
        destinationL3_meanval{k} = destinationL3.mean;
        destinationL3_medianval{k} = destinationL3.median;
    end
    
end
%% plotting likelihood ratio and p-value

if plotLRon
    startplotid = find(isnan(normRtest(1:end-1)));
    normRtest(startplotid) = 0;
    RL1L2test(startplotid) = 0;
    Ptest(startplotid) = 0;
    normRval(startplotid) = 0;
    RL1L2val(startplotid) = 0;
    Pval(startplotid) = 0;
    
    figure; subplot(2,1,1); yyaxis left
    
    if strcmp(timertype,'window')
        plot(windomain,normRtest,'DisplayName',['rho = ' num2str(rho) '\alpha= ' ...
            num2str(alevy-1)]); ylabel(['normalized LR ']); grid on
        if insampfrac<1
            hold on;
            plot(windomain,normRval,'DisplayName','validate');
        end
    elseif strcmp(timertype,'cumulative')
        plot(datecutset,normRtest,'DisplayName',['rho = ' num2str(rho) '\alpha= ' ...
            num2str(alevy-1)]); ylabel(['normalized LR ']); grid on
        if insampfrac<1
            hold on;
            plot(datecutset,normRval,'DisplayName','validate');
        end
    end
    legend toggle;
    axis([0,550,-10,10]);
    yyaxis right
    if strcmp(timertype,'window')
        plot(windomain,cellfun('length',keepidtest),'DisplayName',['rho = ' num2str(rho)...
            '\alpha= ' num2str(alevy-1)]); ylabel('linkage count'); grid on
        if insampfrac<1
            hold on;
            plot(windomain,cellfun('length',keepidval),'DisplayName','validate');
        end
    elseif strcmp(timertype,'cumulative')
        plot(datecutset,cellfun('length',keepidtest),'DisplayName',['rho = ' num2str(rho)...
            '\alpha= ' num2str(alevy-1)]); ylabel('linkage count'); grid on
        if insampfrac<1
            hold on;
            plot(datecutset,cellfun('length',keepidval),'DisplayName','validate');
        end
    end
    title([phylotype, ' of linkages N=',num2str(size(testsetid,1)),' ',chiefpop,' ',levyfit]);
    if strcmp(timertype,'window')
        xlabel(['center of ',num2str(twinsize),'-day window (days after 18-Mar-2014)']);
    elseif strcmp(timertype,'cumulative')
        xlabel('days after 18-Mar-2014');
    else
        error('timertype wrong');
    end
    
    subplot(2,1,2);
    if strcmp(timertype,'window')
        semilogy(windomain,10.^(-1.*Ptest),'DisplayName',['rho = ' num2str(rho) '\alpha= '...
            num2str(alevy-1)]); ylabel('p value'); grid on
        if insampfrac<1
            hold on;
            semilogy(windomain,10.^(-1.*Pval),'DisplayName','validate');
        end
    elseif strcmp(timertype,'cumulative')
        semilogy(datecutset,10.^(-1.*Ptest),'DisplayName',['rho = ' num2str(rho) '\alpha= '...
            num2str(alevy-1)]); ylabel('p value'); grid on
        if insampfrac<1
            hold on;
            semilogy(datecutset,10.^(-1.*Pval),'DisplayName','validate');
        end
    end
    legend toggle;
    axis([0,550,10^-40,10^0])
end