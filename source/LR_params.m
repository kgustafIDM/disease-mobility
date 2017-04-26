%% various fixed parameters for the likelihood ratio calculation
% written by Kyle B. Gustafson at the Institute for Disease Modeling, Bellevue, WA 
% between June 2016 and April 2017
%% for cross-validation
% fraction of data in the test set
insampfrac = 1;
% not currently used
crossval = 'testset';

% default 153
numchiefdomSLE = 153;

% prune tree for shortest links
prunetime = 1;
% 
nonlinfits = 0;
filtertime = 15;

gravity_type = 'input'; % 'input', 'classical' or 'grav2param' or 'grav3param'
levyfit = 'fixalpha'; % 'fixalpha' 'clauset' 'modelfun'

% this uses the MLE estimates for tau2 and rho, which must be found
% separately
MLE_use = 0;
% this must be on when running the MLE scan
MLE_scan = 0;

if MLE_use
    MLE_use
elseif MLE_scan
    alevy 
    rho 
    tau2 
    tau1
else
% here, alevy is alpha+1 since q(x) = C|x|^(-(alpha+1))
    alevy = 1.6; %linspace(1,3,21);
    rho   = 1;
    tau2 = 1;
    tau1 = 1;
end

% option3 is used in the publication results
resting = 'option3'; 
chiefpop = 'median';

xcutvalue = 1;
xminvalue = 1;
xmaxvalue = 500;

rad_chiefdom = 10;

timertype = 'window'; % 'cumulative' or 'window'
cumultpts = 1;
twinsize = 50;
dayslide = 50;

% datebounddefault = 0;
% if datebounddefault==1
%     dateboundlow =  735676;%min(ebolaTime);
%     dateboundhigh = 736220;%max(ebolaTime);
% else
%     dateboundlow = 735676+300;
%     dateboundhigh = 735676+550;
% end

xminkm = 1;
% this is roughly the average radius of a chiefdom if they are round and
% equally-sized, which they are not
