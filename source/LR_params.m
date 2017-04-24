insampfrac = 1;
crossval = 'testset';
numchiefdomSLE = 153;

prunetime = 1;
nonlinfits = 0;
filtertime = 15;

gravity_type = 'input'; % 'input', 'classical' or 'grav2param' or 'grav3param'
levyfit = 'fixalpha'; % 'fixalpha' 'clauset' 'modelfun'

MLE_use = 1; % this uses the MLE estimates for tau2 and rho
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
