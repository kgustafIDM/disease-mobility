% written by Kyle B. Gustafson at the Institute for Disease Modeling, Bellevue, WA 
% between June 2016 and April 2017
%%
load ../data/ebolaTree_POTN_SierraLeone.mat

%% for road distances on the Park et al data
load ../data/drive_durdist_fullto_kk152.mat drive_dur_seconds_full drive_dist_meters_full
drive_dist_meters_fullsym = drive_dist_meters_full + drive_dist_meters_full.';
drive_dur_seconds_fullsym = drive_dur_seconds_full + drive_dur_seconds_full.';
%% preliminary, including name hacks for source list of chiefdoms

longlatL3SLE = ebolaAdminL3.numberAttribute;
latlongL3SLE = [];
latlongL3SLE(:,1) = longlatL3SLE(:,2);
latlongL3SLE(:,2) = longlatL3SLE(:,1);

WGS_84_earthradius = 6378137;

nameL3SLE = ebolaAdminL3.textAttribute;
nameL3SLE(20) = {'KoyaE'};
nameL3SLE(80) = {'KoyaN'};
nameL2SLE = ebolaAdminL2.textAttribute;

alllatSLEL3 = ebolaAdminL3.numberAttribute(:,2);
alllonSLEL3 = ebolaAdminL3.numberAttribute(:,1);

alllatSLEL2 = ebolaAdminL2.numberAttribute(:,2);
alllonSLEL2 = ebolaAdminL2.numberAttribute(:,1);

eAdL3inL2 = cell2mat(ebolaAdminL3.parentShapefileID);

% hack away all these name differences
nameL2SLE{strcmp(nameL2SLE,'Port Loko')}='PortLoko';
nameL2SLE{strcmp(nameL2SLE,'Western Urban')}='WesternUrban';
nameL2SLE{strcmp(nameL2SLE,'Western Rural')}='WesternRural';

%% find sorted population of nodes

% THESE NAME CHANGES only need to be done again if you've lost the MAT file
% SierraLeonePop.mat

% import SierraLeonePop.csv as a Table
% change name of second column to chiefdom
% change multiple Koya to KoyaE and KoyaN
% SierraLeonePop(1,:) = [];
% SierraLeonePop.chiefdom{77} = 'KoyaN';
% SierraLeonePop.chiefdom{7} = 'KoyaE';
% SierraLeonePop.chiefdom{74} = 'WaraWaraYagala';
% SierraLeonePop.chiefdom{75} = 'WaraWaraBafodia';
% SierraLeonePop.chiefdom{86} = 'TinkatupaMakamaSafroko';
%
% for k = 1:size(SierraLeonePop,1)
%     tmp=get_tokens(SierraLeonePop.chiefdom{k},' ');
%     if size(tmp,1)>1
%         SierraLeonePop2.chiefdom{k} = strcat(tmp(1),tmp(2));
%     end
% end
% SierraLeonePop = SierraLeonePop2;
% for k = 1:size(SierraLeonePop,1)
%         SierraLeonePop.chiefdom{k} = char(SierraLeonePop.chiefdom{k});
% end

% save SierraLeonePop.mat SierraLeonePop
load('C:\Users\kgustafson\EVD\kyle\SierraLeonePop.mat');
[SLPnames, nid] = sort(SierraLeonePop.chiefdom);
SLPpop = SierraLeonePop.poptotal;
SLPpop = SLPpop(nid);

sortednameL3SLE = sort(nameL3SLE);
size(intersect(sortednameL3SLE,SLPnames))