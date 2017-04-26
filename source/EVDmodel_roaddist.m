%% RUN from here
% written by Kyle B. Gustafson at the Institute for Disease Modeling, Bellevue, WA 
% between June 2016 and April 2017

%%
load('../data/ebolaAdmin.mat', 'ebolaAdminL2')
load('../data/ebolaAdmin.mat', 'ebolaAdminL3')
%%

% yes/no to plotting the likelihood ratio
plotLRon = 1;
% yes/no to loading the data
needtoload = 1;

% choice between 'tree' (pruned out redundant generations) and full 'network'
phylotype = 'tree'
% this is the dataset from https://github.com/ebov/space-time
dataset = 'Rambaut'

cleanEboladata;

date0 = 735676;
% 
firstday = [0,151,301];
lastday = [150,300,550];
% firstday = [0];
% lastday = [550];
% firstday = [0];
% lastday = [550];

% firstday = 51:slidedays:550;
% lastday  = firstday-1+windowdays;

gg0 = 0; gp1 = 0; gp2 = 0; alevy=0;

LR_linkage_plot

