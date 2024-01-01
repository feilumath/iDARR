%% 


folder = fileparts(which(mfilename));
addpath(genpath(folder));
clear folder

parentpath = pwd;
addpath([parentpath '/eigenspace_analysis/']);
addpath([parentpath '/regu_LSE/']);
addpath([parentpath '/iDARR/']);
addpath([parentpath '/DARTR/']);
addpath([parentpath '/plotsFn/']);
addpath(genpath([parentpath '/PRcodes/']));
addpath([parentpath '/regu_tools_Hansen/']);

%% save data to your local DIR:  
% making DIR to your root DIR
if ispc
    SAVE_DIR = [userpath, '\DataAnalysis\iDarr\output_disc\'];
else
    SAVE_DIR = [getenv('HOME'),'/DataAnalysis/iDarr/output_disc/'];     
end
if ~exist(SAVE_DIR,'dir'), mkdir(SAVE_DIR); end 
fig_dir = [SAVE_DIR,'figures/']; 
if ~exist(fig_dir,'dir'), mkdir(fig_dir); end   
addpath(SAVE_DIR);


