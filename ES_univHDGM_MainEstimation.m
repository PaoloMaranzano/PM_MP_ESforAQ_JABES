%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% %%%%%%%%%% Spatio-temporal Event Studies with univariate HDGM %%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%% Part A: HDGM estimation and CV

%%%%% Application: Lockdown in Lombardy and effect on NO2
%%%%% Journal: JABES (METMA X 2022 conference)

clear
clc
close all

%% %%%%% Set working directory
cd('C:\Users\paulm\OneDrive\Documenti\GitHub\PM_MP_ESforAQ_JABES\')
addpath(genpath('C:\Users\paulm\OneDrive\Documenti\GitHub\PM_MP_ESforAQ_JABES'));

%% %%%%% Auxiliary functions
% addpath(genpath('H:/.shortcut-targets-by-id/1hRMzSqAcE5AcmrsOu5k3nw1rmUY1NHyO/QA&COVID19/DSTEM_Code_General/DSTEM_software'));
% addpath(genpath('H:/.shortcut-targets-by-id/1hRMzSqAcE5AcmrsOu5k3nw1rmUY1NHyO/QA&COVID19/DSTEM_Code_General/Matlab_auxfuns'));
% addpath(genpath('H:/.shortcut-targets-by-id/1hRMzSqAcE5AcmrsOu5k3nw1rmUY1NHyO/QA&COVID19/SPASTA2021/Data'));

%%% Output dir
out_path = 'Output_r1_lags12365_typeY_FourieryearY_FourierweekN_BSplTempY/';

%% Loading data
data_building = 1;
if data_building == 1
    run('ES_univHDGM_DataManagement.m');
else
    load('ES_univHDGM_data.mat');
end

%% Model estimation
run('ES_univHDGM_Estimation.m');

%% Results export
run('ES_univHDGM_ExportResults.m');
[DSTEM_ext] = DSTEM_extraction(obj_stem_model);
DSTEM_ext.SpatTime_pars
DSTEM_InfoCrit(obj_stem_model)
writetable(DSTEM_ext.Reg_pars.Beta_tab.NO2_ground,'Beta_HDGM.xlsx');



