%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% %%%%%%%%%% Spatio-temporal Event Studies with univariate HDGM %%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%% Part A: Data management

%%%%% Application: Lockdown in Lombardy and effect on NO2
%%%%% Journal: JABES (METMA X 2022 conference)

%% %%%%%
seas_weather = 0;
stat_type = 1;
lags_weather = {1,2,365};
Fourier_y = 1;
Fourier_w = 0;
BSpline_covs = {'Temperature','WindU','WindV'};
% BSpline_covs = {};
Ground = Reshape_long_to_HDGM_LombardyData('H:\.shortcut-targets-by-id\1hRMzSqAcE5AcmrsOu5k3nw1rmUY1NHyO\QA&COVID19\SPASTA2021\Data', ...
    'data_long.mat','day',{'NO2'},0,1,lags_weather,seas_weather,0,0,0,0,...
    Fourier_y,Fourier_w,BSpline_covs);

%%% (Estimation) Sub-period selection: 01/03/2017 - 01/03/2020
d_estim_start = datetime(2018,01,01,'TimeZone', 'Z');
d_estim_end = datetime(2020,01,31,'TimeZone', 'Z');
%%% (Forecasting) Sub-period selection: 01/03/2020 - 01/06/2020
d_event_start = datetime(2020,02,01,'TimeZone', 'Z');
d_event_end = datetime(2020,05,31,'TimeZone', 'Z');

%%% Covariates selection
if ~isempty(BSpline_covs)
    % Non-linear covariates
    vars_names = setdiff({'cons','Temperature','Rainfall','Pressure','RelHumid',...
        'WindU','WindV',...
        'WindU_max','WindV_max',...
        'VegetationHigh','VegetationLow','GeopotHeight',...
        'weekend','Holidays'},BSpline_covs);
    vars_names = [vars_names, strcat(BSpline_covs,repmat({'_B'},1,length(BSpline_covs)))];
else
    vars_names = {'cons','Temperature','Rainfall','Pressure','RelHumid',...
        'WindU','WindV',...
        'WindU_max','WindV_max',...
        'VegetationHigh','VegetationLow','GeopotHeight','weekend','Holidays'};
end
if stat_type == 1
    % Seasonal covariates
    vars_names = [vars_names, {'ARPA_MetrArea','ARPA_UrbPlain','ARPA_Mountain'}];
    vars_names = [vars_names, {'Ind','Traf','Rural'}];
end
if seas_weather == 1
    % Seasonal covariates
    vars_names = [vars_names, {'Winter','Summer','Autumn'}];
end
if ~isempty(lags_weather)
    % Lagged covariates
    vars_names = [vars_names, 'Lag'];
end
if Fourier_y == 1
    % Lagged covariates
    vars_names = [vars_names, '_365'];
end
if Fourier_w == 1
    % Lagged covariates
    vars_names = [vars_names, '_7'];
end

vars_names = Ground.X_names(contains(Ground.X_names,vars_names));
covs_idx = find(contains(Ground.X_names,vars_names));
for i = 1:length(vars_names)
    vars_names{i} = ['X_beta_' vars_names{i}];
end
Ground.vars_names = vars_names;
Ground.X_names = Ground.X_names(covs_idx);

%%% Defining time indices
t = datetime(Ground.date_time, 'ConvertFrom', 'datenum',...
    'Format', 'yyyy-MM-dd HH:mm','TimeZone', 'Z');
date_idx = find(t >= d_estim_start & t < d_event_end);
Ground.date_time = Ground.date_time(date_idx);
Ground.time_stamps = length(Ground.date_time);


%%% Filtering original HDGM format (full period)
for p = 1:length(Ground.poll)
    %%% Store original data
    Ground.([Ground.poll{p} '_original']) = Ground.(Ground.poll{p});
    Ground.(['X_' Ground.poll{p} '_original']) = Ground.(['X_' Ground.poll{p}])(:,covs_idx,:);
    %%% Filtering covariates and dates
    Ground.(Ground.poll{p}) = Ground.(Ground.poll{p})(:,date_idx);
    Ground.(['X_' Ground.poll{p}]) = Ground.(['X_' Ground.poll{p}])(:,covs_idx,date_idx);
end


%%% Filtering original HDGM format according to the subperiod
t = datetime(Ground.date_time, 'ConvertFrom', 'datenum',...
    'Format', 'yyyy-MM-dd HH:mm','TimeZone', 'Z');
date_full_idx = find(t >= d_estim_start & t < d_estim_end);
date_estim_idx = find(t >= d_estim_start & t <= d_estim_end);
date_event_idx = find(t >= d_event_start & t <= d_event_end);
for p = 1:length(Ground.poll)
    %%% Subsetting Y
    Ground.([Ground.poll{p} '_full']) = Ground.(Ground.poll{p});
    Ground.([Ground.poll{p} '_event']) = Ground.(Ground.poll{p})(:,date_event_idx);
    Ground.([Ground.poll{p} '_estim']) = Ground.(Ground.poll{p})(:,date_estim_idx);
    %%% Subsetting X
    Ground.(['X_' Ground.poll{p} '_full']) = Ground.(['X_' Ground.poll{p}])(:,date_full_idx);
    Ground.(['X_' Ground.poll{p} '_event']) = Ground.(['X_' Ground.poll{p}])(:,date_event_idx);
    Ground.(['X_' Ground.poll{p} '_estim']) = Ground.(['X_' Ground.poll{p}])(:,date_estim_idx);
    %%% Putting event period as NaN
    Ground.(Ground.poll{p})(:,date_event_idx) = NaN;
end


%%% Remove duplicate tables
% clearvars -except Ground;


save_dataset = 0;
if save_dataset == 1
    save('ES_univHDGM_data.mat','Ground','-v7.3')
end