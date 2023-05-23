%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% %%%%%%%%%% Spatio-temporal Event Studies with univariate HDGM %%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%% Part C: extraction of quantities for ES

%%%%% Application: Lockdown in Lombardy and effect on NO2
%%%%% Journal: METMA X 2022

%% Loading data
% load([out_path 'Full_model.mat']);

%% Extraction of data 
Y_obs = Ground.NO2_full;
Y_obs_estim = Ground.NO2_estim;
Y_obs_event = Ground.NO2_event;
Y_hat = cell2mat(obj_stem_model.stem_EM_result.y_hat_back);
Y_hat_estim = Y_hat(:,1:size(Y_obs_estim,2));
Y_hat_event = Y_hat(:,(size(Y_obs_estim,2)+1:end));
if log_transform == 1
    Y_obs_log = log(Y_obs);
    Y_obs_log_estim = Y_obs_log(:,1:size(Y_obs_estim,2));
    Y_obs_log_event = Y_obs_log(:,(size(Y_obs_estim,2)+1:end));
    Y_hat_log = cell2mat(obj_stem_model.stem_EM_result.y_hat);
    Y_hat_log_estim = Y_hat_log(:,1:size(Y_obs_estim,2));
    Y_hat_log_event = Y_hat_log(:,(size(Y_obs_estim,2)+1:end));
end

%%% Datetime
t = datetime(Ground.date_time', 'ConvertFrom', 'datenum',...
    'Format', 'yyyy-MM-dd HH:mm','TimeZone', 'Z');
t = array2table(t);
t.Properties.VariableNames = {'Date'};

%%% Treatment (Estim or Event)
Window = [repelem("Estim",size(Ground.NO2_estim,2),1) ; repelem("Event",size(Ground.NO2_event,2),1)];
Window = array2table(Window);
Window.Properties.VariableNames = {'Window'};
Window.Window = categorical(Window.Window);

%%% Data
for st = 1:Ground.sites
    for cov = 1:length(Ground.vars_names)
        % Covariates
        vx = Ground.X_NO2(st,cov,:);
        vcov(:,cov) = vx(:);
        % Observed Y
        vy = Y_obs(st,:)';
        % Estimated Y (original scale)
        vy_hat = Y_hat(st,:)';
        if log_transform == 1
            % Observed log(Y)
            vy_log = Y_obs_log(st,:)';
            % Estimated Y (log scale)
            vy_hat_log = Y_hat_log(st,:)';
        end
        % Station metadata
        st_code = repelem(Ground.ARPA_stats_reg.New_cod_stz(st),length(vy_hat),1);
        st_type = repelem(Ground.ARPA_stats_reg.Tipology(st),length(vy_hat),1);
        st_type_rec = repelem(Ground.ARPA_stats_reg.Tipology_rec(st),length(vy_hat),1);
        st_zone_rec = repelem(Ground.ARPA_stats_reg.ARPA_zone_rec(st),length(vy_hat),1);
        st_name = repelem(Ground.ARPA_stats_reg.NameStation(st),length(vy_hat),1);
        st_long = repelem(Ground.ARPA_stats_reg.Longitude(st),length(vy_hat),1);
        st_lat = repelem(Ground.ARPA_stats_reg.Latitude(st),length(vy_hat),1);
        st_alt = repelem(Ground.ARPA_stats_reg.Altitude(st),length(vy_hat),1);
    end
    tab_char = array2table([st_code , st_name , st_type, st_type_rec, st_zone_rec]);
    tab_char.Properties.VariableNames = {'Stz_Code','Stz_Name',...
        'Stz_Type','Stz_Type_rec','Stz_ARPA_zone_rec'};
    if log_transform == 1
        tab_num = array2table([st_long, st_lat, st_alt, vy, vy_log, ...
            vy_hat, vy_hat_log, vcov]);
        tab_num.Properties.VariableNames = {'Stz_Long','Stz_Lat','Stz_Alt',...
            Ground.var_names{:},...
            [Ground.var_names{:} '_log'],...
            [Ground.var_names{:} '_hat_HDGM'],...
            [Ground.var_names{:} '_hat_log_HDGM'],...
            Ground.X_names{:}};
    else
        tab_num = array2table([st_long, st_lat, st_alt, vy, vy_hat, vcov]);
        tab_num.Properties.VariableNames = {'Stz_Long','Stz_Lat','Stz_Alt',...
            Ground.var_names{:},...
            [Ground.var_names{:} '_hat_HDGM'],...
            Ground.X_names{:}};
    end
    tab_output = [t, Window, tab_char , tab_num];
    output_cell{st} = tab_output;
end
output_tab = vertcat(output_cell{:});

%%%%% Export csv
writetable(output_tab,[out_path 'HDGM_output.csv'])

