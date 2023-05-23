%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% %%%%%%%%%% Spatio-temporal Event Studies with univariate HDGM %%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%% Part B: HDGM estimation and CV

%%%%% Application: Lockdown in Lombardy and effect on NO2
%%%%% Journal: JABES (METMA X 2022 conference)

data_preparation = 1;
standardization = 0;
log_transform = 1;
save_all = 0;
full_estimation = 1;
full_varcov = 1;
CV_step = 0;
save_parts = 1;

parfor_int_1loop_flag = 0;
randK_cv_parfor_int_flag = 0;
stratK_cv_parfor_int_flag = 0;
loso_cv_parfor_int_flag = 0;

time_begin = datetime('now')



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%      Step 0. Data definition     %%
%%    ECMWF + ARPA Lombardia data   %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if data_preparation == 1
    sprintf('Step 0: data definition')
    time_data_definition_begin = datetime('now')
    %% Data building
    poll = Ground.poll;
    n = cell(length(poll),1);
    obj_stem_gridlist_p = stem_gridlist();
    for p = 1:length(poll)
        %%% Dependent variables
        ground.Y{p} = Ground.([poll{p}]);
        ground.Y_name{p} = [poll{p} '_ground'];
        n{p,1} = size(ground.Y{p}, 1);
        %%% Loading covariates for the selected pollutants at each monitoring stations
        ground.X_beta{p} = Ground.(['X_' poll{p}]);
        ground.X_beta_name{p} = Ground.vars_names;
        ground.X_z{p} = ones(n{p}, 1);
        ground.X_z_name{p} = {['constant_' poll{p}]};
        %%% Coordinates grid
        ground.coordinates{p} = [Ground.(['Coords_' poll{p}]).Latitude, ...
            Ground.(['Coords_' poll{p}]).Longitude];
        obj_stem_grid = cell(length(poll),1);
        obj_stem_grid{p} = stem_grid(ground.coordinates{p}, 'deg', 'sparse', 'point');
        %%% stem_gridlist
        obj_stem_gridlist_p.add(obj_stem_grid{p});
    end
    T = Ground.time_stamps;
    %% DSTEM settings
    %%% DSTEM_varset
    obj_stem_varset_p = stem_varset(ground.Y, ground.Y_name, [], [], ...
        ground.X_beta, ground.X_beta_name, ...
        ground.X_z, ground.X_z_name);
    %%% DSTEM_data
    obj_stem_datestamp = stem_datestamp(datetime(datestr(min(Ground.date_time))),...
        datetime(datestr(max(Ground.date_time))),T);
    shape = [];
    obj_stem_modeltype = stem_modeltype('HDGM');
    obj_stem_data = stem_data(obj_stem_varset_p, obj_stem_gridlist_p, ...
        [], [], obj_stem_datestamp, [], ...
        obj_stem_modeltype, shape);
    %%% DSTEM_par
    obj_stem_par_constr_red=stem_par_constraints();
    obj_stem_par_constr_red.time_diagonal=0;
    obj_stem_par = stem_par(obj_stem_data, 'exponential',obj_stem_par_constr_red);
    obj_stem_model = stem_model(obj_stem_data, obj_stem_par);
    %%% Data transform
    if log_transform == 1
        obj_stem_model.stem_data.log_transform;
    end
    if standardization == 1
        obj_stem_model.stem_data.standardize;
    end
    time_data_definition_end = datetime('now')
end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Step 1. Full model estimation %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if full_estimation == 1
    sprintf('Step 1: full model estimation')
    obj_stem_model_full = obj_stem_model;
    %%% Starting values
    obj_stem_par.beta = obj_stem_model_full.get_beta0();
    obj_stem_par.theta_z = km2deg(50);         % Kilometers
    obj_stem_par.v_z = [1];
    % obj_stem_par.v_z = [
    %     1 0.6 0.6;
    %     0.6 1 0.9
    %     0.6 0.9 1];             % Cross-correlation between multiple variables
    obj_stem_par.sigma_eta = diag([0.2]);
    % obj_stem_par.sigma_eta = diag([0.2 0.2 0.2]);
    obj_stem_par.G = diag([0.8]);
    % obj_stem_par.G = diag([0.8 0.8 0.8]);
    obj_stem_par.sigma_eps = diag([0.3]);
    % obj_stem_par.sigma_eps = diag([0.3 0.3 0.3]);
    obj_stem_model_full.set_initial_values(obj_stem_par);
    %%% Parameters estimation
    obj_stem_EM_options = stem_EM_options();
    obj_stem_EM_options.exit_tol_loglike = 0.0001;
    obj_stem_EM_options.exit_tol_par = 0.0001;
    obj_stem_EM_options.max_iterations = 400;
    time_full_estim_begin = datetime('now')
    obj_stem_model_full.EM_estimate(obj_stem_EM_options);
    time_full_estim_end = datetime('now')
    if full_varcov == 1
        %%% Exact VarCov matrix of estimated pars
        time_full_varcov_begin = datetime('now')
        obj_stem_model_full.set_varcov;
        time_full_varcov_end = datetime('now')
        %%% Log-likelihood of the data
        time_full_loglik_begin = datetime('now')
        obj_stem_model_full.set_logL;
        time_full_loglik_end = datetime('now')
        %%% Extract main objects from full model
        [Full_model_MLE,Full_model_pef_met,Full_model_setup] = DSTEM_extraction(...
            obj_stem_model_full,[],[]);
    end
    save([out_path 'Full_model.mat'])
end


if CV_step == 1
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% Step 2. Cross-validation partitioning %%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    sprintf('Step 2: cross-validation partitioning')
    time_CVpart_begin = datetime('now')
    %%% Number of folds
    K = 10;
    %%% Partitioning data
    % Stratified K-fold CV algorithm
    cv_part_stratK = CVpart_strat_Kfold(Ground,'Tipology_rec',K);
    % Random K-fold CV algorithm
    cv_part_randK = CVpart_random_Kfold(Ground,K);
    % Leave-one-station-out CV algorithm
    cv_part_LOSO = CVpart_LeaveOneStatOut(Ground);
    time_CVpart_end = datetime('now')
    if save_parts == 1
        save('CV_part.mat','cv_part_stratK','cv_part_randK','cv_part_LOSO',...
            'time_CVpart_begin','time_CVpart_end')
    end



    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% Step 3. Spatio-temporal CV %%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % Definition ActiveSets (solo full model)
    ActiveSets = Full_model_MLE.Reg_pars.Reg_names_tab;
    ActiveSets.Properties.VariableNames = {'Variable','Response_idx'};
    ActiveSets.ActiveSets1 = ones(length(Ground.vars_names),1);
    nAS = sum(contains(ActiveSets.Properties.VariableNames,'ActiveSets'));

    %%%%% Parfor interno con un loop solo per tuttr le tipologie di CV
    if parfor_int_1loop_flag == 1
        time_CV_begin = datetime('now')
        clear MSE_lambdas RMSE_lambdas MAE_lambdas R2_lambdas
        for ASi = 1:nAS
            sprintf('Estimating Active Set number %d of %d', ASi,nAS)
            %%% Spatio-temporal stratified K-fold CV
            if stratK_cv_parfor_int_flag == 1
                sprintf('Step 4a: spatio-temporal stratified K-fold CV for each active set')
                Perf_CV_metrics = DSTEM_KfoldCV_par(cv_part_stratK,obj_stem_model_full,ActiveSets,ASi);
                stratK_CV_metrics.AIC_ASi(:,ASi) = Perf_CV_metrics.AIC;
                stratK_CV_metrics.LogL_ASi(:,ASi) = Perf_CV_metrics.LogL;
                for p = 1:length(poll)
                    stratK_CV_metrics.MSE_ASi{p}(:,ASi) = Perf_CV_metrics.MSE(p,:)';
                    stratK_CV_metrics.RMSE_ASi{p}(:,ASi) = Perf_CV_metrics.RMSE(p,:)';
                    stratK_CV_metrics.MAE_ASi{p}(:,ASi) = Perf_CV_metrics.MAE(p,:)';
                    stratK_CV_metrics.R2_ASi{p}(:,ASi) = Perf_CV_metrics.R2(p,:)';
                end
            end
            %%% Spatio-temporal random K-fold CV
            sprintf('Step 4b: spatio-temporal random K-fold CV for each active set')
            if randK_cv_parfor_int_flag == 1
                Perf_CV_metrics = DSTEM_KfoldCV_par(cv_part_randK,obj_stem_model_full,ActiveSets,ASi);
                randK_CV_metrics.AIC_ASi(:,ASi) = Perf_CV_metrics.AIC;
                randK_CV_metrics.LogL_ASi(:,ASi) = Perf_CV_metrics.LogL;
                for p = 1:length(poll)
                    randK_CV_metrics.MSE_ASi{p}(:,ASi) = Perf_CV_metrics.MSE(p,:)';
                    randK_CV_metrics.RMSE_ASi{p}(:,ASi) = Perf_CV_metrics.RMSE(p,:)';
                    randK_CV_metrics.MAE_ASi{p}(:,ASi) = Perf_CV_metrics.MAE(p,:)';
                    randK_CV_metrics.R2_ASi{p}(:,ASi) = Perf_CV_metrics.R2(p,:)';
                end
            end
            %%% Leave-one-stat-out CV
            if loso_cv_parfor_int_flag == 1
                sprintf('Step 4c: spatial leave-one-stat-out CV for each active set')
                Perf_CV_metrics = DSTEM_KfoldCV_par(cv_part_LOSO,obj_stem_model_full,ActiveSets,ASi);
                LOSO_CV_metrics.AIC_ASi(:,ASi) = Perf_CV_metrics.AIC;
                LOSO_CV_metrics.LogL_ASi(:,ASi) = Perf_CV_metrics.LogL;
                for p = 1:length(poll)
                    LOSO_CV_metrics.MSE_ASi{p}(:,ASi) = Perf_CV_metrics.MSE(p,:)';
                    LOSO_CV_metrics.RMSE_ASi{p}(:,ASi) = Perf_CV_metrics.RMSE(p,:)';
                    LOSO_CV_metrics.MAE_ASi{p}(:,ASi) = Perf_CV_metrics.MAE(p,:)';
                    LOSO_CV_metrics.R2_ASi{p}(:,ASi) = Perf_CV_metrics.R2(p,:)';
                end
            end
        end
        time_CV_end = datetime('now')
        if save_parts == 1
            save('CV_results.mat','stratK_CV_metrics','randK_CV_metrics','LOSO_CV_metrics',...
                'time_CV_begin','time_CV_end')
        end
    end

    %%%%% Fine
    time_end = datetime('now');

    %%%%% Saving complete output
    if save_all == 1
        save('output_ST_ES_univHDGM_stream.mat')
    end
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%      Output saving      %%
%%%%%%%%%%%%%%%%%%%%%%#%%%%%%
if save_all == 1
    date_begin_save = datetime('now')
    save([out_path 'Output_ST_ES_UnivHDGM_' char(poll) '_std'  num2str(standardization) '_log' num2str(log_transform) '_' date '.mat']);
    date_end = datetime('now')
end
