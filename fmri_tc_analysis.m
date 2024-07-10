%% NHP project2 %%
% L Priestley
% Adapted 16/08/22

clear all; fclose all; clc;

stat = 1;
GLM = 6;
feat_name = 'GLM0102';

timecourse_dir = '/Volumes/Lacie/nature_neuro_submission/data/fmri/tc_data/';
addpath /Volumes/Lacie/nature_neuro_submission/fmt; 
tc_graph_dir = '/Volumes/Lacie/nature_neuro_submission/data/fmri/tc_output';

behav = readtable('/Volumes/Lacie/nature_neuro_submission/data/fmri/fmri_behaviour_with_glmhmm.csv');
behav = table2cell(behav);

writedata = 1; % write out timecourse regression coefficients
%%

animals = {'Ultra', 'Ulrich', 'Vampire', 'Winky'};

sessions={{'MI00543P' 'MI00546P' 'MI00548P' 'MI00555P' 'MI00559P' 'MI00562P' 'MI00566P' 'MI00570P' 'MI00574P' 'MI00577P' 'MI00581P' 'MI00585P' 'MI00591P' 'MI00593P' 'MI00595P' 'MI00597P'};...
    {'MI00565P' 'MI00569P' 'MI00573P' 'MI00576P' 'MI00580P' 'MI00584P' 'MI00592P' 'MI00594P' 'MI00596P' 'MI00598P' 'MI00600P' 'MI00601P' 'MI00602P' 'MI00604P' 'MI00605P' 'MI00606P'};...
    {'MI00615P' 'MI00616P' 'MI00617P' 'MI00618P' 'MI00619P' 'MI00620P' 'MI00621P' 'MI00622P' 'MI00623P' 'MI00624P' 'MI00625P' 'MI00627P'};...
    {'MI00485P' 'MI00486P' 'MI00487P' 'MI00488P' 'MI00489P' 'MI00491P' 'MI00493P' 'MI00496P' 'MI00505P' 'MI00510P' 'MI00513P' 'MI00521P' 'MI00527P' 'MI00531P' 'MI00539P'}};

%roi = {'HB_LP_decision' 'DRN_decision' 'VTA_decision', 'NB_decision', 'SN_decision', 'LC_decision', 'AI_func_decision', 'SMA_decision', 'ACC_func_decision'};
% roi = {'HB_LP_decision' 'DRN_decision' 'VTA_decision', 'NB_decision', 'SN_decision', 'LC_decision', 'MRN_decision', 'vent_decision'};

% roi = {}; % for GLM
% seed_roi = {}; roi = {}; % for PPI

%%

for r = 1:numel(roi)

    z = 0;

    for a = [1,2,3,4]

        sess = 0;

        for s = 1:numel(sessions{a})

            session = sessions{a}{s};
            animal = animals{a};

            z = z+1;

            % prepare behaviour
            REG.totalTime = cell2mat(behav(strcmp(behav(:, 1), animal) & cell2mat(behav(:,2))==s, 3));
            REG.response = cell2mat(behav(strcmp(behav(:, 1), animal) & cell2mat(behav(:,2))==s, 6));
            REG.outcome = cell2mat(behav(strcmp(behav(:, 1), animal) & cell2mat(behav(:, 2))==s, 8));

            REG.rewMag = cell2mat(behav(strcmp(behav(:, 1), animal) & cell2mat(behav(:,2))==s, 4));
            REG.rewProb = cell2mat(behav(strcmp(behav(:, 1), animal) & cell2mat(behav(:, 2))==s, 5));
            REG.rewAvg = cell2mat(behav(strcmp(behav(:, 1), animal) & cell2mat(behav(:, 2))==s, 23));
            REG.rewSd = cell2mat(behav(strcmp(behav(:, 1), animal) & cell2mat(behav(:, 2))==s, 25));
            REG.richness = cell2mat(behav(strcmp(behav(:, 1), animal) & cell2mat(behav(:, 2))==s, 15));
            REG.stochasticity = cell2mat(behav(strcmp(behav(:, 1), animal) & cell2mat(behav(:, 2))==s, 17));
            REG.ev = cell2mat(behav(strcmp(behav(:, 1), animal) & cell2mat(behav(:, 2))==s, 18));
            REG.preResp = cell2mat(behav(strcmp(behav(:, 1), animal) & cell2mat(behav(:, 2))==s, 20));
            REG.preRew = cell2mat(behav(strcmp(behav(:, 1), animal) & cell2mat(behav(:, 2))==s, 22));
            REG.evAvg = cell2mat(behav(strcmp(behav(:, 1), animal) & cell2mat(behav(:, 2))==s, 24));

            REG.pupil_iti = cell2mat(behav(strcmp(behav(:, 1), animal) & cell2mat(behav(:, 2))==s, 12));
            REG.pupil_decision = cell2mat(behav(strcmp(behav(:, 1), animal) & cell2mat(behav(:, 2))==s, 11));
            REG.RT = log(cell2mat(behav(strcmp(behav(:, 1), animal) & cell2mat(behav(:, 2))==s, 7)));
            REG.viterbi_state = cell2mat(behav(strcmp(behav(:, 1), animal) & cell2mat(behav(:, 2))==s, 37));
            REG.state_difference = cell2mat(behav(strcmp(behav(:, 1), animal) & cell2mat(behav(:, 2))==s, 38));
            REG.state_transition = cell2mat(behav(strcmp(behav(:, 1), animal) & cell2mat(behav(:, 2))==s, 39));
            REG.state_increase = cell2mat(behav(strcmp(behav(:, 1), animal) & cell2mat(behav(:, 2))==s, 40));
            REG.state_decrease = cell2mat(behav(strcmp(behav(:, 1), animal) & cell2mat(behav(:, 2))==s, 41));

            REG.constant = ones(length(REG.response),1);

            if GLM == 1  %% GLM3.3

                % load time-series data
                load([timecourse_dir,animal,'/',session,'/',feat_name,'/',([roi{r},'_epoched'])]);

                % regressor design matrix (don't forget constant!)
                dmat = [REG.response, REG.totalTime, REG.constant];

                % contrast design matrix
                contrasts = diag(ones(size(dmat,2),1));

                % remove trials with no response
                dmat(isnan(REG.response),:)=[];
                trial_data(isnan(REG.response),:)=[];

                if    (a==1 && s==8) || (a==3 && s==7) %offerOnset
                    dmat(end-1:end,:)=[];
                end

                % normalise data
                dmat(:,1:(size(dmat,2)-1)) = zscore(dmat(:,1:(size(dmat,2)-1)));

                % beta X time output
                betaOut = ols(trial_data,dmat,contrasts);
                allBeta.(sprintf('%s',feat_name)).(sprintf('%s',roi{r}))(:,:,z) = betaOut;


            elseif GLM == 2  %% GLM3.1

                % load time-series data
                load([timecourse_dir,animal,'/',session,'/',feat_name,'/',([roi{r},'_epoched'])]);

                % regressor design matrix (don't forget constant!)
                dmat = [REG.rewMag, REG.rewProb, REG.rewAvg, REG.rewSd, REG.preResp, REG.rewAvg.*REG.preResp, REG.rewSd.*REG.preResp, REG.pupil_iti, REG.totalTime, REG.constant];
                %dmat = [REG.rewMag, REG.rewProb, REG.rewAvg, REG.rewSd, REG.preResp, REG.pupil_iti, REG.totalTime, REG.constant];
                %dmat = [REG.richness, REG.totalTime, REG.constant];

                % contrast design matrix
                contrasts = diag(ones(size(dmat,2),1));

                % remove trials with missing data in the design matrix
                dmat(isnan(REG.pupil_iti)|isnan(REG.rewSd)|isnan(REG.preResp),:)=[];
                trial_data(isnan(REG.pupil_iti)|isnan(REG.rewSd)|isnan(REG.preResp),:)=[];

                if    (a==1 && s==8) || (a==3 && s==7) %offerOnset
                    dmat(end-1:end,:)=[];
                end

                % normalise data
                dmat(:,1:(size(dmat,2)-1)) = zscore(dmat(:,1:(size(dmat,2)-1)));

                % beta X time output
                betaOut = ols(trial_data,dmat,contrasts);
                allBeta.(sprintf('%s',feat_name)).(sprintf('%s',roi{r}))(:,:,z) = betaOut;

            elseif GLM == 3  %% GLM3.2a (behavioural_history == pursuit)

                % load time-series data
                load([timecourse_dir,animal,'/',session,'/',feat_name,'/',([roi{r},'_epoched'])]);

                % regressor design matrix (don't forget constant!)
                dmat = [REG.rewMag, REG.rewProb, REG.rewAvg, REG.rewSd, REG.pupil_iti, REG.totalTime, REG.constant];

                %Responded
                if  (a==1 && s==8) || (a==3 && s==7) %offerOnset
                    dmat(end-1:end,:)=[];
                    dmat = dmat(REG.preResp(1:end-2)==1,:);
                    trial_data=trial_data(REG.preResp(1:end-2)==1,:);

                    dmat(isnan(REG.rewSd(REG.preResp(1:end-2)==1))|isnan(REG.pupil_iti(REG.preResp(1:end-2)==1)),:)=[];
                    trial_data(isnan(REG.rewSd(REG.preResp(1:end-2)==1))|isnan(REG.pupil_iti(REG.preResp(1:end-2)==1)),:)=[];

                else
                    dmat = dmat(REG.preResp==1,:);
                    trial_data=trial_data(REG.preResp==1,:);

                    dmat(isnan(REG.rewSd(REG.preResp==1))|isnan(REG.pupil_iti(REG.preResp==1)),:)=[];
                    trial_data(isnan(REG.rewSd(REG.preResp==1))|isnan(REG.pupil_iti(REG.preResp==1)),:)=[];
                end

                % contrast design matrix
                contrasts = diag(ones(size(dmat,2),1));

                % normalise data
                dmat(:,1:(size(dmat,2)-1)) = zscore(dmat(:,1:(size(dmat,2)-1)));

                % beta X time output
                betaOut = ols(trial_data,dmat,contrasts);
                allBeta.(sprintf('%s',feat_name)).(sprintf('%s',roi{r}))(:,:,z) = betaOut;

            elseif GLM == 4  %GLM3.2b (behavioural_history == reject)

                % load time-series data
                load([timecourse_dir,animal,'/',session,'/',feat_name,'/',([roi{r},'_epoched'])]);

                % regressor design matrix (don't forget constant!)
                dmat = [REG.rewMag, REG.rewProb, REG.rewAvg, REG.rewSd, REG.pupil_iti, REG.totalTime, REG.constant];

                %NotResponded
                if  (a==1 && s==8) || (a==3 && s==7) %offerOnset
                    dmat(end-1:end,:)=[];
                    dmat = dmat(REG.preResp(1:end-2)==0,:);
                    trial_data=trial_data(REG.preResp(1:end-2)==0,:);

                    dmat(isnan(REG.rewSd(REG.preResp(1:end-2)==0))|isnan(REG.pupil_iti(REG.preResp(1:end-2)==0)),:)=[];
                    trial_data(isnan(REG.rewSd(REG.preResp(1:end-2)==0))|isnan(REG.pupil_iti(REG.preResp(1:end-2)==0)),:)=[];
                else
                    dmat = dmat(REG.preResp==0,:);
                    trial_data=trial_data(REG.preResp==0,:);

                    dmat(isnan(REG.rewSd(REG.preResp==0))|isnan(REG.pupil_iti(REG.preResp==0)),:)=[];
                    trial_data(isnan(REG.rewSd(REG.preResp==0))|isnan(REG.pupil_iti(REG.preResp==0)),:)=[];
                end

                % contrast design matrix
                contrasts = diag(ones(size(dmat,2),1));

                % normalise data
                dmat(:,1:(size(dmat,2)-1)) = zscore(dmat(:,1:(size(dmat,2)-1)));

                % beta X time output
                betaOut = ols(trial_data,dmat,contrasts);
                allBeta.(sprintf('%s',feat_name)).(sprintf('%s',roi{r}))(:,:,z) = betaOut;

            elseif GLM == 5  %% GLM3.4

                load([timecourse_dir,animal,'/',session,'/',feat_name,'/',([roi{r},'_epoched'])]);

                k=find(REG.state_transition==1);
                pre_window = 3;
                post_window = 3;

                for i = 1:length(k)
                    range_start = max(1, k(i) - pre_window);    % Adjusted to stay within vector bounds
                    range_end = min(length(REG.state_transition), k(i) + post_window);
                    REG.state_transition(range_start:range_end) = 1; % Adjusted to stay within vector bounds
                end

                % regressor design matrix (don't forget constant!)
                dmat = [REG.state_transition, REG.viterbi_state, REG.pupil_iti, REG.totalTime, REG.constant];

                % contrast design matrix
                contrasts = diag(ones(size(dmat,2),1));

                % remove trials with no response
                dmat(isnan(REG.state_transition)|isnan(REG.viterbi_state)|isnan(REG.pupil_iti)|isnan(REG.response),:)=[];
                trial_data(isnan(REG.state_transition)|isnan(REG.viterbi_state)|isnan(REG.pupil_iti)|isnan(REG.response),:)=[];

                if    (a==1 && s==8) || (a==3 && s==7) %offerOnset
                    dmat(end-1:end,:)=[];
                end

                % normalise data
                dmat(:,1:(size(dmat,2)-1)) = zscore(dmat(:,1:(size(dmat,2)-1)));

                % beta X time output
                betaOut = ols(trial_data,dmat,contrasts);
                allBeta.(sprintf('%s',feat_name)).(sprintf('%s',roi{r}))(:,:,z) = betaOut;

            elseif GLM == 6  %% GLM3.4a (state_transition == high-to-low)

                load([timecourse_dir,animal,'/',session,'/',feat_name,'/',([roi{r},'_epoched'])]);

                k=find(REG.state_decrease==1);
                pre_window = 3;
                post_window = 3;

                for i = 1:length(k)
                    range_start = max(1, k(i) - pre_window);    % Adjusted to stay within vector bounds
                    range_end = min(length(REG.state_decrease), k(i) + post_window);
                    REG.state_decrease(range_start:range_end) = 1; % Adjusted to stay within vector bounds
                end

                % regressor design matrix (don't forget constant!)
                dmat = [REG.state_decrease, REG.pupil_iti, REG.totalTime, REG.constant];

                % contrast design matrix
                contrasts = diag(ones(size(dmat,2),1));

                % remove trials with no response
                dmat(isnan(REG.state_decrease)|isnan(REG.viterbi_state)|isnan(REG.pupil_iti)|isnan(REG.response),:)=[];
                trial_data(isnan(REG.state_decrease)|isnan(REG.viterbi_state)|isnan(REG.pupil_iti)|isnan(REG.response),:)=[];

                if    (a==1 && s==8) || (a==3 && s==7) %offerOnset
                    dmat(end-1:end,:)=[];
                end

                % normalise data
                dmat(:,1:(size(dmat,2)-1)) = zscore(dmat(:,1:(size(dmat,2)-1)));

                % beta X time output
                betaOut = ols(trial_data,dmat,contrasts);
                allBeta.(sprintf('%s',feat_name)).(sprintf('%s',roi{r}))(:,:,z) = betaOut;

            elseif GLM == 7  %% GLM3.4b (state_transition == low-to-high)

                load([timecourse_dir,animal,'/',session,'/',feat_name,'/',([roi{r},'_epoched'])]);

                k=find(REG.state_increase==1);

                for i = 1:length(k)
                    range_start = max(1, k(i) - 3);    % Adjusted to stay within vector bounds
                    range_end = min(length(REG.state_increase), k(i) + 3);
                    REG.state_increase(range_start:range_end) = 1; % Adjusted to stay within vector bounds
                end

                % regressor design matrix (don't forget constant!)
                dmat = [REG.state_increase, REG.pupil_iti, REG.totalTime, REG.constant];

                % contrast design matrix
                contrasts = diag(ones(size(dmat,2),1));

                % remove trials with no response
                dmat(isnan(REG.state_increase)|isnan(REG.viterbi_state)|isnan(REG.pupil_iti)|isnan(REG.response),:)=[];
                trial_data(isnan(REG.state_increase)|isnan(REG.viterbi_state)|isnan(REG.pupil_iti)|isnan(REG.response),:)=[];

                if    (a==1 && s==8) || (a==3 && s==7) %offerOnset
                    dmat(end-1:end,:)=[];
                end

                % normalise data
                dmat(:,1:(size(dmat,2)-1)) = zscore(dmat(:,1:(size(dmat,2)-1)));

                % beta X time output
                betaOut = ols(trial_data,dmat,contrasts);
                allBeta.(sprintf('%s',feat_name)).(sprintf('%s',roi{r}))(:,:,z) = betaOut;

            elseif GLM == 8  %% GLM3.5

                % load time-series data
                load([timecourse_dir,animal,'/',session,'/',feat_name,'/',([roi{r},'_epoched'])]);

                % regressor design matrix (don't forget constant!)
                dmat = [REG.ev, REG.evAvg, REG.pupil_iti, REG.totalTime, REG.constant];

                % contrast design matrix
                contrasts = diag(ones(size(dmat,2),1));

                % remove trials with missing data in the design matrix
                dmat(isnan(REG.pupil_iti)|isnan(REG.rewSd)|isnan(REG.preResp),:)=[];
                trial_data(isnan(REG.pupil_iti)|isnan(REG.rewSd)|isnan(REG.preResp),:)=[];

                if    (a==1 && s==8) || (a==3 && s==7) %offerOnset
                    dmat(end-1:end,:)=[];
                end

                % normalise data
                dmat(:,1:(size(dmat,2)-1)) = zscore(dmat(:,1:(size(dmat,2)-1)));

                % beta X time output
                betaOut = ols(trial_data,dmat,contrasts);
                allBeta.(sprintf('%s',feat_name)).(sprintf('%s',roi{r}))(:,:,z) = betaOut;

            elseif GLM == 9 % GLM3.8

                % load time-series data
                seed_TC = load([timecourse_dir,animal,'/',session,'/',feat_name,'/',([seed_roi{1},'_epoched'])]);
                roi_TC  = load([timecourse_dir,animal,'/',session,'/',feat_name,'/',([roi{r},'_epoched'])]);
                seed_TC.trial_data(isnan(REG.response)|isnan(REG.pupil_iti),:)=[];

                % run GLM
                for i = 1:size(seed_TC.trial_data,2)

                    % load ROI time-series data
                    REG.TC  = roi_TC.trial_data(:,i);

                    dmat = [REG.response, REG.totalTime, REG.constant];

                    if    (a==1 && s==8) || (a==3 && s==7) %offerOnset
                        dmat(end-1:end,:)=[];
                    end
                    dmat = [REG.TC, dmat];

                    % remove trials with no response
                    dmat(isnan(REG.response)|isnan(REG.pupil_iti),:)=[];

                    % normalise data
                    dmat(:,1:(size(dmat,2)-1)) = zscore(dmat(:,1:(size(dmat,2)-1)));

                    % create PPI regressor
                    REG.PPI = zscore (dmat(:,1) .* dmat(:,2));
                    dmat = [REG.PPI, dmat];

                    % contrast design matrix
                    contrasts = diag(ones(size(dmat,2),1));

                    % beta X time output
                    betaOut(:,i) = ols(seed_TC.trial_data(:,i),dmat,contrasts);
                    clear dmat contrasts REG.TC REG.PPI

                end

                allBeta.(sprintf('%s',feat_name)).(sprintf('%s',roi{r}))(:,:,z) = betaOut;

            elseif GLM == 10 % GLM 3.6a (behavioural_history == response)

                % load time-series data
                seed_TC = load([timecourse_dir,animal,'/',session,'/',feat_name,'/',([seed_roi{1},'_epoched'])]);
                roi_TC  = load([timecourse_dir,animal,'/',session,'/',feat_name,'/',([roi{r},'_epoched'])]);

                %Responded
                if  (a==1 && s==8) || (a==3 && s==7) %offerOnset

                    seed_TC.trial_data=seed_TC.trial_data(REG.preResp(1:end-2)==1,:);
                    seed_TC.trial_data(isnan(REG.rewSd(REG.preResp(1:end-2)==1))|isnan(REG.pupil_iti(REG.preResp(1:end-2)==1)),:)=[];

                    for i = 1:size(seed_TC.trial_data,2)

                        % load ROI time-series data
                        REG.TC  = roi_TC.trial_data(:,i);
                        % regressor design matrix (don't forget constant!)
                        dmat = [REG.rewMag, REG.rewProb, REG.rewAvg, REG.rewSd, REG.pupil_iti, REG.totalTime, REG.constant];

                        dmat(end-1:end,:)=[];
                        dmat = dmat(REG.preResp(1:end-2)==1,:);
                        REG.TC=REG.TC(REG.preResp(1:end-2)==1,:);

                        dmat(isnan(REG.rewSd(REG.preResp(1:end-2)==1))|isnan(REG.pupil_iti(REG.preResp(1:end-2)==1)),:)=[];
                        REG.TC(isnan(REG.rewSd(REG.preResp(1:end-2)==1))|isnan(REG.pupil_iti(REG.preResp(1:end-2)==1)),:)=[];

                        dmat = [REG.TC, dmat];

                        % normalise data
                        dmat(:,1:(size(dmat,2)-1)) = zscore(dmat(:,1:(size(dmat,2)-1)));

                        % create PPI regressor
                        REG.PPI = zscore (dmat(:,1) .* dmat(:,4));
                        dmat = [REG.PPI, dmat];

                        % contrast design matrix
                        contrasts = diag(ones(size(dmat,2),1));

                        % beta X time output
                        betaOut(:,i) = ols(seed_TC.trial_data(:,i),dmat,contrasts);
                        clear dmat contrasts REG.TC REG.PPI
                    end
                else
                    seed_TC.trial_data=seed_TC.trial_data(REG.preResp==1,:);
                    seed_TC.trial_data(isnan(REG.rewSd(REG.preResp==1))|isnan(REG.pupil_iti(REG.preResp==1)),:)=[];
                    for i = 1:size(seed_TC.trial_data,2)

                        % load ROI time-series data
                        REG.TC  = roi_TC.trial_data(:,i);
                        % regressor design matrix (don't forget constant!)
                        dmat = [REG.rewMag, REG.rewProb, REG.rewAvg, REG.rewSd, REG.pupil_iti, REG.totalTime, REG.constant];

                        dmat = dmat(REG.preResp==1,:);
                        REG.TC=REG.TC(REG.preResp==1,:);

                        dmat(isnan(REG.rewSd(REG.preResp==1))|isnan(REG.pupil_iti(REG.preResp==1)),:)=[];
                        REG.TC(isnan(REG.rewSd(REG.preResp==1))|isnan(REG.pupil_iti(REG.preResp==1)),:)=[];

                        dmat = [REG.TC, dmat];

                        % normalise data
                        dmat(:,1:(size(dmat,2)-1)) = zscore(dmat(:,1:(size(dmat,2)-1)));

                        % create PPI regressor
                        REG.PPI = zscore (dmat(:,1) .* dmat(:,4));
                        dmat = [REG.PPI, dmat];

                        % contrast design matrix
                        contrasts = diag(ones(size(dmat,2),1));

                        % beta X time output
                        betaOut(:,i) = ols(seed_TC.trial_data(:,i),dmat,contrasts);
                        clear dmat contrasts REG.TC REG.PPI
                    end
                end
                allBeta.(sprintf('%s',feat_name)).(sprintf('%s',roi{r}))(:,:,z) = betaOut;

            elseif GLM == 11 % GLM3.6b (behavioural_history == reject)
                % load time-series data
                seed_TC = load([timecourse_dir,animal,'/',session,'/',feat_name,'/',([seed_roi{1},'_epoched'])]);
                roi_TC  = load([timecourse_dir,animal,'/',session,'/',feat_name,'/',([roi{r},'_epoched'])]);

                %Responded
                if  (a==1 && s==8) || (a==3 && s==7) %offerOnset

                    seed_TC.trial_data=seed_TC.trial_data(REG.preResp(1:end-2)==0,:);
                    seed_TC.trial_data(isnan(REG.rewSd(REG.preResp(1:end-2)==0))|isnan(REG.pupil_iti(REG.preResp(1:end-2)==0)),:)=[];

                    for i = 1:size(seed_TC.trial_data,2)

                        % load ROI time-series data
                        REG.TC  = roi_TC.trial_data(:,i);
                        % regressor design matrix (don't forget constant!)
                        dmat = [REG.rewMag, REG.rewProb, REG.rewAvg, REG.rewSd, REG.pupil_iti, REG.totalTime, REG.constant];

                        dmat(end-1:end,:)=[];
                        dmat = dmat(REG.preResp(1:end-2)==0,:);
                        REG.TC=REG.TC(REG.preResp(1:end-2)==0,:);

                        dmat(isnan(REG.rewSd(REG.preResp(1:end-2)==0))|isnan(REG.pupil_iti(REG.preResp(1:end-2)==0)),:)=[];
                        REG.TC(isnan(REG.rewSd(REG.preResp(1:end-2)==0))|isnan(REG.pupil_iti(REG.preResp(1:end-2)==0)),:)=[];

                        dmat = [REG.TC, dmat];

                        % normalise data
                        dmat(:,1:(size(dmat,2)-1)) = zscore(dmat(:,1:(size(dmat,2)-1)));

                        % create PPI regressor
                        REG.PPI = zscore (dmat(:,1) .* dmat(:,4));
                        dmat = [REG.PPI, dmat];

                        % contrast design matrix
                        contrasts = diag(ones(size(dmat,2),1));

                        % beta X time output
                        betaOut(:,i) = ols(seed_TC.trial_data(:,i),dmat,contrasts);
                        clear dmat contrasts REG.TC REG.PPI
                    end
                else
                    seed_TC.trial_data=seed_TC.trial_data(REG.preResp==0,:);
                    seed_TC.trial_data(isnan(REG.rewSd(REG.preResp==0))|isnan(REG.pupil_iti(REG.preResp==0)),:)=[];
                    for i = 1:size(seed_TC.trial_data,2)

                        % load ROI time-series data
                        REG.TC  = roi_TC.trial_data(:,i);
                        % regressor design matrix (don't forget constant!)
                        dmat = [REG.rewMag, REG.rewProb, REG.rewAvg, REG.rewSd, REG.pupil_iti, REG.totalTime, REG.constant];

                        dmat = dmat(REG.preResp==0,:);
                        REG.TC=REG.TC(REG.preResp==0,:);

                        dmat(isnan(REG.rewSd(REG.preResp==0))|isnan(REG.pupil_iti(REG.preResp==0)),:)=[];
                        REG.TC(isnan(REG.rewSd(REG.preResp==0))|isnan(REG.pupil_iti(REG.preResp==0)),:)=[];

                        dmat = [REG.TC, dmat];

                        % normalise data
                        dmat(:,1:(size(dmat,2)-1)) = zscore(dmat(:,1:(size(dmat,2)-1)));

                        % create PPI regressor
                        REG.PPI = zscore (dmat(:,1) .* dmat(:,4));
                        dmat = [REG.PPI, dmat];

                        % contrast design matrix
                        contrasts = diag(ones(size(dmat,2),1));

                        % beta X time output
                        betaOut(:,i) = ols(seed_TC.trial_data(:,i),dmat,contrasts);
                        clear dmat contrasts REG.TC REG.PPI
                    end
                end
                allBeta.(sprintf('%s',feat_name)).(sprintf('%s',roi{r}))(:,:,z) = betaOut;

            elseif GLM == 12 % GLM3.7a (PPI as a function of high-to-low state-transition)
                
                k=find(REG.state_decrease==1);
                pre_window = 3;
                post_window = 3;
                
                for j = 1:length(k)
                    range_start = max(1, k(j) - pre_window);    % Adjusted to stay within vector bounds
                    range_end = min(length(REG.state_transition), k(j) + post_window);
                    REG.state_decrease(range_start:range_end) = 1; % Adjusted to stay within vector bounds
                end

                k=find(REG.state_increase==1);
                pre_window = 3;
                post_window = 3;
                
                for j = 1:length(k)
                    range_start = max(1, k(j) - pre_window);    % Adjusted to stay within vector bounds
                    range_end = min(length(REG.state_transition), k(j) + post_window);
                    REG.state_increase(range_start:range_end) = 1; % Adjusted to stay within vector bounds
                end
                
                % load time-series data
                seed_TC = load([timecourse_dir,animal,'/',session,'/',feat_name,'/',([seed_roi{1},'_epoched'])]);
                roi_TC  = load([timecourse_dir,animal,'/',session,'/',feat_name,'/',([roi{r},'_epoched'])]);
                seed_TC.trial_data(isnan(REG.state_decrease)|isnan(REG.state_increase)|isnan(REG.viterbi_state)|isnan(REG.pupil_iti),:)=[];
                
                % run GLM
                for i = 1:size(seed_TC.trial_data,2)
                    
                    % load ROI time-series data
                    REG.TC  = roi_TC.trial_data(:,i);
                    
                    dmat = [REG.state_decrease, REG.state_increase, REG.viterbi_state, REG.pupil_iti, REG.totalTime, REG.constant];
                    
                    if    (a==1 && s==8) || (a==3 && s==7) %offerOnset
                        dmat(end-1:end,:)=[];
                    end
                    
                    dmat = [REG.TC, dmat];
                    dmat(isnan(REG.pupil_iti)|isnan(REG.state_decrease)|isnan(REG.state_increase)|isnan(REG.viterbi_state),:)=[];
                    
                    % normalise data
                    dmat(:,1:(size(dmat,2)-1)) = zscore(dmat(:,1:(size(dmat,2)-1)));
                    
                    % create PPI regressor
                    REG.PPI = zscore (dmat(:,1) .* dmat(:,2));
                    dmat = [REG.PPI, dmat];
                    
                    % contrast design matrix
                    contrasts = diag(ones(size(dmat,2),1));
                    
                    % beta X time output
                    betaOut(:,i) = ols(seed_TC.trial_data(:,i),dmat,contrasts);
                    clear dmat contrasts REG.TC REG.PPI
                end
                allBeta.(sprintf('%s',feat_name)).(sprintf('%s',roi{r}))(:,:,z) = betaOut;
                       
            elseif GLM == 13 % GLM3.7b (PPI as a function of low-to-high state-transition)
                
                k=find(REG.state_decrease==1);
                pre_window = 3;
                post_window = 3;
                
                for j = 1:length(k)
                    range_start = max(1, k(j) - pre_window);    % Adjusted to stay within vector bounds
                    range_end = min(length(REG.state_transition), k(j) + post_window);
                    REG.state_decrease(range_start:range_end) = 1; % Adjusted to stay within vector bounds
                end

                k=find(REG.state_increase==1);
                pre_window = 3;
                post_window = 3;
                
                for j = 1:length(k)
                    range_start = max(1, k(j) - pre_window);    % Adjusted to stay within vector bounds
                    range_end = min(length(REG.state_transition), k(j) + post_window);
                    REG.state_increase(range_start:range_end) = 1; % Adjusted to stay within vector bounds
                end
                
                % load time-series data
                seed_TC = load([timecourse_dir,animal,'/',session,'/',feat_name,'/',([seed_roi{1},'_epoched'])]);
                roi_TC  = load([timecourse_dir,animal,'/',session,'/',feat_name,'/',([roi{r},'_epoched'])]);
                seed_TC.trial_data(isnan(REG.state_decrease)|isnan(REG.state_increase)|isnan(REG.viterbi_state)|isnan(REG.pupil_iti),:)=[];
                
                % run GLM
                for i = 1:size(seed_TC.trial_data,2)
                    
                    % load ROI time-series data
                    REG.TC  = roi_TC.trial_data(:,i);
                    
                    dmat = [REG.state_decrease, REG.state_increase, REG.viterbi_state, REG.pupil_iti, REG.totalTime, REG.constant];
                    
                    if    (a==1 && s==8) || (a==3 && s==7) %offerOnset
                        dmat(end-1:end,:)=[];
                    end
                    
                    dmat = [REG.TC, dmat];
                    dmat(isnan(REG.pupil_iti)|isnan(REG.state_decrease)|isnan(REG.state_increase)|isnan(REG.viterbi_state),:)=[];
                    
                    % normalise data
                    dmat(:,1:(size(dmat,2)-1)) = zscore(dmat(:,1:(size(dmat,2)-1)));
                    
                    % create PPI regressor
                    REG.PPI = zscore (dmat(:,1) .* dmat(:,3));
                    dmat = [REG.PPI, dmat];
                    
                    % contrast design matrix
                    contrasts = diag(ones(size(dmat,2),1));
                    
                    % beta X time output
                    betaOut(:,i) = ols(seed_TC.trial_data(:,i),dmat,contrasts);
                    clear dmat contrasts REG.TC REG.PPI
                end
                allBeta.(sprintf('%s',feat_name)).(sprintf('%s',roi{r}))(:,:,z) = betaOut;
            
            elseif GLM == 14 % GLM3.7c (PPI as a function of state level)
                
                k=find(REG.state_transition==1);
                pre_window = 3;
                post_window = 3;
                
                for j = 1:length(k)
                    range_start = max(1, k(j) - pre_window);    % Adjusted to stay within vector bounds
                    range_end = min(length(REG.state_transition), k(j) + post_window);
                    REG.state_transition(range_start:range_end) = 1; % Adjusted to stay within vector bounds
                end
                
                % load time-series data
                seed_TC = load([timecourse_dir,animal,'/',session,'/',feat_name,'/',([seed_roi{1},'_epoched'])]);
                roi_TC  = load([timecourse_dir,animal,'/',session,'/',feat_name,'/',([roi{r},'_epoched'])]);
                seed_TC.trial_data(isnan(REG.state_transition)|isnan(REG.viterbi_state)|isnan(REG.pupil_iti),:)=[];
                
                % run GLM
                for i = 1:size(seed_TC.trial_data,2)
                    
                    % load ROI time-series data
                    REG.TC  = roi_TC.trial_data(:,i);
                    
                    dmat = [REG.state_transition, REG.viterbi_state, REG.pupil_iti, REG.totalTime, REG.constant];
                    
                    if    (a==1 && s==8) || (a==3 && s==7) %offerOnset
                        dmat(end-1:end,:)=[];
                    end
                    
                    dmat = [REG.TC, dmat];
                    dmat(isnan(REG.pupil_iti)|isnan(REG.state_transition)|isnan(REG.viterbi_state),:)=[];
                    
                    % normalise data
                    dmat(:,1:(size(dmat,2)-1)) = zscore(dmat(:,1:(size(dmat,2)-1)));
                    
                    % create PPI regressor
                    REG.PPI = zscore (dmat(:,1) .* dmat(:,3));
                    dmat = [REG.PPI, dmat];
                    
                    % contrast design matrix
                    contrasts = diag(ones(size(dmat,2),1));
                    
                    % beta X time output
                    betaOut(:,i) = ols(seed_TC.trial_data(:,i),dmat,contrasts);
                    clear dmat contrasts REG.TC REG.PPI
                end
                allBeta.(sprintf('%s',feat_name)).(sprintf('%s',roi{r}))(:,:,z) = betaOut;
            end
        end
    end
end


%% plot group data

contrast_to_plot = 1;

pre_win = 1;
post_win = 5;
upsample = 15;
TR = 1.282;
nsamples = round(((pre_win+post_win)./TR)*upsample);
time = -pre_win:(pre_win+post_win)/nsamples:post_win;


for r = 1:numel(roi)

    beta = squeeze(allBeta.(sprintf('%s',feat_name)).(sprintf('%s',roi{r}))(contrast_to_plot,:,:));

    if writedata
        tc_data = {};
        tc_data(:, 1) = num2cell(mean(beta'));
        tc_data(:, 2) = num2cell(std(beta')/sqrt(size(beta',1)));
        tc_data(:, 3) = roi(r);
        tc_data(:, 4) = {['GLM_0', num2str(GLM)]};
        tc_data(:, 5) = {['contr_', num2str(contrast_to_plot)]};
        if GLM == 9||GLM == 10||GLM == 11
            filename = [tc_graph_dir, '/GLM_0', num2str(GLM),'_contr_', num2str(contrast_to_plot),'_',seed_roi{1},'_to_',roi{r},'_PPI.csv'];
        else
            filename = [tc_graph_dir, '/GLM_0', num2str(GLM),'_contr_', num2str(contrast_to_plot),'_',roi{r},'.csv'];
        end
        writetable(cell2table(tc_data), filename);
        clear tc_data
    end

end


%% stats

clearvars window LOOT peak
numsession = 1:59;
numsamples = nsamples;
start=24;
win_end=70;

reg = contrast_to_plot;

% find peak in a specified window
for r = 1:numel(roi)

    for s= 1:length(numsession)

        numsession(s) = [];
        window  = mean(squeeze(allBeta.(sprintf('%s',feat_name)).(sprintf('%s',roi{r}))(reg,:,numsession)),2);

        [m,i]  = max(abs(window(start:win_end)));
        i      = i + (start-1);
        LocPeak(s,r)=(i*6)/length(window);
        peak(s,r) = allBeta.(sprintf('%s',feat_name)).(sprintf('%s',roi{r}))(reg,i,s);

        numsession = 1:59;
        clear window fw beta

    end
    [h,p(r),ci,stats]= ttest(peak(:,r));
    t(r) = stats.tstat;
end

% Holm Bonferroni test
[bonf_p, bonf_h] = bonf_holm(p)
roi

if writedata
    tb = num2cell(peak);
    tb(:, size(tb,2)+1) = {['GLM_0', num2str(GLM)]};
    tb(:, size(tb,2)+1) = {['contr', num2str(contrast_to_plot)]};
    tb = cell2table(tb, 'VariableNames', [roi, 'GLM', 'contrast']);
    writetable(tb,[tc_graph_dir, '/GLM_0', num2str(GLM),'_contr_', num2str(contrast_to_plot),'_peaks.csv']);
    clear tc_data
end