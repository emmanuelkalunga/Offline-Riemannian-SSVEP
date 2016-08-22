%**************************************************************************
%File: GIS2015
%Main file for GIS2015 results
%Description: Compare the classification prformance while data filtered  
%with FIR vs IRR
%Author: Emmanuel K. Kalunga
%**************************************************************************
tLen = 4;
delay = 2;
filters = [1,2]; %1 for fir and 2 for iir
ac = zeros(12, 5, numel(filters));
%load filters
%Hd = {Hd_fir;Hd_iir};
for f =filters
for sub = 6:17
    clear x_all H_all P X Pm PSD
    %% Load data
    disp('********************************************************');
    disp(['Load data subject ', num2str(sub)]);
    [S_all, H_all] = loaddata(sub); %Returns cells of data from all available sessions
    Fs = H_all{1}.SampleRate;
    nbrSessions = length(S_all);
    sessions = 1:nbrSessions;
    %% Preprocessing of all available sessions (Same for training and test data)
    disp(['Preprocessing subject ', num2str(sub)]);
    % 1) Band pass filter and return super trials
    for session = 1:nbrSessions
         x_all{session} = bandpass_ext([12.95 13.05], [16.9 17.1], [20.9 21.1], S_all{session}, H_all{session},f);
    end
    % 2) Rearange data per trial
    X = get_trials(x_all, H_all, tLen, delay);
    %% Covariance matrices of all trials
    disp(['Covariance matrices for subject ', num2str(sub)]);
    Nt = size(X{1},3); %Number of trial
    for k = 1:Nt %loop for evrey trial
        for cl = 1:4
            P{cl}(:,:,k) = shcovft((X{cl}(:,:,k))'); % J. Schaefer Shrinkage covariance from Barachant toolbox
        end
    end
    disp('---------------------------------------------------------------')
    disp(['Start cross validation subject ', num2str(sub)]);
    for testSession = 1:nbrSessions
        disp(['Validation of session ', num2str(testSession)]);
        trials = 1:size(P{1},3);
        trialPerSession = size(P{1},3)/nbrSessions;
        testTrials = (trialPerSession*testSession-trialPerSession+1):(trialPerSession*testSession);
        trainTrials = setxor(trials, testTrials);
        %% GET COV MAT AND CORRESPONDING LABELS FOR TRAINING DATA
        trainSessions = setxor(sessions, testSession);
        COVtrain = cat(3, P{1}(:,:,trainTrials), P{2}(:,:,trainTrials), P{3}(:,:,trainTrials), P{4}(:,:,trainTrials));
        Ytrain = [zeros(1,length(trainTrials)) ones(1,length(trainTrials)) 3*ones(1,length(trainTrials)) 2*ones(1,length(trainTrials))];

        %%                  EVALUATION PHASE                               **
        %********************************************************************
        labels = [zeros(1,trialPerSession) ones(1, trialPerSession) 3*ones(1, trialPerSession) 2*ones(1, trialPerSession)];
        COVtest = cat(3, P{1}(:,:,testTrials), P{2}(:,:,testTrials), P{3}(:,:,testTrials), P{4}(:,:,testTrials));
        % Classification by Remannian Distance
        [Ytest d C] = mdm(COVtest,COVtrain,Ytrain);
        ac(sub-5, testSession, f) = sum((labels-Ytest)==0)/(trialPerSession*4);
    end
end
end
subId = zeros(1,size(ac,1));
subNbrOfSess = zeros(1,size(ac,1));
subAcMean = zeros(1,size(ac,1));
subVar = zeros(1,size(ac,1));
resMatrix = zeros(12,4,numel(filters));
resMean = zeros(numel(filters),3);
for f = 1:numel(filters)
    for i = 1:size(ac,1)
        acSi = ac(i,:,f);
        acSi = acSi(acSi~=0);
        subId(i) = i+5;
        subNbrOfSess(i) = length(acSi);
        subAcMean(i) = mean(acSi);
        subVar(i) = var(acSi);
    end
    resMatrix(:,:,f) = [subId' subNbrOfSess' subAcMean' subVar'];
    resMeanTmp = mean(resMatrix(:,:,f));
    resMeanTpm(2) = sum(resMatrix(:,2,f));
    resMean(f,:) = resMeanTmp(2:end);
end