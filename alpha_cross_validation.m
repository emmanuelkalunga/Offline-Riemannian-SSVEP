%**************************************************************************
%File: GIS2015
%Main file for GIS2015 results
%Description: Compare the classification prformance of alpha-divergence 
%(log determinant) with diferent alpha values
%Author: Emmanuel K. Kalunga
%**************************************************************************
alpha = -1:0.1:1;
tLen = 4;
delay = 2;
ac = zeros(12, 5, numel(alpha));
e = zeros(12, 5, numel(alpha));
t1 = cputime;
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
        x_all{session} = bandpass_filter_ext([12.95 13.05], [16.9 17.1], [20.9 21.1], S_all{session}, H_all{session});
    end
    % 2) Rearange data per trial
    X = get_trials(x_all, H_all, tLen, delay);
    %% Covariance matrices of all trials
    disp(['Covariance matrices ', num2str(sub)]);
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
        
        for a = 1:length(alpha)
            disp(['Alpha: ' num2str(alpha(a))]);
            t = cputime;
            [Ytest d C] = mdm_alpha(COVtest,COVtrain,Ytrain, 'ld', 'ld',alpha(a));
            e(sub-5, testSession, a) = cputime-t;
            ac(sub-5, testSession, a) = sum((labels-Ytest)==0)/(trialPerSession*4);
        end
    end
end
e_cross = cputime-t1;
subId = zeros(1,size(ac,1));
subNbrOfSess = zeros(1,size(ac,1));
subAcMean = zeros(1,size(ac,1));
subVar = zeros(1,size(ac,1));
resMatrix = zeros(12,5,numel(alpha));
resMean = zeros(numel(alpha),4);
for a = 1:numel(alpha)
    for i = 1:size(ac,1)
        acSi = ac(i,:,a);
        acSi = acSi(acSi~=0);
        subId(i) = i+5;
        subNbrOfSess(i) = length(acSi);
        subAcMean(i) = mean(acSi);
        subVar(i) = var(acSi);
        
        timeSi = e(i,:,a);
        timeSi = timeSi(timeSi~=0);
        subTimeMean(i) = mean(timeSi);
    end
    resMatrix(:,:,a) = [subId' subNbrOfSess' subAcMean' subVar' subTimeMean'];
    resMeanTmp = mean(resMatrix(:,:,a));
    resMeanTpm(2) = sum(resMatrix(:,2,a));
    resMean(a,:) = resMeanTmp(2:end);
end
[~, alphaBestIndex] = max(resMean(:,2));
alphaBest = alpha(alphaBestIndex);

%--LOAD value stored in mat file here 
%-- Add standard devition on the plot-------------------------------------
for a = 1:numel(alpha)
    resStdTmp = std(resMatrix(:,:,a));
    resStdTpm(2) = sum(resMatrix(:,2,a));
    resStd(a,:) = resStdTmp(2:end);
end
%--------------------------------------------------------------------------

% figure, plot(alpha, resMean(:,2),'--s')
% xlabel('alpha values (\alpha) ')
% ylabel('Classification accuracy')
% set(gca,'FontSize',12,'fontWeight','normal')
% set(findall(gcf,'type','text'),'FontSize',12,'fontWeight','normal')
% 
% figure, plot(alpha, resMean(:,3),'--s')
% xlabel('alpha values (\alpha) ')
% ylabel('cpu time')
% set(gca,'FontSize',12,'fontWeight','normal')
% set(findall(gcf,'type','text'),'FontSize',12,'fontWeight','normal')
figure
ylim([20 100])
[hAx,hLine1,hLine2] = plotyy(alpha,(round(1000*resMean(:,2)))/10,alpha,resMean(:,4),'plot');
set(get(hAx(1),'Ylabel'),'String','Accuracy (%)') 
set(get(hAx(2),'Ylabel'),'String','CPU time (s)') 
set(hLine1,'LineStyle','--','Marker', 's','LineWidth', 3)
set(hLine2,'LineStyle',':','Marker', 'o','LineWidth', 3)
xlabel('Alpha values (\alpha) ')
set(gca,'FontSize',14,'fontWeight','normal')
set(hAx(2),'FontSize',14,'fontWeight','normal')
set(findall(gcf,'type','text'),'FontSize',14,'fontWeight','normal')
%--PLot error bars
hold(hAx(1),'on')
errorbar(hAx(1), alpha, (round(1000*resMean(:,2)))/10, (round(1000*resStd(:,2)))/10, 's' )
hold(hAx(1),'off')
hold(hAx(2),'on')
errorbar(hAx(2), alpha, resMean(:,4), resStd(:,4), 'o' )
ylim(hAx(1),[20 max((round(1000*resMean(:,2)))/10+ (round(1000*resStd(:,2)))/10)])
% ylim(hAx(2),[0 max(resMean(:,4)+ resStd(:,4))])
ylim(hAx(2),[0 0.8])

%spaceplots; % Does not work with plotyy
% save('alpha_cross_final.mat', 'resMatrix','resMean','alpha','alphaBest')