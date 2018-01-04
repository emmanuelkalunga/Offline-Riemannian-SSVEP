%**************************************************************************
%File: GIS2015
%Main file for GIS2015 results
%Description: Compare the classification prformance of riemannian metrics,
%euclidean, log-euclidean, and alpha-divergence (log determinant)
%Author: Emmanuel K. Kalunga
%**************************************************************************
% methodMean = {'arithmetic','harmonic','geometric','riemann','logeuclid','kullback','sdivergence','ld','bhat','opttransp'};
% methodDist = {'euclid','euclid','euclid','riemann','logeuclid','kullback','sdivergence','ld','bhat','opttransp'};
methodMean = {'arithmetic','riemann','logeuclid','sdivergence'};
methodDist = {'euclid','riemann','logeuclid','sdivergence'};
tLen = 4;
delay = 2;
ac = zeros(12, 5, length(methodMean));
e = zeros(12, 5, length(methodMean));
for sub = 6:17
    clear x_all H_all P X Pm PSD
    %% Load data
    disp('********************************************************');
    disp(['Load data subject ', num2str(sub)]);
    [S_all, H_all] = loaddata(sub);%Returns cells of data from all available sessions
    Fs = H_all{1}.SampleRate;
    nbrSessions = length(S_all);
    sessions = 1:nbrSessions;
     
    trow = zeros(1,8);
    tcol = zeros(1,8);
    %tcol(1) = 0.2;
    %trow(1) = 0.2;
    tcol(2) = 0.3;
    tcol(3) = 0.3;
    tcol(5) = 0.5;
    trow = tcol;
    
    shiftMatrix = toeplitz(tcol, trow); %Shit (transformation) matrix
    
    for sess = sessions
        S_all{sess} = (shiftMatrix * S_all{sess}')' * shiftMatrix';
        S_all{sess} = S_all{sess}(:,2:end);
    end 
    
    %% Preprocessing of all available sessions (Same for training and test data)
    disp(['Preprocessing subject ', num2str(sub)]);
    % 1) Band pass filter and return super trials
    for session = 1:nbrSessions
        x_all{session} = bandpass_filter_ext([12.95 13.05], [16.9 17.1], [20.9 21.1], S_all{session}, H_all{session});
    end
    % 2) Rearange data per trial
    X = get_trials(x_all, H_all, tLen, delay);
    %% Covariance matrices of all trials
    Nt = size(X{1},3); %Number of trial
%     trow = zeros(1,24);
%     tcol = zeros(1,24);
%     tcol(1) = 0.2;
%     trow(1) = 0.2;
%     tcol(2:3) = 0.5;
%     shiftMatrix = toeplitz(tcol, trow); %Shit (transformation) matrix
    for k = 1:Nt %loop for evrey trial
        for cl = 1:4
            %P_trans{cl}(:,:,k) = shcovft((shiftMatrix'*X{cl}(:,:,k))'); % J. Schaefer Shrinkage covariance from Barachant toolbox
            P{cl}(:,:,k) = shcovft((X{cl}(:,:,k))'); % J. Schaefer Shrinkage covariance from Barachant toolbox
        end
    end
    %% Get determinant of all trials per class for each subject and the determinant of means computed with defferent means. And plot the averages covmat
    % THIS IS DONE ONLY FOR ANALYSIS PURPURSES (CAN BE TAKEN TO A DIFFERENT
    % FILE
    %-- Trials determinants
    for k = 1:Nt %loop for evrey trial
        for cl = 1:4
            determinants(k,cl,sub-5) = det(P{cl}(:,:,k));
        end
    end
    %-- Means determinant
    COVMAT = cat(3, P{:});
    LABELS = [zeros(1,Nt) ones(1,Nt) 3*ones(1,Nt) 2*ones(1,Nt)];
    labels = unique(LABELS);
    for method = 1:length(methodMean)
        for cl = 1:4
            Cm{sub-5}(:,:,cl,method) = mean_covariances_alpha(COVMAT(:,:,LABELS==labels(cl)),methodMean{method},0.6); %alpha = 0.6
            detMean(cl,method,sub-5) = det(squeeze(Cm{sub-5}(:,:,cl,method)));
        end
    end
    %-- Plot mean cov
    if (sub-5) == 10
    for method = 1:length(methodMean)
        figure
        for cl = 1:4
            subplot(2,2,cl)
            imagesc(squeeze(Cm{sub-5}(:,:,cl,method)));colormap(flipud(hot));
            title(['method ' methodMean{method} ', class' num2str(cl)])
            set(gca,'FontSize',12,'fontWeight','normal')
            set(findall(gcf,'type','text'),'FontSize',12,'fontWeight','normal')
            set(gca, 'XTick', []);
            set(gca, 'YTick', []);
        end
        %spaceplots;
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
%         COVtrain = cat(3, P_trans{1}(:,:,trainTrials), P_trans{2}(:,:,trainTrials), P_trans{3}(:,:,trainTrials), P_trans{4}(:,:,trainTrials));
        
%         for mat = 1:size(COVtrain,3)
%             COVtrain(:,:,mat) = shiftMatrix*COVtrain(:,:,mat)*shiftMatrix';
%         end
        Ytrain = [zeros(1,length(trainTrials)) ones(1,length(trainTrials)) 3*ones(1,length(trainTrials)) 2*ones(1,length(trainTrials))];

        %%                  EVALUATION PHASE                               **
        %********************************************************************
        labels = [zeros(1,trialPerSession) ones(1, trialPerSession) 3*ones(1, trialPerSession) 2*ones(1, trialPerSession)];
        COVtest = cat(3, P{1}(:,:,testTrials), P{2}(:,:,testTrials), P{3}(:,:,testTrials), P{4}(:,:,testTrials));
%         COVtest = cat(3, P_trans{1}(:,:,testTrials), P_trans{2}(:,:,testTrials), P_trans{3}(:,:,testTrials), P_trans{4}(:,:,testTrials));
%         for mat = 1:size(COVtest,3)
%             COVtest(:,:,mat) = shiftMatrix*COVtest(:,:,mat)*shiftMatrix';
%         end
        % Classification by Remannian Distance
        %--For optimal alpha-divergence
        alpha = 0.6; %Alpha is set to 0.6 through cross validation (alpha_cross_validation.m)
        %--For Bhattacharyya divergence and mean 
        %alpha = 0; %Bhattacharyya metrics are similar to alpha divergence (with log-determinant function) with alpha being 0
        
        %Priority(3); %Set highest priority
        for method = 1:length(methodMean)
            disp(['Method ' num2str(method)]);       
            if strcmp(methodMean{method}, 'kullback')
                alpha = 1;
                meanmethod = 'arithmetic';
                distmethod = 'kullback';
            elseif strcmp(methodMean{method}, 'bhat')
                alpha = 0;
                meanmethod = 'ld';
                distmethod = 'ld';
            else
                alpha = 0.6;
                meanmethod = methodMean{method};
                distmethod = methodDist{method};
            end
            t = cputime;
            [Ytest d C] = mdm_alpha(COVtest,COVtrain,Ytrain, meanmethod, distmethod, alpha);
            e(sub-5, testSession, method) = cputime-t;
            
            if e(sub-5, testSession, method) == 0
                e(sub-5, testSession, method) = 0.0001;
            end
            ac(sub-5, testSession, method) = sum((labels-Ytest)==0)/(trialPerSession*4);
        end
        %Priority(0); %Reset priority to default
    end
end
subId = zeros(1,size(ac,1));
subNbrOfSess = zeros(1,size(ac,1));
subAcMean = zeros(1,size(ac,1));
subTimeMean = zeros(1,size(ac,1));
subVar = zeros(1,size(ac,1));
resMatrix = zeros(12,5,length(methodMean));
resMean = zeros(length(methodMean),4);
for method = 1:length(methodMean)
    for i = 1:size(ac,1)
        acSi = ac(i,:,method);
        acSi = acSi(acSi~=0);
        subId(i) = i+5;
        subNbrOfSess(i) = length(acSi);
        subAcMean(i) = mean(acSi);
        subVar(i) = var(acSi);
        
        timeSi = e(i,:,method);
        timeSi = timeSi(timeSi~=0);
        subTimeMean(i) = mean(timeSi);
    end
    resMatrix(:,:,method) = [subId' subNbrOfSess' subAcMean' subVar' subTimeMean'];
    resMeanTmp = mean(resMatrix(:,:,method));
    resMeanTpm(2) = sum(resMatrix(:,2,method));
    resMean(method,:) = resMeanTmp(2:end);
end
% save('gsi.mat', 'resMatrix','resMean');

%resMatrix2 = resMatrix(:,:,[2 3 1 4]);
for method = 1:length(methodMean)
    resMeanTmp = mean(resMatrix(:,:,method));
    resMeanTpm(2) = sum(resMatrix(:,2,method));
    resMean2(method,:) = resMeanTmp(2:end);
end
% save('entropy.mat', 'resMean2', 'resMean', 'resMatrix');
figure
[hAx,hLine1,hLine2] = plotyy(1:length(methodMean),(round(1000*resMean2(:,2)))/10,1:length(methodMean),resMean2(:,4),'plot');
set(get(hAx(1),'Ylabel'),'String','Accuracy (%)') 
set(get(hAx(2),'Ylabel'),'String','cpu time (s)') 
set(hLine1,'LineStyle','--','Marker', 's')
set(hLine2,'LineStyle',':','Marker', 's')
xlabel('Method')
set(gca,'FontSize',12,'fontWeight','normal')
set(findall(gcf,'type','text'),'FontSize',12,'fontWeight','normal')
