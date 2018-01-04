%**************************************************************************
%File: main.m
%Main file for Riemannian diastances and divergences for ssvep eeg
%Description: Compare the classification prformance of riemannian metrics,
%euclidean, and cca based. Results are stored in variables resMatrix and resMean
%Author: Emmanuel K. Kalunga
%**************************************************************************

methodMean = {'arithmetic','harmonic','riemann','logeuclid','kullback','sdivergence','ld','bhat','wasserstein', 'jeffreys'};
methodDist = {'euclid','harmonic','riemann','logeuclid','kullback','sdivergence','ld','bhat','wasserstein', 'jeffreys'};

tLen = 4;
delay = 2;
ac = zeros(12, 5, length(methodMean));
e = zeros(12, 5, length(methodMean));
det_mean = zeros(12, 5, length(methodMean));
for sub = 6:17
    clear x_all H_all P X Pm PSD
    %% Load data
    disp('********************************************************');
    disp(['Load data subject ', num2str(sub)]);
    % Note: Ensure that the function loaddata points to the data directory
    [S_all, H_all] = loaddata(sub);%Returns cells of data from all available sessions
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
    disp(['Preprocessing subject ', num2str(sub)]);
    Nt = size(X{1},3); %Number of trial
    for k = 1:Nt %loop for evrey trial
        for cl = 1:4
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
    %-- Get distances of trial from riemannian mean
    COVMAT = cat(3, P{:});
    LABELS = [zeros(1,Nt) ones(1,Nt) 3*ones(1,Nt) 2*ones(1,Nt)];
    labels = unique(LABELS);
    
    for cl = 1:4
        covmat = COVMAT(:,:,LABELS==labels(cl));
        %covmat = COVMAT;
        riem_mean = riemann_mean(covmat);
        dist = [];
        for m = 1:size(covmat,3)
            dist(m) = distance_riemann(covmat(:,:,m),riem_mean);
        end
        geom_mean = exp(mean(log(dist)));
        arit_mean = mean(dist);
        wiener_2(cl, sub-5) = geom_mean/arit_mean;
    end
    %--Fonctionnel de l'optimisation de la moyenne
%     if (sub-5) == 12
%         methodOpt = methodMean([3 6 7 8 9]);
%         for method = 1:length(methodOpt)
%             if strcmp(methodOpt{method}, 'bhat')
%                 alpha = 0;
%                 meanmethod = 'ld';
%             else
%                 alpha = 0.6;
%                 meanmethod = methodOpt{method};
%             end
% 
%             [M,sum_dis{method}] = mean_analysis(COVMAT, alpha,meanmethod);
%             figure, plot(sum_dis{method}(2:end),'*-','LineWidth', 1,'DisplayName', 'trial det')  
% 
%             xlabel('Iterations')
%             ylabel('$ \sum d^2(\Sigma_i,\Sigma)$', 'Interpreter', 'latex', 'FontSize',14, 'fontWeight','normal');
%             set(gca,'FontSize',14,'fontWeight','normal')
%             set(findall(gcf,'type','text'),'FontSize',14,'fontWeight','normal')
%         end  
%     end
    
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
        Ytrain = [zeros(1,length(trainTrials)) ones(1,length(trainTrials)) 3*ones(1,length(trainTrials)) 2*ones(1,length(trainTrials))];

        %%                  EVALUATION PHASE                               **
        %********************************************************************
        labels = [zeros(1,trialPerSession) ones(1, trialPerSession) 3*ones(1, trialPerSession) 2*ones(1, trialPerSession)];
        COVtest = cat(3, P{1}(:,:,testTrials), P{2}(:,:,testTrials), P{3}(:,:,testTrials), P{4}(:,:,testTrials));
        % Classification by Remannian Distance
        %--For optimal alpha-divergence
        alpha = 0.6; %Alpha is set to 0.6 through cross validation (alpha_cross_validation.m)
        %--For Bhattacharyya divergence and mean 
        %alpha = 0; %Bhattacharyya metrics are similar to alpha divergence (with log-determinant function) with alpha being 0
        
        labels_unique = unique(labels);
        for cl = 1:4
            covmat = COVtrain(:,:,Ytrain==labels_unique(cl));
            %covmat = COVMAT;
            riem_mean = riemann_mean(covmat);
            dist = [];
            for m = 1:size(covmat,3)
                dist(m) = distance_riemann(covmat(:,:,m),riem_mean);
            end
            geom_mean = exp(mean(log(dist)));
            arit_mean = mean(dist);
            wiener_3(sub-5,testSession,cl) = geom_mean/arit_mean; 
        end
        
        
        %Priority(3); %Set highest priority
        for method = 1:length(methodMean)
            disp(['Method ' num2str(method)]);       
            if strcmp(methodMean{method}, 'kullback')
                alpha = 1;
                %meanmethod = 'arithmetic';
                meanmethod = 'kullback';
                distmethod = 'kullback';
            elseif strcmp(methodMean{method}, 'bhat')
                alpha = 0;
                meanmethod = 'ld';
                distmethod = 'ld';
            elseif strcmp(methodMean{method}, 'mahalanobis')
                alpha = 1; %don't care
                meanmethod = 'arithmetic';
                distmethod = 'mahalanobis';
	    %elseif strcmp(methodMean{method}, 'wasserstein')
	    %	alpha = 0.5;
	%	meanmethod = 'riemann';
	%	distmethod = 'wasserstein';
            else
                alpha = 0.5;
                meanmethod = methodMean{method};
                distmethod = methodDist{method};
            end
            t = cputime;
            [Ytest d C] = mdm_alpha(COVtest,COVtrain,Ytrain, meanmethod, distmethod, alpha);
            
            
            det_mean(sub-5, testSession, method) = det(mean_covariances_alpha(COVtrain,meanmethod,alpha));
            e(sub-5, testSession, method) = cputime-t;
            
            if e(sub-5, testSession, method) == 0
                e(sub-5, testSession, method) = 0.0001;
            end
            ac(sub-5, testSession, method) = sum((labels-Ytest)==0)/(trialPerSession*4);
        end
        %Priority(0); %Reset priority to default
    end
end
g_means = squeeze(det_mean(:,:,2));
a_means = squeeze(det_mean(:,:,1));
wiener = g_means./a_means;
wiener_mean = mean(wiener,2);
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
    resMeanTmp(2) = sum(resMatrix(:,2,method));
    resMean(method,:) = resMeanTmp(2:end);
end
impr = ((resMatrix(:,3,2) - resMatrix(:,3,1))./resMatrix(:,3,1))*100;
[-10*log10(wiener_2)' impr]
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
