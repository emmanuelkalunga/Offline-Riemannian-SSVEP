%**************************************************************************
%File: swelling effect
%Description: analyse swelling effect in different averaging,
%Author: Emmanuel K. Kalunga
%**************************************************************************
%methodMean = {'riemann','arithmetic','logeuclid','ld'};
%methodDist = {'riemann','euclid','logeuclid','ld'};

% methodMean = {'arithmetic','logeuclid','ld','riemann','ld'};
% methodDist = {'euclid','logeuclid','ld','riemann','ld'};
methodMean = {'arithmetic','harmonic','geometric','riemann','logeuclid','sdivergence','ld','bhat','opttransp'};
methodDist = {'euclid','euclid','euclid','riemann','logeuclid','sdivergence','ld','bhat','opttransp'};

tLen = 4;
delay = 2;
ac = zeros(12, 5, length(methodMean));
e = zeros(12, 5, length(methodMean));
clear Cm detMean trMean
for sub = 11+6 %6:17
    clear x_all H_all P X Pm PSD
    %% Load data
    disp('********************************************************');
    disp(['Load data subject ', num2str(sub)]);
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
    clear P
    for k = 1:Nt %loop for evrey trial
        for cl = 1:4
            P{cl}(:,:,k) = shcovft((X{cl}(:,:,k))'); % J. Schaefer Shrinkage covariance from Barachant toolbox
        end
    end
    %% Get determinant of all trials per class for each subject and the determinant of means computed with defferent means. And plot the averages covmat
    % THIS IS DONE ONLY FOR ANALYSIS PURPURSES (CAN BE TAKEN TO A DIFFERENT
    % FILE
    %-- Trials determinants and traces
    determinants = [];
    traces = [];
    for k = 1:Nt %loop for evrey trial
        for cl = 1:4
            determinants(k,cl,sub-5) = det(P{cl}(:,:,k));
            traces(k,cl,sub-5) = trace(P{cl}(:,:,k));
        end
    end
    %-- Means determinants and traces
    COVMAT = cat(3, P{:});
    LABELS = [zeros(1,Nt) ones(1,Nt) 3*ones(1,Nt) 2*ones(1,Nt)];
    labels = unique(LABELS);
    for method = 1:length(methodMean)
        
        
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
        for cl = 1:4
            Cm{sub-5}(:,:,cl,method) = mean_covariances_alpha(COVMAT(:,:,LABELS==labels(cl)),meanmethod,alpha); %alpha = 0.6
            detMean(cl,method,sub-5) = det(squeeze(Cm{sub-5}(:,:,cl,method)));
            trMean(cl,method,sub-5) = trace(squeeze(Cm{sub-5}(:,:,cl,method)));
        end              
%         if method == 3 %Do Bhattacharyya which is alpha-div with alpha = 0
%             for cl = 1:4
%                 Cm{sub-5}(:,:,cl,method) = mean_covariances_alpha(COVMAT(:,:,LABELS==labels(cl)),methodMean{method},0); %alpha = 0 for Bhattacharryya
%                 detMean(cl,method,sub-5) = det(squeeze(Cm{sub-5}(:,:,cl,method)));
%             end
%         else
%             for cl = 1:4
%                 Cm{sub-5}(:,:,cl,method) = mean_covariances_alpha(COVMAT(:,:,LABELS==labels(cl)),methodMean{method},0.6); %alpha = 0.6
%                 detMean(cl,method,sub-5) = det(squeeze(Cm{sub-5}(:,:,cl,method)));
%             end
%         end
    end
    %-- Plot mean cov
    %if (sub-5) == 10
    detm = detMean(:,:,sub-5)';
    trm = trMean(:,:,sub-5)';
    %-- plot log det value of class 2 trials (subject 10)
    figure, plot(log(determinants(:,2,sub-5)),':+','LineWidth', 3,'DisplayName', 'trial det')
    legend('-DynamicLegend');
    hold all;
    
    %-- plot the means obtained with different techniques. 
    plot(repmat(log(detm(1,2)),1,Nt),'-r', 'LineWidth', 2,'DisplayName', 'Arithmetic')
    plot(repmat(log(detm(2,2)),1,Nt),'--m', 'LineWidth', 2,'DisplayName','Harmonic')
    %plot(repmat(log(detm(3,2)),1,Nt),'--k', 'LineWidth', 2,'DisplayName','Geometric')
    %--------------------------------------------   
    plot(repmat(log(detm(4,2)),1,Nt),'-c', 'LineWidth', 2,'DisplayName','AIR')
    plot(repmat(log(detm(5,2)),1,Nt),'-.k', 'LineWidth', 2,'DisplayName','Log-Euclidean')
    plot(repmat(log(detm(6,2)),1,Nt),'-.m', 'LineWidth', 2,'DisplayName','S-divergence')
    plot(repmat(log(detm(7,2)),1,Nt),'--g', 'LineWidth', 2,'DisplayName','log-det \alpha-div')% or methodMean{4} (ld)
    plot(repmat(log(detm(8,2)),1,Nt),'--b', 'LineWidth', 2,'DisplayName','Bhattacharyya')
    plot(repmat(log(detm(9,2)),1,Nt),'--y', 'LineWidth', 2,'DisplayName','Wasserstein')
     
    set(gca,'FontSize',14,'fontWeight','normal')
    set(findall(gcf,'type','text'),'FontSize',14,'fontWeight','normal')
    xlabel('Trials')
    ylabel('Log-det values of covariance matrix')
    spaceplots;
    ylabels = {{'Resting class'}; {'13 Hz class'};{'21 Hz class'};{'17 Hz class'}};
    
    yticklabel = cell(24,1);
    yticklabel_tmp = cat(2,repmat(H_all{1}.Label,3,1), num2cell(repmat('-',24,1),2), reshape(repmat(mat2cell(['13Hz';'21Hz';'17Hz'],[1 1 1],[4]),1,8)',24,1) );
    for k = 1:size(yticklabel_tmp,1)
        yticklabel{k} = cat(2,yticklabel_tmp{k,:});
    end  
    for method = 1:4%1:length(methodMean)
        tmp = (Cm{sub-5}(:,:,:,method));
        clims = [min(tmp(:)) max(tmp(:))];
        figure
        for cl = 1:4
            subplot(2,2,cl)
            imagesc(squeeze(Cm{sub-5}(:,:,cl,method)), clims);colormap(flipud(hot));%colorbar;
            %title(['method ' methodMean{method} ', class' num2str(cl)])
            title(ylabels{cl})
            set(gca,'FontSize',14,'fontWeight','normal')
            set(findall(gcf,'type','text'),'FontSize',14,'fontWeight','normal')
            set(gca, 'XTick', []);
            %if cl == 1
                set(gca, 'YTick', 1:24);
                set(gca, 'YTickLabel', yticklabel);
            %else
            %    set(gca, 'YTick', []);
            %end
        end
        %spaceplots;
    end
    %end
end
