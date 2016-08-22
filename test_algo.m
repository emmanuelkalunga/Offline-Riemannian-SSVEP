%Test Mean and distance algorithms

methodMean = 'opttransp';
methodDist = 'opttransp';
tLen = 4;
delay = 2;
sub = 6

%Get Data
[S_all, H_all] = loaddata(sub);
Fs = H_all{1}.SampleRate;
nbrSessions = length(S_all);
sessions = 1:nbrSessions;
% Preprocessing of all available sessions (Same for training and test data)
disp(['Preprocessing subject ', num2str(sub)]);
for session = 1:nbrSessions
    x_all{session} = bandpass_filter_ext([12.95 13.05], [16.9 17.1], [20.9 21.1], S_all{session}, H_all{session});
end
X = get_trials(x_all, H_all, tLen, delay);
%Covariance matrices of all trials
disp(['Preprocessing subject ', num2str(sub)]);
Nt = size(X{1},3); %Number of trial
for k = 1:Nt %loop for evrey trial
    for cl = 1:4
        P{cl}(:,:,k) = shcovft((X{cl}(:,:,k))'); % J. Schaefer Shrinkage covariance from Barachant toolbox
    end
end
trainTrials = 9:16;
COVMAT = cat(3, P{1}(:,:,trainTrials), P{2}(:,:,trainTrials), P{3}(:,:,trainTrials), P{4}(:,:,trainTrials));
% %Save data for python comparison
% pyCOVMAT = permute(COVMAT, [3 1 2]);
% save('pyCOVMAT.mat', 'pyCOVMAT');

%Test Wasserstein mean
MW1 = opttransp_mean(COVMAT(:,:,1:8));
MW2 = opttransp_mean(COVMAT(:,:,9:16));
MW3 = opttransp_mean(COVMAT(:,:,17:24));
MW4 = opttransp_mean(COVMAT(:,:,25:32));

%Test S-Div 
for i = 1:size(COVMAT,3)
    sdiv(i) = distance_sdiv(COVMAT(:,:,20),COVMAT(:,:,i));
end

%Test S-Div mean
[MS,niter, crit, Sn] = sdiv_mean(COVMAT(:,:,1:8));