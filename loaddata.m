% function [S, H] = loaddata(subject, session)
% channels = 0;
% direct = ['../../Data/Online/subject',num2str(subject),'/Training'];
% fnames = dir(fullfile(direct));
% Filename = ['../../Data/Online/subject',num2str(subject),'/Training/',fnames(session+2).name];
% [S, H] = sload(Filename, channels, 'OVERFLOWDETECTION:OFF');


function [S, H] = loaddata(subject, session)
channels = 0;
direct = ['../../../Data/Online2/subject',num2str(subject),'/Training'];
fnames = dir(fullfile(direct));
nbrSessions = length(fnames) - 2;
if nargin == 1 %Load all session
    for session = 1:nbrSessions
        Filename = ['../../../Data/Online2/subject',num2str(subject),'/Training/',fnames(session+2).name];
        [S{session}, H{session}] = sload(Filename, channels, 'OVERFLOWDETECTION:OFF');
    end
else
    Filename = ['../../../Data/Online2/subject',num2str(subject),'/Training/',fnames(session+2).name];
    [S, H] = sload(Filename, channels, 'OVERFLOWDETECTION:OFF');
end