for sub = 6:17
    %% Load data
    disp(['Load data subject ', num2str(sub)]);
    [S_all, H_all] = loaddata(sub);%Returns cells of data from all available sessions
    nbrSessions = length(S_all);
    %S = cat(3,S_all{:});
    %H = H_all{1};
    subject = sub-5;
    for sess = 1:nbrSessions
        S = S_all{sess};
        H = H_all{sess};
        if subject > 9
           fn = ['./extracted_data/S',num2str(subject),'_',num2str(sess),'.mat'];
        else
           fn = ['./extracted_data/S0',num2str(subject),'_',num2str(sess),'.mat'];
        end
        labels = H.Label;
        dimensions = H.PhysDim;
        sample_rate = H.AS.SampleRate;
        physical_min = H.PhysMin;
        physical_max = H.PhysMax;
        digital_min = H.DigMin;
        digital_max = H.DigMax;
        transducer = H.Transducer;
        prefilter = H.PreFilt;
        threshold = H.THRESHOLD;
        event_duration = H.EVENT.DUR;
        event_position = H.EVENT.POS;
        event_type = H.EVENT.TYP;
        save(fn,'S','labels', 'dimensions', 'sample_rate', 'physical_min', 'physical_max', 'digital_min', 'digital_max', 'threshold', 'event_duration','event_position', 'event_type');
    end
end
