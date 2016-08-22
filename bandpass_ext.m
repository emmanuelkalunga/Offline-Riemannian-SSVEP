function x = bandpass_ext(f1, f2, f3, S, H, ftype)
Fs = H.SampleRate;
Fn = Fs/2;
if nargin == 3
    n_butter = 8;
    [b, a] = butter(n_butter, f1./Fn, 'bandpass');
    x = filtfilt(b, a, S); %Band pass
elseif nargin >= 5
    Rp = 3; Rs = 10;
    Wp1 = f1/Fn; Ws1 = [f1(1)-1 f1(2)+1]/Fn;
    Wp2 = f2/Fn; Ws2 = [f2(1)-1 f2(2)+1]/Fn;
    Wp3 = f3/Fn; Ws3 = [f3(1)-1 f3(2)+1]/Fn;
    
    if ftype == 2 %IIR
        [n1,Wn1] = buttord(Wp1,Ws1,Rp,Rs);
        [n2,Wn2] = buttord(Wp2,Ws2,Rp,Rs);
        [n3,Wn3] = buttord(Wp3,Ws3,Rp,Rs);

        [b1, a1] = butter(n1, Wn1, 'bandpass');
        [b2, a2] = butter(n2, Wn2, 'bandpass');
        [b3, a3] = butter(n3, Wn3, 'bandpass');

        x1 = filtfilt(b1, a1, S); %Band pass
        x2 = filtfilt(b2, a2, S); %Band pass
        x3 = filtfilt(b3, a3, S); %Band pass
    else%FIR
        mags = [0 1 0];
        devs = [0.01 0.05 0.01];
        
        fcuts1 = [f1(1)-0.5 f1(1) f1(2) f1(2)+0.5];
        [n1,Wn1,beta1,ftype1] = kaiserord(fcuts1,mags,devs,Fs);
        n1 = n1 + rem(n1,2);
        b1 = fir1(n1,Wn1,ftype1,kaiser(n1+1,beta1),'noscale');
        %-------------------------------
        fcuts2 = [f2(1)-0.5 f2(1) f2(2) f2(2)+0.5];
        [n2,Wn2,beta2,ftype2] = kaiserord(fcuts2,mags,devs,Fs);
        n2 = n2 + rem(n2,2);
        b2 = fir1(n2,Wn2,ftype2,kaiser(n2+1,beta1),'noscale');
        %--------------------------------
        fcuts3 = [f3(1)-0.5 f3(1) f3(2) f3(2)+0.5];
        [n3,Wn3,beta3,ftype3] = kaiserord(fcuts3,mags,devs,Fs);
        n3 = n3 + rem(n3,2);
        b3 = fir1(n3,Wn3,ftype3,kaiser(n3+1,beta3),'noscale');
        x1 = filtfilt(b1, 1, S); %Band pass
        x2 = filtfilt(b2, 1, S); %Band pass
        x3 = filtfilt(b3, 1, S); %Band pass
            
    end
    x = [x1 x2 x3];
end


% function X = bandpass_ext(S, H, Hd)
%     for k = 1:3
%         x{k} = filter(Hd{1},S);
%     end
%     X = [x{1} x{2} x{3}];
% end