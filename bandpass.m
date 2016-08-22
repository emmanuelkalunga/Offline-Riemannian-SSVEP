function X = bandpass_ext(S, H, Hd)
    for k = 1:3
        x{k} = filter(Hd{1},S)
    end
    X = [x{1} x{2} x{3}];
end