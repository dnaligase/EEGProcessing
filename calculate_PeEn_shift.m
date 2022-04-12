function [PeEn,Prob]= calculate_PeEn_shift (EEG, fs, windowlength, shift,art)

% calculate_PSD calculates the PSD of the input signal
% EEG:          input signal

% fhighpass:    highpass cutoff in Hz
% fs:           sample rate
% windowlength: window in s the power is calculated over
% shift:        shift in s:
% PSD:          power spectral density matrix
% frequency_resolution: freq info

%h = waitbar(0,'Progress of PSD calculation');
%NFFT=128;
cnt =1;
start=1:fs*shift:length(EEG);%-windowlength*fs;
stop=start+(windowlength*fs-1);
PeEn(1,1)= NaN;
Prob(1,1:6) = NaN;
for a=1:length(start)
    %         start = a;
    %         stop = a+windowlength*fs-1;
    if stop(a)<=length(EEG)
        segment_art = art(start(a):stop(a));
        if any(segment_art)
            PeEn(a,1)= NaN;
            Prob(a,1:6) = NaN; 
        else
            [PeEn(a,1),~,Prob(a,1:6)]=peen(EEG(start(a):stop(a)),3,1);
        end

    end
    %cnt = cnt+1;
end
