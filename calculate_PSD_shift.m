function [PSD freq]= calculate_PSD_shift (EEG, fs, NFFT, windowlength, shift,rel)

% calculate_PSD calculates the PSD of the input signal
% EEG:          input signal
% NFFT:         Unused (Standard values from matlab will be used)
% fs:           sample rate
% windowlength: window in s the power is calculated over
% shift:        shift in s:
% PSD:          power spectral density matrix
% frequency_resolution: freq info
%rel:           boolean for relative PSD

%h = waitbar(0,'Progress of PSD calculation');
%NFFT=128;
cnt =1;

for a=1:fs*shift:length(EEG)-windowlength*fs
        start = a;
        stop = a+windowlength*fs-1;
        section = EEG(start:stop,:);
        [BUF2, freq] = pwelch(section,[],[],NFFT,fs,'onesided');
        if rel == 1
        sum_BUF2 = sum(BUF2,1);
        BUF2 = BUF2./sum_BUF2;
        end
        
        PSD(:,cnt)=(BUF2);
        clear BUF BUF2
        cnt = cnt+1;
        
        
end     
