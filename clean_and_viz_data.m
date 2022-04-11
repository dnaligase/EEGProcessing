% initialize EEGLAB variables
if ~exist("ALLEEG", "var")
    eeglab
end

load ch_names.mat
fnames = dir("**/*_eyesclosed_*.hdf5");
ch_id = listdlg('PromptString',{'Select a channel to plot.',...
    'Only one can be selected at a time.',''},...
    'SelectionMode','single','ListString',channels);

to_avg = 128;
ALLEEG = [];
for i = length(fnames)

    path = fullfile(fnames(i).folder, fnames(i).name);
    patientName = split(fnames(i).folder, filesep);
    datastruct = ghdf5read(path);
    
    fs = datastruct.RawData.AcquisitionTaskDescription.SamplingFrequency;
    data = datastruct.RawData.Samples(1:64, :)';
    
    data = hifi(data, 1e6/fs, 0.5);
    data = lofi(data, 1e6/fs, 40);

    data = data';

    EEG = pop_importdata('dataformat','array','nbchan',64,'subject',patientName(end),'data','data', ...
        'srate',250,'pnts',0,'xmin',0, ...
        'chanlocs','/Volumes/EXTSTORAGE/Labrotation/64_channel.ced');
    EEG = eeg_checkset( EEG );
    EEG = pop_clean_rawdata(EEG, 'FlatlineCriterion','off', ...
        'ChannelCriterion','off','LineNoiseCriterion','off', ...
        'Highpass','off','BurstCriterion',20,'WindowCriterion','off', ...
        'BurstRejection','off','Distance','Euclidian');
    EEG = eeg_checkset( EEG );
    data_cleaned = EEG.data';
    data = data';
    
    figure("Name", "Patient #" + patientName(end) + ", " + channels(ch_id));

    EEG = pop_reref( EEG, []);
    [ALLEEG, EEG, CURRENTSET] = pop_newset(ALLEEG, EEG, 2,'gui','off');
    
    subplot(3, 2, 5)
    [p_non, f_non] = pwelch(data(:, :), [], [], 256, fs);
    topoplot(sum_freq_band(p_non, f_non, 8, 13), EEG.chanlocs)
    title("8-13 Hz")
    xlabel([to_avg, " points averaged per channel"])

    subplot(3, 2, 6)
    [p, f] = pwelch(data_cleaned(:, :), [], [], 256, fs);
    topoplot(sum_freq_band(p, f, 8, 13), EEG.chanlocs)

    title("8-13 Hz [ASR-Filtered]")
    xlabel([to_avg, " points averaged per channel"])

    [psd, freq] = pwelch(data(:, ch_id), [], [], 256, fs, "onesided");
    [psd_clean, freq_clean] = pwelch(data_cleaned(:, ch_id), [], [], 256, fs, "onesided");
    
%     LoHi-pass filtered
    subplot(3, 2, 2) 
    plot(freq, psd, "LineWidth", 2)
    hold on 
%     ASR-filtered
    plot(freq_clean, psd_clean, "LineWidth", 2)
    title("Power Spectral Density")
    xlabel("Frequency, Hz")
    legend("LoHi-pass filtered", "ASR-filtered")
    xlim([0, 50])
%     ASR-filtered Spectrogram
    subplot(3, 2, 4)
    [psds, freqs] = calculate_PSD_shift(data_cleaned(:, ch_id), 250, 512, 4, 2, 0);
    plot_single_DSA( ...
        psds, ...
        freqs, ...
        [1:size(psds, 2)] ...
        );

    subplot(3, 2, 1)
    plot(data(:, ch_id))
    hold on
    plot(data_cleaned(:, ch_id))
    title("EEG Trace")
    xlabel("Time, datapoints")
    legend("Raw signal", "ASR-filtered signal")

    subplot(3, 2, 3)
    [psds, freqs] = calculate_PSD_shift(data(:, ch_id), 250, 512, 4, 2, 0);
    plot_single_DSA( ...
        psds, ...
        freqs, ...
        [1:size(psds, 2)] ...
        );

%     axes_blueprint(patient, channels, 40, 60) % to update accordingly
        

end
