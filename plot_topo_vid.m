function plot_topo_vid(fname, patient, fs, time_window, noverlap, fps, l_freq, h_freq)
% PLOT_TOPO_VID  Create an animated topoplot through time.
% Usage:
%       >> plot_topo_vid([], patient, 250, 5, 3, 10); % create a video using
%                                                   overlapping windows

% Required Inputs:
%       fname – a name of video file to save as
%         patient – a struct from patientLAB (patient.data of size [chs, datapoints])
%          fs – sampling frequency
% time_window – window to slide signal with, s
%    noverlap – overlap between consecutive windows, s
%         fps – frame rate of output video
%      l_freq – lower frequency boundary
%      h_freq – upper frequency boundary

% Output:
%   .avi video file


arguments
    fname
    patient
    fs = patient.fs
    time_window = 5
    noverlap = 3
    fps double {mustBeNumeric} = 10
    l_freq = 8
    h_freq = 13
end

if isempty(fname)
    fname = strcat(patient.subject, '_',  ...
        string(datetime("now", "Format", 'yyyy-MM-dd')), '.avi');
    warning("No filename specified. Inferring filename from data.")
end

assert(fs*noverlap <= fs*time_window - 1, ...
    "noverlap must not exceed time_window! " + ...
    "Expected noverlap less than " + time_window + ...
    "s, received " + noverlap + "s.")

% count number of frames for future video given time_window and noverlap
framesNo = fix((size(patient.data,2) - fs*time_window) / (fs*time_window - fs*noverlap)) + 1;

% allocate memory for future figures
F = struct('cdata', cell(1,framesNo), 'colormap', cell(1,framesNo));
F_gridplot = struct('cdata', cell(1,framesNo), 'colormap', cell(1,framesNo));

[psd_whole,freq_whole] = pwelch(patient.data',[],[],256,fs);
sum_powers = sum_freq_band(psd_whole, freq_whole, l_freq, h_freq);
power_lim = max(abs(sum_powers)); 

channels = {patient.chanlocs(:).labels};

pb = CmdLineProgressBar(["Getting "+framesNo+" frames using "+noverlap+"s overlap and "+time_window+"s window size... "]);
for i = 1:framesNo
    pb.print(i,framesNo)

    patient.data = patient.data';
    [psd,freq] = pwelch(patient.data((time_window*fs - noverlap*fs) * (i-1) + 1:time_window*fs + (time_window*fs - noverlap*fs) * (i-1) + 1, :), [],[],256,fs);
    patient.data = patient.data';

    figure('Visible','off');
    topoplot(sum_freq_band(psd, freq, l_freq, h_freq), patient.chanlocs, 'electrodes', 'labels');
    
    title([num2str(time_window + (time_window - noverlap) * (i-1)), ' s'])
    
    % display ch names exceeding the std of summed powers
%     band_sum = sum_freq_band(psd, freq);
%     idxs = band_sum > std(band_sum);
%     chs = channels(idxs);
%     annotation('textbox',[0.05 0 1.0000 .95], ...
%     'String',chs', EdgeColor='none');

    % colorbar settings
    c = colorbar;
    c.Label.String = 'dB';
    caxis([-1 1] * power_lim)

    F(i) = getframe(gcf);
    
    tmin = floor(((time_window*fs - noverlap*fs) * (i-1) + 1)/fs);
    tmax = floor((time_window*fs + (time_window*fs - noverlap*fs) * (i-1) + 1)/fs);
    axes_blueprint(patient, channels, tmin, tmax, "off")
    F_gridplot(i) = getframe(gcf);
end

writerObj = VideoWriter(fname);
writerObj.FrameRate = fps;
disp("Creating " + length(F_gridplot) / writerObj.FrameRate + " s video...");
% open the video writer
open(writerObj);
% write the frames to the video
for i=1:length(F)
    % convert the image to a frame
    frame_grid = F_gridplot(i); 
    frame_topo = F(i);
    cdata = cat(2, frame_grid.cdata, frame_topo.cdata);
    frame = struct('cdata', cdata, 'colormap', []);
    writeVideo(writerObj, frame);
end

% close the writer object
close(writerObj);
disp('Done.')
disp(strcat('Find it as ', pwd, filesep, fname))
end