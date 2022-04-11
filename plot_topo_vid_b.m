function plot_topo_vid_b(fname, matrix, chanlocs, ...
    time_window, noverlap, power_lim, fs, fps, plot_grid, options)
% PLOT_TOPO_VID  Create an animated topoplot through time.
% Usage:
%       >> plot_topo_vid([], EEG, 250, 5, 3, 10); % create a video using
%                                                   overlapping windows

% Required Inputs:
%       fname – a name of video file to save as
%      matrix – a matrix of powers of size of size [time_windows, chs]
%    chanlocs – EEGLAB struct
% time_window – window to slide signal with, s
%    noverlap – overlap between consecutive windows, s
%   power_lim – limit for colorbar display, dB
%          fs – sampling frequency
%         fps – frame rate of output video
%   plot_grid – whether to plot raw eeg topomap (single patient only)
%     options – specify patient instance of BaseRaw class, optional

% Output:
%   .avi video file


arguments
    fname
    matrix
    chanlocs
    time_window
    noverlap
    power_lim
    fs = 250
    fps = 10
    plot_grid logical = false
    options.patient
    
end

% count number of frames for future video given time_window and noverlap
framesNo = size(matrix, 1);

% allocate memory for future figures
F = struct('cdata', cell(1,framesNo), 'colormap', cell(1,framesNo));

if plot_grid
    F_gridplot = struct('cdata', cell(1,framesNo), 'colormap', cell(1,framesNo));
end

pb = CmdLineProgressBar(("Preparing "+framesNo+" frames using "+noverlap+"s overlap and "+time_window+"s window size... "));
for i = 1:framesNo
    pb.print(i,framesNo)

    figure('Visible','off');
    topoplot(matrix(i, :), chanlocs, 'electrodes', 'labels');
    
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
    if plot_grid
        axes_blueprint(options.patient, channels, tmin, tmax, "off")
        F_gridplot(i) = getframe(gcf);
    end
end

writerObj = VideoWriter(fname);
writerObj.FrameRate = fps;
disp("Creating " + length(F) / writerObj.FrameRate + " s video...");
% open the video writer
open(writerObj);
% write the frames to the video
for i=1:length(F)
    % convert the image to a frame
    frame_topo = F(i);
    
    if plot_grid
        frame_grid = F_gridplot(i); 
        cdata = cat(2, frame_grid.cdata, frame_topo.cdata);
        frame = struct('cdata', cdata, 'colormap', []);
    else
        frame = frame_topo;
    end

    writeVideo(writerObj, frame);
end

% close the writer object
close(writerObj);
disp('Done.')
disp(strcat('Find it as ', pwd, filesep, fname))
end