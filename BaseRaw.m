classdef BaseRaw < handle
    properties
        data
        subject
        times
        psd
        freq
        raw
        no_chan
        DSA
    end
    properties (SetAccess = private)
        fs
        chanlocs
        first_samp = 1
        last_samp
    end
    properties (Hidden = true)
        time_as_index
        id = 1
    end
%% Methods for Basic EEG Processing and Visualization
    methods
        function obj = BaseRaw(EEG, picks)
            arguments
                EEG
                picks = (1: size(EEG.data, 1))
            end

            if nargin >= 1
                obj.data = EEG.data(picks, :);
                obj.raw = EEG.data(picks, :);
                obj.fs = EEG.srate;
                obj.no_chan = size(obj.data,1);

                if  ~isfield(EEG, 'subject') || isempty(EEG.subject)
                    obj.subject = string(obj.id);
                    increment_id = @(x) x + 1;
                    obj.id = increment_id(obj.id);
                else
                    obj.subject = EEG.subject;
                end

                if  ~isfield(EEG, 'times') || isempty(EEG.times)
                    obj.times = obj.create_times();
                else
                    % convert EEGLAB times to ms
                    obj.times = EEG.times / 1000;
                end

                [obj.psd, obj.freq] = pwelch(obj.data',[],[],256,obj.fs);
                if isfield(EEG, 'chanlocs')
                    obj.chanlocs = EEG.chanlocs(1,picks);
                else
                    obj.chanlocs = [];
                end
                obj.last_samp = length(obj);
                obj.time_as_index = (1 : length(obj));
            end

        end

        function filterEEG(obj,hi,lo)
            % Apply Basic High and Low Pass filter
            if ~isempty(hi) && ~isempty(lo)
                obj.data = hifi(obj.raw', 1e6/obj.fs, hi)';
                obj.data = lofi(obj.data', 1e6/obj.fs, lo)';
            elseif ~isempty(lo)
                obj.data = lofi(obj.raw', 1e6/obj.fs, lo)';
            elseif ~isempty(hi)
                obj.data = hifi(obj.raw', 1e6/obj.fs, hi)';
            end
        end

        function crop(obj, tmin, tmax)
            % Crop at specific timepoint 
            % To Do: Crop at index, not at second
            arguments
                obj
                tmin {mustBeGreaterThanOrEqual(tmin, 0)}
                tmax {mustBeReal} = obj.times(end)
            end

            start = tmin * obj.fs + 1;
            stop = tmax * obj.fs + 1;

            obj.first_samp = obj.time_as_index(1) + start - 1;
            obj.last_samp = obj.time_as_index(1) + stop - 1;

            obj.time_as_index = obj.time_as_index(1, ...
                start:stop);

            obj.data = obj.data(:, start:stop);
            times_vec = obj.times(:, start:stop);
            obj.times = times_vec - min(times_vec);

        end

        function [r, r1, meta] = windowedPower(obj,  time_window, noverlap, ...
                lfreq, hfreq, verbose)
            % calculate powers in a windowed fashion
            arguments
                obj
                time_window (1,1) double
                noverlap (1,1) double
                lfreq (1,1) = 8
                hfreq (1,1) = 13
                verbose logical = false
            end

            rowsNo = fix((size(obj.data, 2) - obj.fs*time_window) / ...
                (obj.fs*time_window - obj.fs*noverlap)) + 1;
            matrix = zeros(rowsNo, obj.no_chan);
            matrix_DSA = zeros(length(obj.freq), obj.no_chan,rowsNo);
            for j = 1:rowsNo
                begin = (time_window*obj.fs - noverlap*obj.fs) * (j-1) + 1;
                stop = time_window*obj.fs + ...
                    (time_window*obj.fs - noverlap*obj.fs) * (j-1) + 1;

                [powers,cur_DSA] = obj.sum_power_segment((begin:stop), lfreq, hfreq);
                matrix(j, :) = powers;
                matrix_DSA(:,:,j) = cur_DSA;
            end

            r = matrix;
            r1 = matrix_DSA;
            meta = struct('time_window', time_window, 'noverlap', noverlap);

            if verbose
                disp(meta);
            end
        end

        function r = sum_power(obj, l_freq, h_freq)
            r = obj.sum_freq_band(obj.psd, l_freq, h_freq);
        end

        function [r,r1] = sum_power_segment(obj, segment, l_freq, h_freq)
            transposed = obj.data';
            % segment in datapoints as ``double``
            art_mat = abs(transposed(segment,:)) > 100;
            art_vec = ~any(art_mat);
            psd_window = NaN(129,obj.no_chan);
            freq_window = NaN(129,1);
            r = NaN(1,obj.no_chan);
            if sum(art_vec) > 0
                [psd_window(:,art_vec), freq_window] = pwelch(transposed(segment, art_vec), ...
                    [],[],256,obj.fs);
                r = obj.sum_freq_band(psd_window, l_freq, h_freq);
            end
            r1 = psd_window;
        end

   %Test functions
        function matrix_DSA = create_DSA(obj,time_window, noverlap)
            %Creates DSA of Power Spectrum for window and noverlap
            rowsNo = fix((size(obj.data, 2) - obj.fs*time_window) / ...
                (obj.fs*noverlap)) + 1;
            matrix_DSA = zeros(length(obj.freq), obj.no_chan,rowsNo);
            cnt = 1;
            for j=1:obj.fs*noverlap:length(obj.data)-time_window*obj.fs
                begin = j;                              
                stop = j+time_window*obj.fs-1;
                [powers] = obj.power_segment((begin:stop));
                matrix_DSA(:,:,cnt) = powers;
                cnt = cnt+1;
            end             
        end

        function [r, meta] = windowedPower1(obj,  time_window, noverlap, ...
                lfreq, hfreq, verbose)
            % calculate powers in a windowed fashion
            arguments
                obj
                time_window (1,1) double
                noverlap (1,1) double
                lfreq (1,1) = 8
                hfreq (1,1) = 13
                verbose logical = false
            end

            matrix_DSA = obj.create_DSA(time_window, noverlap);            
            r = obj.sum_freq_band(matrix_DSA, lfreq, hfreq);
            meta = struct('time_window', time_window, 'noverlap', noverlap);

            if verbose
                disp(meta);
            end
        end

        function [psd_window] = power_segment(obj, segment)
            transposed = obj.data';
            % segment in datapoints as ``double``
            art_mat = abs(transposed(segment,:)) > 100;
            art_vec = ~any(art_mat);
            psd_window = NaN(129,obj.no_chan);
            if sum(art_vec) > 0
                [psd_window(:,art_vec), ~] = pwelch(transposed(segment, art_vec), ...
                    [],[],256,obj.fs);
            end
        end

        function datavec = sum_freq_band(obj, psd, l_freq, h_freq)
            assert(size(psd, 1) == size(obj.freq, 1), "psd and freqs dims don't match!")

            idxs = obj.freq >= l_freq & obj.freq < h_freq;
            if size(psd,3) == 1
                datavec = sum(psd(idxs, :), 1);
            elseif size(psd,3) > 1                  % Works on 3D matrix
                datavec = sum(psd(idxs, :,:), 1);
                datavec = squeeze(datavec)';
            end
        end

        function plt = plot_single_channel(obj,ch,time_window, noverlap,xax)
            % Plot single channel DSA and Raw EEG 
            % To Do: 
            % check time_vector and position of the colorbar
            plt = figure('WindowState','maximized','Name',obj.subject);
            ch_name = obj.chanlocs(ch).labels;
            if isempty(xax)
                xax = 0:noverlap:size(obj.data,2)/obj.fs-time_window;
            end
            sgtitle(ch_name,'FontName','Arial','FontSize',12,'FontWeight','Bold')
            % Calculate DSA
            matrix_DSA = create_DSA(obj,time_window, noverlap);
            subplot(2,1,1)
            pcolor(xax,obj.freq,10*log10(squeeze(matrix_DSA(:,ch,:))))
            shading flat
            shading interp
            colormap turbo
            colorbar
            caxis([-20 20])
            ylim([0 47])
            xlim([xax(1),xax(end)])
            ylabel('Frequency [Hz]')
            set(gca,'FontName','Arial','FontSize',12,'FontWeight','Bold', 'LineWidth', 1)
            box(gca,'off')
            subplot(2,1,2)
            plot(obj.times,obj.data(ch,:))
            ylabel('Amplitude [mV]')
            xlabel('Time [s]')
            set(gca,'FontName','Arial','FontSize',12,'FontWeight','Bold', 'LineWidth', 1)
            box(gca,'off')
        end

        function tg = plot_all_channels(obj,matrix_DSA,time_window, noverlap,xax)
            % Plot channels DSA and Raw EEG 
            % To Do: 
            % Display Channel name, check time_vector and position of
            % colorbar
            if isempty(xax)
                xax = 0:noverlap:size(obj.data,2)/obj.fs-time_window;
            end
            plt = figure('WindowState','maximized','Name',obj.subject);
            % Calculate DSA
            if isempty(matrix_DSA)
            matrix_DSA = create_DSA(obj,time_window, noverlap);
            end
            tg = uitabgroup(plt); % tabgroup
            for i = 1:obj.no_chan
                ch_name = obj.chanlocs(i).labels;
                thistab = uitab(tg,"Title",ch_name); % build iith tab
                axes('Parent',thistab); % somewhere to plot
                sgtitle(ch_name,'FontName','Arial','FontSize',12,'FontWeight','Bold')
                subplot(2,1,1)
                pcolor(xax,obj.freq,10*log10(squeeze(matrix_DSA(:,i,:))))
                shading flat
                shading interp
                colormap turbo
                colorbar
                caxis([-20 20])
                ylim([0 47])
                xlim([xax(1),xax(end)])
                ylabel('Frequency [Hz]')
                set(gca,'FontName','Arial','FontSize',12,'FontWeight','Bold', 'LineWidth', 1)
                box(gca,'off')
                subplot(2,1,2)
                plot(obj.times,obj.data(i,:))
                ylabel('Amplitude [mV]')
                xlabel('Time [s]')
                set(gca,'FontName','Arial','FontSize',12,'FontWeight','Bold', 'LineWidth', 1)
                box(gca,'off')
            end
        end

        function tg = plot_power_bands(obj,matrix_DSA,time_window,noverlap,xax)
            % User can supply DSA_matrix (e.g. Averaged DSA over all
            % patients), if nothing is supplied function operates on single
            % patient basis
            % X-Axis can be supplied,otherwise xax will be calculated based
            % on window and noverlap
            if isempty(xax)
                xax = 0:noverlap:size(obj.data,2)/obj.fs-time_window;
            end
            if isempty(matrix_DSA)
            matrix_DSA = create_DSA(obj,time_window, noverlap);
            end
            plt = figure('WindowState','maximized','Name',obj.subject);
            l = obj.sum_freq_band(matrix_DSA, 1, 6.99);
            a = obj.sum_freq_band(matrix_DSA, 7, 12.99);
            b = obj.sum_freq_band(matrix_DSA, 13, 29.99);
            g = obj.sum_freq_band(matrix_DSA, 30, 46);
            % Create 4 plots with bandpower over time
            tg = uitabgroup(plt); % tabgroup
            sz = 5;
            for i = 1:obj.no_chan
                ch_name = obj.chanlocs(i).labels;
                thistab = uitab(tg,"Title",ch_name); % build iith tab
                axes('Parent',thistab); % somewhere to plot
                sgtitle(ch_name,'FontName','Arial','FontSize',12,'FontWeight','Bold')
                ax(1) = subplot(2,2,1);
                scatter(xax,l(:,i),sz,'filled')
                ylabel('Power')
                title('Low Frequency Band')
                ax(2) = subplot(2,2,2);
                scatter(xax,a(:,i),sz,'filled')
                ylabel('Power')
                title('Alpha Band')
                ax(3) = subplot(2,2,3);
                scatter(xax,b(:,i),sz,'filled')
                ylabel('Power')
                title('Beta Band')
                ax(4) = subplot(2,2,4);
                scatter(xax,g(:,i),sz,'filled')
                ylabel('Power')
                title('Gamma Band')
                set(ax,'FontName','Arial','FontSize',12,'FontWeight','Bold', 'LineWidth', 1)
                box(ax,'off')
                xlim(ax,[xax(1),xax(end)])
            end            
        end

        function fig = plot_artefact_matrix(obj,time_window, threshold,time_vec)

            size_data = size(obj.data,2);
            no_epochs = floor(size_data(1)/(time_window*obj.fs));
            M_woA = zeros(obj.no_chan, no_epochs); % Matrix to contain time series without artefacts
            if isempty(time_vec)
                time_vec = 1:time_window:no_epochs*time_window;
            end
            for ch = 1:obj.no_chan   % loop over (the 64) channels

                channel_data = obj.data(ch, :); % data set for channel under consideration

                for j = 1:no_epochs % loop over (the 4s) epochs

                    first_sample = (j-1)*time_window*obj.fs + 1;
                    last_sample = j*time_window*obj.fs;
                    epoch_data = abs(channel_data(first_sample:last_sample));
                    artifacts = epoch_data >= threshold;
                    if (any(artifacts)) % set values for Matrix M = 1 if artifacts found
                        M_woA(ch, j) = 1;
                    end
                end
            end

            % Plot matrix M
            figure('WindowState','maximized','Name',obj.subject)
            [r, c] = size(M_woA);                               % Get the matrix size
            imagesc((1:c)+0.5, (1:r)+0.5, M_woA);               % Plot the image
            hold on;
            colormap(gray); % Use a gray colormap
            xlabel('Number of Epoch')
            axis equal
            for i = 1:obj.no_chan
                locs(i) = {obj.chanlocs(i).labels};
            end
            %y_labels = {channels{1}:channels{r}};
            set(gca,'XTick', 1:2:c, 'YTick', 1:r, ...        % Change axes properties
                'YTickLabel', locs,...
                'XLim', [1 c+1], 'YLim', [1 r+1], ...
                'GridLineStyle', '-', 'XGrid', 'on', 'YGrid', 'on','FontSize',6);
            hold off
        end

        % Test over
        function r = apply_function(obj, func, channel_wise, inplace)
            arguments
                obj
                func
                channel_wise logical = false
                inplace logical = true
            end

            if inplace
                obj.data = obj.apply_for_signal(func, channel_wise);
            else
                r = obj.apply_for_signal(func, channel_wise);
            end
        end

        function r = plot_channels_3d(obj, return_figure)
            arguments
                obj
                return_figure logical = false
            end
            assert(~isempty(obj.chanlocs), "Channel locations not provided.")
            figure;
            plot3([obj.chanlocs.X],[obj.chanlocs.Y],[obj.chanlocs.Z], 'o');
            for i = 1:obj.no_chan
                locs(i) = {obj.chanlocs(i).labels};
            end
            text([obj.chanlocs.X],[obj.chanlocs.Y],[obj.chanlocs.Z],locs,'VerticalAlignment','bottom','HorizontalAlignment','right','FontSize',6)
            title(('3D channel location for #' + string(obj.subject)))
            xlabel('X')
            ylabel('Y')
            zlabel('Z')
            if return_figure; r = gcf(); end
        end
        function r = length(obj)
            r = length(obj.data);
        end
    end

    methods (Access = protected)
        function r = create_times(obj)
            time_step = 1/obj.fs;
            endpoint = length(obj.data)/obj.fs;
            r = [time_step:time_step:endpoint];
        end

        function r = apply_for_signal(obj, func, channel_wise)

            if channel_wise
                n = 1;
                template = zeros(size(obj.data));
                while n ~= size(obj.data, 1)
                    template(n, :) = func(obj.data(n, :));
                    n = n + 1;
                end
                r = template;
            else
                r = func(obj.data);
            end
        end

    end
end