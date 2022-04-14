classdef BaseRaw < handle
    properties
        data
        subject
        times
        psd
        freq
        raw
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

    methods
        function obj = BaseRaw(EEG, picks)
            arguments
                EEG
                picks = (1: size(EEG.data, 1))
            end

            if nargin >= 1
                obj.data = EEG.data(picks, :);
                obj.fs = EEG.srate;

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
                    obj.chanlocs = EEG.chanlocs;
                else
                    obj.chanlocs = [];
                end
                obj.last_samp = length(obj);
                obj.time_as_index = (1 : length(obj));
            end

        end
        
        function [r, meta] = windowedPower(obj,  time_window, noverlap, ...
                lfreq, hfreq, verbose)
            % calculate powers in a windowed fashion
            arguments
                obj
                time_window (1,1) double
                noverlap (1,1) double
                lfreq (1,1) = 8
                hfreq (1,1) = 13
                verbose logical = true
            end
            
            rowsNo = fix((size(obj.data, 2) - obj.fs*time_window) / ...
                (obj.fs*time_window - obj.fs*noverlap)) + 1;
            matrix = zeros(rowsNo, size(obj.data,1));

            for j = 1:rowsNo
                begin = (time_window*obj.fs - noverlap*obj.fs) * (j-1) + 1;
                stop = time_window*obj.fs + ...
                    (time_window*obj.fs - noverlap*obj.fs) * (j-1) + 1;
                
                powers = obj.sum_power_segment((begin:stop), lfreq, hfreq);
                matrix(j, :) = powers;
            end

            r = matrix;
            meta = struct('time_window', time_window, 'noverlap', noverlap);

            if verbose
                disp(meta);
            end
        end

        function filterEEG(obj,hi,lo)
            obj.raw = obj.data;
            if ~isempty(hi)
                   obj.data = hifi(obj.data', 1e6/obj.fs, hi)';
            end
            if ~isempty(lo)
                   obj.data = lofi(obj.data', 1e6/obj.fs, lo)';
            end
        end


        
        function [r, meta] = windowedPower(obj,  time_window, noverlap, ...
                lfreq, hfreq, verbose)
            % calculate powers in a windowed fashion
            arguments
                obj
                time_window (1,1) double
                noverlap (1,1) double
                lfreq (1,1) = 8
                hfreq (1,1) = 13
                verbose logical = true
            end
            
            rowsNo = fix((size(obj.data, 2) - obj.fs*time_window) / ...
                (obj.fs*time_window - obj.fs*noverlap)) + 1;
            matrix = zeros(rowsNo, size(obj.data,1));

            for j = 1:rowsNo
                begin = (time_window*obj.fs - noverlap*obj.fs) * (j-1) + 1;
                stop = time_window*obj.fs + ...
                    (time_window*obj.fs - noverlap*obj.fs) * (j-1) + 1;
                
                powers = obj.sum_power_segment((begin:stop), lfreq, hfreq);
                matrix(j, :) = powers;
            end

            r = matrix;
            meta = struct('time_window', time_window, 'noverlap', noverlap);

            if verbose
                disp(meta);
            end
        end

        function r = sum_power(obj, l_freq, h_freq)
            r = sum_freq_band(obj.psd, obj.freq, l_freq, h_freq);
        end
        function r = sum_power_segment(obj, segment, l_freq, h_freq)
            transposed = obj.data';
            % segment in datapoints as ``double``
            art_mat = abs(transposed(segment,:)) > 100;
            art_vec = ~any(art_mat);
            psd_window = NaN(129,64);
            freq_window = NaN(129,1);
            r = NaN(1,64);
            if sum(art_vec) > 0
            [psd_window(:,art_vec), freq_window] = pwelch(transposed(segment, art_vec), ...
                [],[],256,obj.fs);
            
            r = sum_freq_band(psd_window, freq_window, l_freq, h_freq);
            end
        end
        function crop(obj, tmin, tmax)
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
            for i = 1:64
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
%             r = reshape(double( 1/obj.fs : length(obj.data) ) * 1/obj.fs, 1, []);
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