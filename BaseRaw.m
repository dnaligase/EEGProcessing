classdef BaseRaw < handle
    properties
        data
        subject
        times
        psd
        freq
    end
    properties (SetAccess = private)
        fs
        chanlocs
        first_samp = 1
        last_samp
    end
    properties (Hidden = true)
        time_as_index
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
                obj.subject = EEG.subject;

                if isempty(EEG.times)
                    obj.times = obj.create_times();
                else
                    % convert EEGLAB times to ms
                    obj.times = EEG.times / 1000;
                end

                [obj.psd, obj.freq] = pwelch(obj.data',[],[],256,obj.fs);
                if ~isempty(EEG.chanlocs)
                    obj.chanlocs = EEG.chanlocs;
                else
                    obj.chanlocs = [];
                end
                obj.last_samp = length(obj);
                obj.time_as_index = (1 : length(obj));
            end

        end

        function r = sum_power(obj, l_freq, h_freq)
            r = sum_freq_band(obj.psd, obj.freq, l_freq, h_freq);
        end
        function r = sum_power_segment(obj, segment, l_freq, h_freq)
            transposed = obj.data';
            % segment in datapoints as ``double``
            [psd_window, freq_window] = pwelch(transposed(segment, :), ...
                [],[],256,obj.fs);
            
            r = sum_freq_band(psd_window, freq_window, l_freq, h_freq);
        end
        function crop(obj, tmin, tmax)
            arguments
                obj
                tmin {mustBeGreaterThanOrEqual(tmin, 0)}
                tmax {mustBeReal}
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
            figure;
            plot3([obj.chanlocs.X],[obj.chanlocs.Y],[obj.chanlocs.Z], 'o');
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
            r = reshape(double( 0 : length(obj.data) ) * 1/obj.fs, 1, []);
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