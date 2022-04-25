fnames = dir('*.mat');
% retrieve size of loading data w/o loading it to memory
data_size = whos('-file', fnames(end).name).size();
matrix = zeros(length(fnames), data_size(1), data_size(2));
%%
ntrials   = 0;     % number of trials
nobs      = size(matrix,3);   % number of observations per trial
regmode   = 'LWR';  % VAR model estimation regression mode ('OLS', 'LWR' or empty for default)
icregmode = 'LWR';  % information criteria regression mode ('OLS', 'LWR' or empty for default)
morder    = 'AIC';  % model order to use ('actual', 'AIC', 'BIC' or supplied numerical value)
momax     = 3;     % maximum model order for model order estimation
acmaxlags = [];   % maximum autocovariance lags (empty for automatic calculation)
tstat     = 'F';     % statistical test for MVGC:  'F' for Granger's F-test (default) or 'chi2' for Geweke's chi2 test
alpha     = 0.05;   % significance level for significance test
mhtc      = 'FDR';  % multiple hypothesis test correction (see routine 'significance')
fs        = 250;    % sample rate (Hz)
fres      = 250;     % frequency resolution (empty for automatic calculation)
seed      = 0;

bands = struct('delta', {1, 4}, 'alpha', {7, 13}, 'beta', {13, 30});
%%
% SET PARAMS
event_time = 120; % in seconds
patient_select = [13];
from_electrodes = [(2:4), 18];
time_window = 2; % in seconds
noverlap = 0; % in seconds
k_neighbors = 5; % k-nearest neighbors for smoothing
band = 'delta'; % choose from ['delta' | 'alpha' | 'beta']
%%
for i = patient_select
    path = fullfile(fnames(i).folder, fnames(i).name);
    matrix(i, :, :) = importdata(path);
end

rowsNo = fix((size(matrix, 3) - fs*time_window) / (fs*time_window - fs*noverlap)) + 1;
matrix = matrix(:, :, 1:size(matrix,3)-rem(size(matrix,3), time_window*fs));
matrix = reshape(matrix, length(fnames), data_size(1), [], rowsNo);
% Z-scale our data across epochs (3rd dim)
matrix = normalize(matrix, 3);
%%
% bivariate granger
event = event_time / time_window;
traces = nan(rowsNo,1);
pairs = nchoosek(from_electrodes, 2);
lfreq = getfield(bands(1), band);
hfreq = getfield(bands(2), band);
f = figure();

for k = 1:size(pairs, 1)
    % pick channel pair
    matrix_pick = matrix(patient_select, pairs(k,:), :, :);

    slices = zeros(rowsNo, size(matrix_pick, 2)^2 - size(matrix_pick, 2));
    grangers_time = zeros(2,2,rowsNo);
    Fints_time = zeros(2,2,rowsNo);
    for slice_id = 1:size(matrix_pick, 4)
        grangers = nan(size(matrix_pick, 2));
        Fints = nan(size(matrix_pick, 2));
        for patient_id = 1:size(matrix_pick, 1)
            data = squeeze(matrix_pick(patient_id, :, :, slice_id));
            [A,SIG] = tsdata_to_var(data,momax);

            f = var_to_spwcgc(A,SIG,fres);
            Fint = smvgc_to_mvgc(f, ...
                [ ...
                lfreq / (fs/2), ...
                hfreq / (fs/2) ...
                ]);

            [F, ~] = var_to_pwcgc(A,SIG,data,regmode,tstat);
            grangers = cat(3, grangers, F);
            Fints = cat(3, Fints, Fint);

        end
        grangers_mean = mean(grangers, 3, 'omitnan');
        grangers_time(:, :, slice_id) = grangers_mean;
        Fints_mean = median(Fints, 3, 'omitnan');
        Fints_time(:, :, slice_id) = Fints_mean;
    end

% PLOTS

    for i = 1:2
        for j = 1:2
            datplot = squeeze(Fints_time(j, i, :));
            % normalize
            datplot_scaled = datplot / prctile(datplot(1:event), 95);
            % apply moving average to smooth signal
            datplot_scaled_movm = movmean(datplot_scaled, k_neighbors);
            % plot log-y
            semilogy(datplot_scaled_movm, Color='black')
            hold on
            if ~all(isnan(datplot_scaled_movm))
                traces = cat(2, traces, datplot_scaled_movm);
            end
        end
    end
end

xline(event)
title(band, ('across '+string(length(patient_select))+' patient(s)'))
plot(median(traces, 2, 'omitnan'), LineWidth=3, Color=[1.0000 0.6471 0])
hold off