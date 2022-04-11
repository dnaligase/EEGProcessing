time_window = 5;
noverlap = 3;

matrix = zeros(rowsNo, size(patient.data, 1), length(ALLEEG));



for i = 1:length(ALLEEG)
    EEG = ALLEEG(i);
    patient = BaseRaw(EEG);
    patient.crop(0, 60);
    rowsNo = fix((size(patient.data,2) - fs*time_window) / (fs*time_window - fs*noverlap)) + 1;

    for j = 1:rowsNo

        begin = (time_window*fs - noverlap*fs) * (j-1) + 1;
        stop = time_window*fs + (time_window*fs - noverlap*fs) * (j-1) + 1;
        powers = patient.sum_power_segment((begin:stop), 8, 13);
        matrix(j, :, i) = powers;
    end
end

matrix_mean = mean(matrix, 3);