function datavec = sum_freq_band(psd, freq, l_freq, h_freq)
assert(size(psd, 1) == size(freq, 1), "psd and freqs dims don't match!")

idxs = freq > l_freq & freq < h_freq;
datavec = sum(psd(idxs, :), 1);
