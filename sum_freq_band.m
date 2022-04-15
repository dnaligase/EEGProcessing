function datavec = sum_freq_band(psd, freq, l_freq, h_freq)
assert(size(psd, 1) == size(freq, 1), "psd and freqs dims don't match!")

idxs = freq > l_freq & freq < h_freq;
if size(psd,3) == 1
datavec = sum(psd(idxs, :), 1);
elseif size(psd,3) > 1                  % Works on 3D matrix
  datavec = sum(psd(idxs, :,:), 1); 
  datavec = squeeze(datavec)';
end