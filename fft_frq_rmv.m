function data = fft_frq_rmv(data_raw,freq_rmv_ind)

data_fft = fft(data_raw);
data_fft_rmv = data_fft.*freq_rmv_ind';
data = ifft(data_fft_rmv,'symmetric'); 

end