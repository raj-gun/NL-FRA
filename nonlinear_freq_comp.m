function [freq_rng_nonlinear_ord,freq_rng_nonlinear_ord_NOFRF] = nonlinear_freq_comp(N,pos_freq_comp,len_fft,Fs)%#codegen

%%
% N - maximum nonlinear order,
% pos_freq_comp - vector containing the positive discrete frequency components of the input
% This code calculates the discrete output frequency components of a nonlinear system with order N
% when input has a discrete frequency components

%%
if size(pos_freq_comp,1) ~= 1 pos_freq_comp=pos_freq_comp.';end %Make pos_freq_comp to a row vector

freq_comp = [-1.*fliplr(pos_freq_comp),pos_freq_comp];
len_freq_comp = length(freq_comp);

[unq_nl_comb] = nl_term_comb(N,len_freq_comp);

freq_rng_nonlinear_ord = cell(N,1);
freq_rng_nonlinear_ord_NOFRF = cell(N,1);

for n = 1:N
    
    if n == 1
        freq_rng_nonlinear_ord{n,1} = pos_freq_comp;
    else
        n_freq_comp = sum(freq_comp(unq_nl_comb{n,1}),2);
        n_freq_comp = sort(n_freq_comp(find(n_freq_comp >= 0)),'ascend');
        n_freq_comp = unique(n_freq_comp);
        freq_rng_nonlinear_ord{n,1} = n_freq_comp;
    end
    
    %calculate the corresponding fft indexe ranges
    if len_fft ~=0 && Fs ~=0  
        freq_rng_nonlinear_ord_NOFRF{n,1} = (round(freq_rng_nonlinear_ord{n,1}.* (len_fft/Fs)) ) + ones(size(freq_rng_nonlinear_ord{n,1}));  
    end

end

end