function [Yn_NOFRF,Y_NOFRF,len_adj_hlf] = nofrf_test(Amp,u,G_LS_2,len_adj,len_adj_hlf,N,Fs,Ts,f1,f2)

%-----------------------------------------------
[freq_rng_nonlinear_ord,freq_rng_nonlinear_ord_NOFRF] = nonlinear_freq_range(N,f2,f1,len_adj,Fs);
rmv_u_mat_full_len = freq_rng_nonlinear_ord_NOFRF{N,1}(end);
rmv_U_mat_full = zeros(rmv_u_mat_full_len,N);
for i = 1:N
    rng_rmv_U_mat = freq_rng_nonlinear_ord_NOFRF{i,1};
    rmv_U_mat_full(rng_rmv_U_mat,i) = 1; % nth order NOFRF valid freq ranges
end

if len_adj_hlf >= rmv_u_mat_full_len 
    len_adj_hlf = rmv_u_mat_full_len;
end
rmv_U_mat = zeros(len_adj_hlf,N);
rmv_U_mat = rmv_U_mat_full(1:len_adj_hlf,:);

freq_rng_all = sum(rmv_U_mat,2) ~= 0; %all NOFRFs valid ranges
%-----------------------------------------------

U_mat = zeros(len_adj_hlf,N);
A_mat = zeros(1,N);
for i = 1:N
    U_fft = fft(u.^i,len_adj).*Ts.*((1/sqrt(i))/((2*pi)^(i-1)));
    U_mat(:,i) = U_fft(1:len_adj_hlf);
    A_mat(i) = Amp^i;
end

Yn_NOFRF = G_LS_2 .* (U_mat .* repmat(A_mat,len_adj_hlf,1) .* rmv_U_mat);
Y_NOFRF = sum(Yn_NOFRF ,2);

end