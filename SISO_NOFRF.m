function [Norm_SSE_abs_LS_2,Norm_SSE_arg_LS_2,Norm_SSE_LS_2,Fe_n_Yn,Fe_p_Yn,Y_NOFRF_MLS_n,Y_NOFRF_LS_2,G_LS_2,sse,Y_model] = ...
    SISO_NOFRF(Fs,Ts,tspan,fftn,f1,f2,u_nofrf,A,Amp,nl_ord_set,gc,harm_inpt,lw,norm,displ,Y, y_test)
%% Evaluate NOFRFs

N = max(nl_ord_set);
[G_LS_2,~,~,~,rmv_U_mat,~,len_adj_hlf,sse] = NOFRF_LS2(u_nofrf,Y,A,N,Fs,Ts,f1,f2,nl_ord_set,fftn,displ(1));
freq_rng_all = sum(rmv_U_mat,2) ~= 0; %all NOFRFs valid ranges
G_LS_2 = G_LS_2.';

%% NOFRFs testing and validation

len_adj = fftn;

[Norm_SSE_abs_LS_2,Norm_SSE_arg_LS_2,Norm_SSE_LS_2,Fe_n_Yn,Fe_p_Yn,Y_NOFRF_MLS_n,Y_NOFRF_LS_2,Y_model] = ...
    NOFRF_Y_pred(f1, f2, u_nofrf, y_test, Ts, Fs, N, Amp, len_adj, len_adj_hlf, G_LS_2, rmv_U_mat, freq_rng_all, gc, harm_inpt, lw, norm, displ(2));

%% NOFRF-plotting

w = [0:Fs/len_adj:(Fs/2)];
w_adj = w(1:len_adj_hlf);

% -------------------- Plot NOFRFs -----------------------------
% NOFRF_plots(G_LS_2, N, w_adj, freq_rng_all, rmv_U_mat, gc, harm_inpt, lw, norm);
NOFRF_plots(G_LS_2, N, nl_ord_set, w_adj, freq_rng_all, rmv_U_mat, gc, harm_inpt, lw, norm);
NOFRF_plplots(G_LS_2, N, rmv_U_mat, gc, lw, norm);

% ----------------------------- U-mat, Y-mat & Valid freq -------------------------------
n_A = length(Amp);
if displ(3) == 1
    % NOFRF_plot_Un_Y_vfreq(tspan, u_nofrf, Y, y_test, Y_NOFRF_LS_2, A, n_A, N, Ts, Fs, len_adj, len_adj_hlf, freq_rng_all, rmv_U_mat, harm_inpt, lw);
    NOFRF_plot_Un_Y_vfreq(tspan, u_nofrf, Y, y_test, Y_NOFRF_LS_2, A, n_A, N, nl_ord_set, Ts, Fs, len_adj, len_adj_hlf, freq_rng_all, rmv_U_mat, harm_inpt, lw)
end
%%
end
