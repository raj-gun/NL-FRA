function [Norm_SSE_abs_LS_2,Norm_SSE_arg_LS_2,Norm_SSE_LS_2,Fe_n_Yn,Fe_p_Yn,Y_NOFRF_MLS_n,Y_NOFRF_LS_2,G_LS_2,sse,Y_model,rmv_U_mat] = ...
    general_nofrf_comp(Fs,Ts,tspan,fftn,pos_freq_comp,u_nofrf,A,Amp,zero_cross,nl_ord_set,gc,harm_inpt,lw,norm,displ,model)
%------------
len = length(u_nofrf);
n_A = length(A);
Y = zeros(len,n_A);

for i = 1:n_A
    Y(:,i) = model(i).';
end

%-------------------------------

if zero_cross == 1
    [~,zc_tspan] = min(abs(tspan));
    %tspan = tspan(zc_tspan:end);
    u_nofrf = u_nofrf(zc_tspan:end);
    Y = Y(zc_tspan:end,:);
end

% len = length(u_nofrf);
%len_adj = len;
len_adj = fftn;
% len_adj_hlf = floor(len_adj/2)+1;
%w = [0:Fs/len_adj:(Fs/2)];%[0:Fs/len_adj:(Fs/2)+(Fs/len_adj)];%

% disp(['Inpt Amps = ',num2str(A')]);
% disp('Simulations complete');
%%
N = max(nl_ord_set);
% displ = 1;
[G_LS_2,~,~,~,rmv_U_mat,~,len_adj_hlf,sse] = NOFRF_LS2comp(u_nofrf,Y,A,N,Fs,Ts,pos_freq_comp,nl_ord_set,fftn,displ(1));

%% ------------- NOFRF valid freq ranges -------------------
%[freq_rng_nonlinear_ord,freq_rng_nonlinear_ord_NOFRF] = nonlinear_freq_range(N,f2,f1,len_adj,Fs);
% rmv_u_mat_full_len = freq_rng_nonlinear_ord_NOFRF{1,end}(end);
% rmv_U_mat_full = zeros(rmv_u_mat_full_len,N);
% for i = 1:N
%     rng_rmv_U_mat = freq_rng_nonlinear_ord_NOFRF{i,1};
%     rmv_U_mat_full(rng_rmv_U_mat,i) = 1; % nth order NOFRF valid freq ranges
% end

% if len_adj_hlf >= rmv_u_mat_full_len
%     len_adj_hlf = rmv_u_mat_full_len;
% end
% rmv_U_mat = zeros(len_adj_hlf,N);
% rmv_U_mat = rmv_U_mat_full(1:len_adj_hlf,:);

freq_rng_all = sum(rmv_U_mat,2) ~= 0; %all NOFRFs valid ranges

w = 0:Fs/len_adj:(Fs/2);
w_adj = w(1:len_adj_hlf);
%% ----------------------------------------------------------
G_LS_2 = G_LS_2.';

% Fe_n_nofrf = sum(abs(G_LS_2),1)./sum(sum(abs(G_LS_2),1));
% Fe_p_nofrf = sum(abs(G_LS_2),1);
%% NOFRFs testing

y_test = model(Amp).';

% disp('Simulations complete');

if zero_cross == 1
    [~,zc_tspan] = min(abs(tspan));
    tspan = tspan(zc_tspan:end);
    u_nofrf = u_nofrf(zc_tspan:end);
    y_test = y_test(zc_tspan:end);
end
% %--------------------------

%% Plots
% -------------------- Y_pred -----------------------------
% figure;plot(tspan_org,y_test,'b');%axis([0 inf -inf inf]);

[Norm_SSE_abs_LS_2,Norm_SSE_arg_LS_2,Norm_SSE_LS_2,Fe_n_Yn,Fe_p_Yn,Y_NOFRF_MLS_n,Y_NOFRF_LS_2,Y_model] = ...
    NOFRF_Y_pred_comp(pos_freq_comp, u_nofrf, y_test, Ts, Fs, N, nl_ord_set, Amp, len_adj, len_adj_hlf, G_LS_2, rmv_U_mat, freq_rng_all, gc, harm_inpt, lw, norm, displ(2));

% -------------------- Plot NOFRFs -----------------------------
NOFRF_plots(G_LS_2, N, nl_ord_set, w_adj, freq_rng_all, rmv_U_mat, gc, harm_inpt, lw, norm);
% NOFRF_plplots(G_LS_2, N, rmv_U_mat, gc, lw, norm);

% ----------------------------- U-mat, Y-mat & Valid freq -------------------------------
if displ(3) == 1
    NOFRF_plot_Un_Y_vfreq(tspan, u_nofrf, Y, y_test, Y_NOFRF_LS_2, A, n_A, N, nl_ord_set, Ts, Fs, len_adj, len_adj_hlf, freq_rng_all, rmv_U_mat, harm_inpt, lw);
end
%%
end
