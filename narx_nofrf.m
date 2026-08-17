function [Norm_SSE_abs_LS_2,Norm_SSE_arg_LS_2,Norm_SSE_LS_2,Fe_n_Yn,Fe_p_Yn,Y_NOFRF_MLS_n,Y_NOFRF_LS_2,G_LS_2,sse,Y_model] = ...
    narx_nofrf(Fs,Ts,tspan,fftn,f1,f2,u_nofrf,A,Amp,zero_cross,nl_ord_set,gc,harm_inpt,lw,norm,displ,model)
%------------
theta_procss=model{1};
trm_chsn_lin=model{2};
trm_chsn_lin_org=model{3};
n_lin_trms_org=model{4};
nl_ord_max=model{5};
trm_chsn_nl=model{6};
bias=model{7};
inpt0=model{8};
min_dyn_ord_u=model{9};
max_dyn_ord_u=model{10};
min_dyn_ord_y=model{11};
max_dyn_ord_y=model{12};
%------------
u_org=u_nofrf;
len = length(u_nofrf);
n_A = length(A);
Y = zeros(len,n_A);

for i = 1:n_A
    %%Information matrix formulation
    %---------------------------------------
    n_terms_y = sum(max_dyn_ord_y);
    %---------------------------------------
    [~,~,X,~,~] = info_mat_sysID(inpt0,max_dyn_ord_u,max_dyn_ord_y, (A(i).*u_nofrf) ,u_nofrf.*0);
    X = [ X(:,min_dyn_ord_y:max_dyn_ord_y), X(:,n_terms_y+min_dyn_ord_u:end) ];
    n_terms_y = (max_dyn_ord_y - min_dyn_ord_y)+1;
    
    %%Simulation - OSP and MPO
    %---------------------------
    Xsim_ID = X;%X_ID;%X_trim_ID2;
    Sim_ini_val = Xsim_ID(1,1:n_terms_y);
    U_delay_mat_sim = Xsim_ID(:,n_terms_y+1:end);
    %---------------------------
    if length(Sim_ini_val) == 1
        if Sim_ini_val == 1 || Sim_ini_val == 0
            y_lag_lin_srt = Xsim_ID(1,1:n_terms_y).*Sim_ini_val;
        else
            y_lag_lin_srt = Sim_ini_val;
        end
    else
        y_lag_lin_srt = Sim_ini_val;
    end
    Y(:,i) = sim_model_reg_2(theta_procss,trm_chsn_lin,trm_chsn_lin_org,U_delay_mat_sim,n_lin_trms_org,y_lag_lin_srt,nl_ord_max,trm_chsn_nl,bias);
    
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
[G_LS_2,~,~,~,rmv_U_mat,~,len_adj_hlf,sse] = NOFRF_LS2(u_nofrf,Y,A,N,Fs,Ts,f1,f2,nl_ord_set,fftn,displ(1));

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

%w = 0:Fs/len_adj:(Fs/2);
% w_adj = w(1:len_adj_hlf);
%% ----------------------------------------------------------
G_LS_2 = G_LS_2.';

% Fe_n_nofrf = sum(abs(G_LS_2),1)./sum(sum(abs(G_LS_2),1));
% Fe_p_nofrf = sum(abs(G_LS_2),1);
%% NOFRFs testing

%--------------------------
%tspan_org = tspan;
u_nofrf = u_org;

%%Information matrix formulation
%---------------------------------------
n_terms_y = sum(max_dyn_ord_y);
%---------------------------------------
[~,~,~,X,~] = info_mat_sysID(inpt0,max_dyn_ord_u,max_dyn_ord_y, (Amp.*u_nofrf) ,u_nofrf.*0);
X = [ X(:,min_dyn_ord_y:max_dyn_ord_y), X(:,n_terms_y+min_dyn_ord_u:end) ];
n_terms_y = (max_dyn_ord_y - min_dyn_ord_y)+1;

%%Simulation - OSP and MPO
%---------------------------
Xsim_ID = X;%X_ID;%X_trim_ID2;
Sim_ini_val = Xsim_ID(1,1:n_terms_y);
U_delay_mat_sim = Xsim_ID(:,n_terms_y+1:end);
%---------------------------
if length(Sim_ini_val) == 1
    if Sim_ini_val == 1 || Sim_ini_val == 0
        y_lag_lin_srt = Xsim_ID(1,1:n_terms_y).*Sim_ini_val;
    else
        y_lag_lin_srt = Sim_ini_val;
    end
else
    y_lag_lin_srt = Sim_ini_val;
end
y_test = sim_model_reg_2(theta_procss,trm_chsn_lin,trm_chsn_lin_org,U_delay_mat_sim,n_lin_trms_org,y_lag_lin_srt,nl_ord_max,trm_chsn_nl,bias);

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
    NOFRF_Y_pred(f1, f2, u_nofrf, y_test, Ts, Fs, N, Amp, len_adj, len_adj_hlf, G_LS_2, rmv_U_mat, freq_rng_all, gc, harm_inpt, lw, norm, displ(2));

% -------------------- Plot NOFRFs -----------------------------
% NOFRF_plots(G_LS_2, N, w_adj, freq_rng_all, rmv_U_mat, gc, harm_inpt, lw, norm);
% NOFRF_plplots(G_LS_2, N, rmv_U_mat, c, lw);

% ----------------------------- U-mat, Y-mat & Valid freq -------------------------------
if displ(3) == 1
    NOFRF_plot_Un_Y_vfreq(tspan, u_nofrf, Y, y_test, Y_NOFRF_LS_2, A, n_A, N, Ts, Fs, len_adj, len_adj_hlf, freq_rng_all, rmv_U_mat, harm_inpt, lw);
end
%%
end
