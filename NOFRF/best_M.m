function [M]=best_M(Mi,Mstp,Me,u_nofrf,model)
u_org=u_nofrf;
Mn=zeros(length(Mi:Mstp:Me),2);

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
cnt=1;
for M = Mi:Mstp:Me
    u_nofrf = u_org.*M;
    %%Information matrix formulation
    %---------------------------------------
    n_terms_y = sum(max_dyn_ord_y);
    %---------------------------------------
    [~,~,X,~,~] = info_mat_sysID(inpt0,max_dyn_ord_u,max_dyn_ord_y, (u_nofrf)' ,u_nofrf'.*0);
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
    
    Mn(cnt,:) = [isnan(sum(y_test)),M];
    cnt=cnt+1;
end

Mn=Mn(Mn(:,1)==0,:);
M=Mn(end,2);


end