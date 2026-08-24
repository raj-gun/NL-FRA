function NOFRF_plplots(G_LS_2, N, rmv_U_mat, c, lw, norm)

N_hlf = floor(N/2);
N_rem = N - 2*N_hlf;


figure;
nl_ord_set = 1:N;
for i = 1:N
    subplot(2, N_hlf+N_rem, i);
    theta = angle(G_LS_2(logical(rmv_U_mat(:,i)),i));
    
    if norm == 1
        rho = abs(G_LS_2(logical(rmv_U_mat(:,i)),i)) ./ sum( abs(G_LS_2(logical(rmv_U_mat(:,i)),i)) );
    else
        rho = abs(G_LS_2(logical(rmv_U_mat(:,i)),i));
    end
    
    polarplot( theta , rho  ,'Color',c,'LineWidth',lw);
    % polarscatter( theta , rho  ,'MarkerFaceColor',c);
    
    title(['G',num2str(i)]);
end
sgtitle('Polar-NOFRF \it G_n(j\omega)\rm');


end