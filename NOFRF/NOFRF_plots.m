function NOFRF_plots(G_LS_2, N, nl_ord_set, w_adj, freq_rng_all, rmv_U_mat, c, harm_inpt, lw, norm)

N_hlf = floor(N/2);
N_rem = N - 2*N_hlf;
w_G = {w_adj(2:end),abs(G_LS_2(2:end,:))};
G_LS_2_log_mag = log(abs(G_LS_2(2:end,:))); G_LS_2_log_mag(isinf(G_LS_2_log_mag)) = 0;
G_LS_2_log_mag_freq = G_LS_2_log_mag(freq_rng_all(2:end),:);
% G_LS_2_mag_freq = abs(G_LS_2);
% G_LS_2_mag_freq = G_LS_2_mag_freq(freq_rng_all,:);
% G_LS_2_mag_freq_norm = G_LS_2_mag_freq./max(G_LS_2_mag_freq')';
G_LS_2_freq = G_LS_2(freq_rng_all(2:end),:);
w_freq = w_G{1}(freq_rng_all(2:end));
x_rng = [0,max(w_freq)];

len_freq = length(w_freq); len_freq_hlf = floor(len_freq/2);
len_freq_rem = len_freq - 2*len_freq_hlf;

% G_LS_2_mag = abs(G_LS_2);
% G_LS_2_phase = angle(G_LS_2);
% G_LS_2_phase_deg = G_LS_2_phase.*(180/pi);
% G_LS_2_phase_deg_sum = sum(G_LS_2_phase_deg,1);

% c= 'r';%[0.4660 0.6740 0.1880];%[0.4940 0.1840 0.5560];%[0.9290 0.6940 0.1250];

% figure;
% for i = 1:N
%     subplot(N,1,i);plot(w_adj,abs(G_LS_2(:,i)));axis([x_rng -inf inf]);
%     if i == 1
%         title('Mod-LS');
%     end
% end

% if harm_inpt == 1
%     
%     figure;
%     for i = 1:len_freq
%         subplot(len_freq_hlf+len_freq_rem, 2 , i);plot(G_LS_2_log_mag_freq(i,:),'Color',c,'Marker','o');axis([0 N -inf inf]);hold on;
%         title(['\omega',num2str(i-1),' = ',num2str(w_freq(i))]);
%     end
%     sgtitle('log( |\it G_n( j n\omega)\rm| ) vs \it n\omega \rm');
%     
%     figure;
%     for i = 1:len_freq
%         subplot(len_freq_hlf+len_freq_rem, 2 , i);stem(abs(G_LS_2_freq(i,:)),'Color',c,'LineWidth',lw);axis([0 N -inf inf]);hold on;
%         title(['\omega',num2str(i-1),' = ',num2str(w_freq(i))]);
%     end
%     sgtitle('|\it G_n( j n\omega)\rm| vs \it n\omega \rm');
%     
%     figure;
%     for i = 1:len_freq
%         subplot(len_freq_hlf+len_freq_rem, 2 , i);stem( angle(G_LS_2_freq(i,:)).*(180/pi) ,'Color',c,'LineWidth',lw);axis([0 N -inf inf]);hold on;
%         title(['\omega',num2str(i-1),' = ',num2str(w_freq(i))]);
%     end
%     sgtitle('\angle( \it G_n( j n\omega)\rm ) vs \it n\omega \rm');
%     
% end



figure;
for i = nl_ord_set%1:N
    subplot(N_hlf+N_rem, 2 , i);
    if harm_inpt == 0
        plot(w_G{1},w_G{2}(:,i),'Color',c);axis([x_rng -inf inf]);hold on;
    else
        if norm == 1
            y_plt_dat = abs(w_G{2}(logical(rmv_U_mat(2:end,i)) ,i )) ./ sum(abs(w_G{2}(logical(rmv_U_mat(2:end,i)) ,i )));
        else
            y_plt_dat = abs(w_G{2}(logical(rmv_U_mat(2:end,i)) ,i ));
        end
        stem(w_G{1}(logical(rmv_U_mat(2:end,i))) ,  y_plt_dat ,'Color' , c,'LineWidth',lw);
%         plot(w_G{1}(logical(rmv_U_mat(2:end,i))) ,  y_plt_dat ,'Color' , c,'LineWidth',lw,'Marker','o','LineStyle',':');
        axis([x_rng -inf inf]);hold on;
        
    end
    title(['G',num2str(i)]);
    if mod(i,2) == 1; ylabel('Magnitude'); end %Every odd nonlinearity
    if i == nl_ord_set(find(mod(nl_ord_set, 2) == 0, 1, 'last')); xlabel('Freqeuncy (Hz)'); end %Last even nonlinerity
    if i == nl_ord_set(find(mod(nl_ord_set, 2) ~= 0, 1, 'last')); xlabel('Freqeuncy (Hz)'); end %Last odd nonlinerity
end
if norm == 1
    sgtitle('Normalised Mag-NOFRF \it G_n(j\omega)\rm');
else
    sgtitle('Mag-NOFRF \it G_n(j\omega)\rm');
end

% figure;
% for i = 1:N
%     subplot(N_hlf+N_rem, 2 , i);
%     if harm_inpt == 0
%         plot(w_G{1},G_LS_2_log_mag(:,i),'Color',c);axis([x_rng -inf inf]);hold on;
%     else
%         if norm == 1
%             y_plt_dat =  abs( G_LS_2_log_mag(logical(rmv_U_mat(2:end,i)) ,i ) ./ max(G_LS_2_log_mag(logical(rmv_U_mat(2:end,i)) ,i )) );%sign(G_LS_2_log_mag(logical(rmv_U_mat(:,i)) ,i )) .*
%         else
%             y_plt_dat = G_LS_2_log_mag(logical(rmv_U_mat(2:end,i)) ,i );
%         end
%         stem(w_G{1}(logical(rmv_U_mat(2:end,i))) ,  y_plt_dat , 'Color' , c,'LineWidth',lw);
% %         plot(w_G{1}(logical(rmv_U_mat(2:end,i))) ,  y_plt_dat , 'Color' , c,'LineWidth',lw,'Marker','o','LineStyle',':');
%         axis([x_rng -inf inf]);hold on;
%     end
%     title(['G',num2str(i)]);
% end
% if norm == 1
%     sgtitle('Normalised Log-Mag-NOFRF \it G_n(j\omega)\rm');
% else
%     sgtitle('Log-Mag-NOFRF \it G_n(j\omega)\rm');
% end

x_rng = [0,max(w_adj(freq_rng_all))];
figure;
for i = nl_ord_set%1:N
    subplot(N_hlf+N_rem, 2 , i);
    if harm_inpt == 0
        plot(w_adj,angle(G_LS_2(:,i)),'Color',c);axis([x_rng -inf inf]);hold on;
    else
        stem(w_adj(logical(rmv_U_mat(:,i))) ,  angle(G_LS_2(logical(rmv_U_mat(:,i)) ,i )).*(180/pi) ,'Color' , c,'LineWidth',lw);axis([x_rng -inf inf]);hold on;
    end
    title(['G',num2str(i)]);
    if mod(i,2) == 1; ylabel('Phase'); end %Every odd nonlinearity
    if i == nl_ord_set(find(mod(nl_ord_set, 2) == 0, 1, 'last')); xlabel('Freqeuncy (Hz)'); end %Last even nonlinerity
    if i == nl_ord_set(find(mod(nl_ord_set, 2) ~= 0, 1, 'last')); xlabel('Freqeuncy (Hz)'); end %Last odd nonlinerity
end
sgtitle('Phase-NOFRF \it G_n(j\omega)\rm');

end