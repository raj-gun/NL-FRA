function [Norm_SSE_abs_LS_2,Norm_SSE_arg_LS_2,Norm_SSE_LS_2,Fe_n_Yn,Fe_p_Yn,Y_NOFRF_MLS_n,Y_NOFRF_LS_2,Y_model] = NOFRF_Y_pred(f1, f2, u, y, Ts, Fs, N, Amp, len_adj, len_adj_hlf, G_LS_2, rmv_U_mat, freq_rng_all, c, harm_inpt, lw, norm, draw_plt)

Y_fft = fft(y,len_adj).*Ts;
Y_vec = Y_fft(1:len_adj_hlf);

[Y_NOFRF_MLS_n,Y_NOFRF_LS_2,len_NOFRF_LS2] = ...
    nofrf_test(Amp,u,G_LS_2,len_adj,len_adj_hlf,N,Fs,Ts,f1,f2);

Y_NOFRF_LS_2 = [Y_NOFRF_LS_2;zeros(len_adj_hlf-len_NOFRF_LS2,1)];
Y_NOFRF_MLS_n = [Y_NOFRF_MLS_n;zeros(len_adj_hlf-len_NOFRF_LS2,N)]; Y_NOFRF_MLS_n=Y_NOFRF_MLS_n(2:end,:);
Error_LS_2 = (Y_vec(freq_rng_all)) - (Y_NOFRF_LS_2(freq_rng_all));
Error_abs_LS_2 = abs(Y_vec) - abs(Y_NOFRF_LS_2);
Error_arg_LS_2 = angle(Y_vec(freq_rng_all)) - angle(Y_NOFRF_LS_2(freq_rng_all));
Norm_SSE_LS_2 = sum(abs(Error_LS_2).^2,1)/var(abs(Y_vec(freq_rng_all)));
Norm_SSE_abs_LS_2 = sum(abs(Error_abs_LS_2).^2,1)/var(abs(Y_vec));
Norm_SSE_arg_LS_2 = sum(abs(Error_arg_LS_2).^2,1)/var(angle(Y_vec));
Y_model = Y_vec(freq_rng_all);

Fe_n_Yn = sum(abs(Y_NOFRF_MLS_n),1)./sum(sum(abs(Y_NOFRF_MLS_n),1));
Fe_p_Yn = sum(abs(Y_NOFRF_MLS_n),1);

%
% disp('--------------------------------------------');
% disp(['Norm SSE abs LS_2  = ',num2str(Norm_SSE_abs_LS_2)]);
% disp(['Norm SSE arg_LS_2  = ',num2str(Norm_SSE_arg_LS_2)]);
% disp(['Norm SSE LS_2      = ',num2str(Norm_SSE_LS_2)]);
% disp('--------------------------------------------');


if draw_plt == 1
    w = [0:Fs/len_adj:(Fs/2)];
    w = w(1:len_adj_hlf);
    w2=w(2:end);
    x_rng = [0,max(w(freq_rng_all))];
    figure;subplot(2,1,1);
    if harm_inpt == 0
        plot(w,abs(Y_vec));hold on;plot(w,abs(Y_NOFRF_LS_2),'r--');%set(gca, 'YScale', 'log');
    else
        stem(w(logical(freq_rng_all)) ,  abs( Y_vec(logical(freq_rng_all)) ) ,'LineWidth',lw);axis([x_rng -inf inf]);hold on;
        stem(w(logical(freq_rng_all)) ,  abs( Y_NOFRF_LS_2(logical(freq_rng_all)) ), 'Color' , 'r','LineWidth',lw);
    end
    title('Mod-LS');
    subplot(2,1,2);
    if harm_inpt == 0
        plot(w,angle(Y_vec).*(180/pi));hold on;plot(w,angle(Y_NOFRF_LS_2).*(180/pi),'r--');
    else
        stem(w(logical(freq_rng_all)) ,  angle( Y_vec(logical(freq_rng_all)) ).*(180/pi) ,'LineWidth',lw);axis([x_rng -inf inf]);hold on;
        stem(w(logical(freq_rng_all)) ,  angle( Y_NOFRF_LS_2(logical(freq_rng_all)) ).*(180/pi), 'Color' , 'r','LineWidth',lw);
    end
    
    Y_NOFRF_MLS_n_t = zeros(size(Y_NOFRF_MLS_n));
    %x_rng = [0,inf];
    
    N_hlf = floor(N/2);
    N_rem = N - 2*N_hlf;
    
    Y_NOFRF_MLS_n_log_mag = log(abs(Y_NOFRF_MLS_n)); Y_NOFRF_MLS_n_log_mag(isinf(Y_NOFRF_MLS_n_log_mag)) = 0;
    %Y_NOFRF_MLS_n_mag_freq = Y_NOFRF_MLS_n_log_mag(freq_rng_all,:);
    
    figure;
    for i = 1:N
        subplot(N_hlf+N_rem, 2 , i);
        if norm == 0
            if harm_inpt == 0
                plot(w2,abs(Y_NOFRF_MLS_n(:,i)), 'Color' , c);axis([x_rng -inf inf]);
            else
                stem(w2(logical(rmv_U_mat(2:end,i))) ,  abs( Y_NOFRF_MLS_n( logical(rmv_U_mat(2:end,i)) ,i) ) , 'Color' , c,'LineWidth',lw);axis([x_rng -inf inf]);hold on;
            end
            title(['Y',num2str(i)]);
            Y_NOFRF_MLS_n_t(:,i) = ifft(Y_NOFRF_MLS_n(:,i),'symmetric');
        else
            Y_NOFRF_MLS_n_norm = abs(Y_NOFRF_MLS_n) ./ max(max(abs(Y_NOFRF_MLS_n)));
            if harm_inpt == 0
                plot(w2,Y_NOFRF_MLS_n_norm(:,i), 'Color' , c);axis([x_rng -inf inf]);
            else
                stem(w2(logical(rmv_U_mat(2:end,i))) , Y_NOFRF_MLS_n_norm( logical(rmv_U_mat(2:end,i)) ,i), 'Color' , c,'LineWidth',lw);axis([x_rng -inf inf]);hold on;
            end
            title(['Y-norm',num2str(i)]);
            Y_NOFRF_MLS_n_t(:,i) = ifft(Y_NOFRF_MLS_n(:,i),'symmetric');
        end
    end
    sgtitle('Mag-OFRF \it Y_n(j\omega)\rm');
    
    figure;
    for i = 1:N
        subplot(N_hlf+N_rem, 2 , i);
        if harm_inpt == 0
            plot(w2,abs(Y_NOFRF_MLS_n(:,i)), 'Color' , c);axis([x_rng -inf inf]);set(gca, 'YScale', 'log');
            %plot(w(logical(rmv_U_mat(:,i))),Y_NOFRF_MLS_n_log_mag(logical(rmv_U_mat(:,i)),i), 'Color' , c);axis([x_rng -inf inf]);hold on;
        else
            stem(w2(logical(rmv_U_mat(2:end,i))) ,  Y_NOFRF_MLS_n_log_mag(logical(rmv_U_mat(2:end,i)) ,i ), 'Color' , c,'LineWidth',lw);axis([x_rng -inf inf]);hold on;
        end
        title(['Y',num2str(i)]);
    end
    sgtitle('Log-Mag-OFRF \it Y_n(j\omega)\rm');
    
    figure;
    for i = 1:N
        subplot(N_hlf+N_rem, 2 , i);
        if harm_inpt == 0
            plot(w2,angle(Y_NOFRF_MLS_n(:,i)), 'Color' , c);axis([x_rng -inf inf]);
        else
            stem(w2(logical(rmv_U_mat(2:end,i))) ,  angle(Y_NOFRF_MLS_n(logical(rmv_U_mat(2:end,i)) ,i )).*(180/pi) , 'Color' , c,'LineWidth',lw);axis([x_rng -inf inf]);hold on;
        end
        title(['Y',num2str(i)]);
    end
    sgtitle('Phase-OFRF \it Y_n(j\omega)\rm');
    
    %figure;plot(tspan,Y);hold on;plot([0:Ts*2:20],sum(Y_NOFRF_RLS_n_t,2),'r--');
    
end

end