function NOFRF_plot_Un_Y_vfreq(u, Y, N, Ts, Fs, len_adj, len_adj_hlf, freq_rng_all, rmv_U_mat, harm_inpt)

w_1 = [0:Fs/len_adj:(Fs/2)];%[0:Fs/len_adj:(Fs/2)+(Fs/len_adj)];%
w = w_1(1:len_adj_hlf);
x_rng = [0,max(w(freq_rng_all))];

%% ----------------------------- U-mat -------------------------------

U_mat = zeros(len_adj_hlf,N);
for i = 1:N
    U_fft = fft(u.^i,len_adj).*Ts.*((1/sqrt(i))/((2*pi)^(i-1)));
    U_mat(:,i) = U_fft(1:len_adj_hlf);
end


N_hlf = floor(N/2);
N_rem = N - 2*N_hlf;

figure;
for i = 1:N
    subplot(N_hlf+N_rem, 2 , i);
    if harm_inpt == 0
        plot(w,abs(U_mat(:,i)));axis([x_rng -inf inf]);
    else
        stem(w(logical(rmv_U_mat(:,i))) ,  abs( U_mat( logical(rmv_U_mat(:,i)) ,i) ) );axis([x_rng -inf inf]);
    end
    title(['U',num2str(i)]);
end
suptitle('Mag \it U_n(j\omega)\rm');

figure;
for i = 1:N
    subplot(N_hlf+N_rem, 2 , i);
    if harm_inpt == 0
        plot(w,(angle(U_mat(:,i)).* rmv_U_mat(:,i)).*(180/pi) );axis([x_rng -inf inf]);
    else
        stem(w(logical(rmv_U_mat(:,i))) ,  angle(U_mat(logical(rmv_U_mat(:,i)) ,i )).*(180/pi) );axis([x_rng -inf inf]);
    end
    title(['U',num2str(i)]);
end
suptitle('Phase \it U_n(j\omega)\rm');

%% ----------------------------- Y-mat -------------------------------

Y_mat = zeros(len_adj_hlf,n_A);% matrix containing ffts of Y1,...,Yn_a

for i = 1:n_A
    Y_fft = fft(Y(:,i),len_adj).*Ts;
    Y_mat(:,i) = Y_fft(1:len_adj_hlf);
end

figure;
for i = 1:n_A
    subplot(n_A,1,i);plot(w_adj,abs(Y_mat(:,i)));title(['Amp = ',num2str(A(i))]);
end

figure;
for i = 1:n_A
    subplot(n_A,1,i);plot(tspan,Y(:,i));title(['Amp = ',num2str(A(i))]);
end

%% -------------------- NOFRF valid freq ranges -------------------

figure;
for i = 1:N
    subplot(N,1,i);area(w,rmv_U_mat(:,i));axis([x_rng -inf inf]);
end

figure;subplot(2,1,1);plot(w_adj,abs(U_mat(1:len_adj_hlf,1)));
subplot(2,1,2);plot(w_adj,abs(Y_mat(:,n_A)));



end