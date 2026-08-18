function [G_LS_2,Y_LS_2,U_mat,Y_mat,rmv_U_mat,freq_rng_nonlinear_ord,len_adj_hlf,sse] = NOFRF_LS2(u,Y,A,N,Fs,Ts,f1,f2,nl_ord_set,fftn,displ)
% All inputs are tall matrices or vectors except for u
warning('off','all');
if f2 > f1
    disp('f1 should be greater than f2');
    G_LS_2=[];Y_LS_2=[];U_mat=[];Y_mat=[];rmv_U_mat=[];
    return;
end

% dat_len = length(u);
len_adj = fftn;%2^nextpow2(len);
len_adj_hlf = floor(len_adj/2)+1;
n_A = length(A);

%------------- NOFRF valid freq ranges -------------------
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
%----------------------------------------------------------
%-------------- Modify rmv_U_mat --------------------------
rmv_U_mat_and = zeros(size(rmv_U_mat));
rmv_U_mat_and(:,nl_ord_set) = 1;
rmv_U_mat_mod = rmv_U_mat & rmv_U_mat_and;
rmv_U_mat = rmv_U_mat_mod;
%----------------------------------------------------------

U_mat = zeros(len_adj_hlf,N);% matrix containing ffts of u,...,u^N
U_mat_real = zeros(len_adj_hlf,N);% matrix containing real part of the ffts of u,...,u^N
U_mat_imag = zeros(len_adj_hlf,N);% matrix containing imaginary part of the ffts of u,...,u^N

Y_mat = zeros(len_adj_hlf,n_A);% matrix containing ffts of Y1,...,Yn_a
Y_mat_real = zeros(len_adj_hlf,n_A);% matrix containing real parts of the ffts of Y1,...,Yn_a
Y_mat_imag = zeros(len_adj_hlf,n_A);% matrix containing imaginay parts of the ffts of Y1,...,Yn_a

A_mat = zeros(n_A,N); % matrix containing A,...,A^N values

G_mat_real = zeros(N,len_adj_hlf);
G_mat_imag = zeros(N,len_adj_hlf);

for i = 1:N
    U_fft = fft(u.^i,len_adj).*Ts.*((1/sqrt(i))/((2*pi)^(i-1)));
    U_mat(:,i) = U_fft(1:len_adj_hlf);
    U_mat_real(:,i) = real(U_mat(:,i));
    U_mat_imag(:,i) = imag(U_mat(:,i));
    A_mat(:,i) = A.^i;
end

for i = 1:n_A
    Y_fft = fft(Y(:,i),len_adj).*Ts;
    Y_mat(:,i) = Y_fft(1:len_adj_hlf);
    Y_mat_real(:,i) = real(Y_mat(:,i));
    Y_mat_imag(:,i) = imag(Y_mat(:,i));
end

Y_LS_2 = zeros(size(Y_mat));

% U_mat = U_mat .* rmv_U_mat;
%% Forming AU matrix and calculating NOFRFs
for i = 1:len_adj_hlf
    AU_11 = A_mat.*repmat(U_mat_real(i,:),n_A,1);
    AU_12 = A_mat.*((-1).*repmat(U_mat_imag(i,:),n_A,1));
    AU_21 = A_mat.*repmat(U_mat_imag(i,:),n_A,1);
    AU_22 = A_mat.*repmat(U_mat_real(i,:),n_A,1);
    
    AU = [AU_11 , AU_12 ; AU_21 , AU_22];
    Y_arng = [Y_mat_real(i,:),Y_mat_imag(i,:)]';
    
    AU = AU(:,logical([rmv_U_mat(i,:),rmv_U_mat(i,:)])); %remove corresponding Un's of the NOFRFs that are not valid in the current freq
    
    G_mat = (AU'*AU)\(AU'*Y_arng);
    
    len_G_mat_hlf = length(G_mat)/2;
    
    G_mat_real([logical(rmv_U_mat(i,:))],i) = G_mat(1:len_G_mat_hlf);
    G_mat_imag([logical(rmv_U_mat(i,:))],i) = G_mat(len_G_mat_hlf+1:end);
end
G_LS_2 = complex(G_mat_real,G_mat_imag);
%% Test the evaluated NOFRFs

for i = 1:n_A
    AU_test = (U_mat .* repmat(A_mat(i,:),len_adj_hlf,1)).';
    Y_LS_2(:,i) = sum( ( G_LS_2 .* AU_test ) ,1).';
end

error = (Y_mat - Y_LS_2) .* repmat(freq_rng_all,1,n_A);
sse_norm = sum(abs(error).^2,1) ./ var(abs(Y_mat));

sse_norm_real = sum(abs(real(error)).^2,1) ./ var(abs(real(Y_mat)));
sse_norm_imag = sum(abs(imag(error)).^2,1) ./ var(abs(imag(Y_mat)));

if displ == 1
    disp('--------------------- NOFRF - M-LS ------------------------------');
    disp(['Norm. SSE Real      = ',num2str(sse_norm_real)]);
    disp(['Norm. SSE Imag      = ',num2str(sse_norm_imag)]);
    disp(['Norm. SSE           = ',num2str(sse_norm)]);
    disp('-----------------------------------------------------------------');
end
sse = [sse_norm, sse_norm_real, sse_norm_imag];
%disp('NOFRFs extracted');
end