clear; clc; close all;

%%
% Transmissibility and NOFRF analysis of the jump phenomenon.
% The example follows the low-amplitude local approximation shown in
% Figs. 5.9 and 5.10 of the MPhil thesis.

%% Paths
addpath('C:\Users\rajin\OneDrive - Coventry University\PhD project\GitHub\NL-FRA\NOFRF\');

%% Frequency range
frq_rng = linspace(1,35,200);
len_frq_rng = length(frq_rng);

Fs = 1000;
Ts = 1/Fs;
fftn = 10000;

% Simulate the transient before t = 0 and use the response after t = 0
% for the NOFRF and transmissibility calculations.
t_pre = 10;
t_record = fftn/Fs;
tspan_full = -t_pre:Ts:t_record;
[~,i_tz] = min(abs(tspan_full));
tspan = tspan_full(i_tz:end);
tspan = tspan(1:fftn+1);
len = length(tspan);

len_adj = fftn;
len_adj_hlf = floor(len_adj/2)+1;

%% NOFRF settings
nl_ord_set = [1:4];
N = max(nl_ord_set);

% NOFRFs are evaluated using the seven input amplitudes in equation (5.20).
A = [0.2,0.18333,0.16667,0.15,0.13333,0.11667,0.1]';
n_A = length(A);

% The evaluated NOFRFs are tested at A = 0.4 as in Fig. 5.5.
Amp = 0.4;

% SISO_NOFRF display options
% displ(1) - LS evaluation information
% displ(2) - NOFRF prediction/validation plots
% displ(3) - U_n, Y and valid-frequency plots
% displ(4) - NOFRF magnitude/phase plots
gc = 'b';
harm_inpt = 1;
lw = 0.5;
norm = 0;
displ = [0,0,0,0];

y0 = [0,0]';

%% Preallocate storage
G_LS2_cell = cell(1,len_frq_rng);

Trns_func = zeros(len_frq_rng,1);
Trns_func_3H = zeros(len_frq_rng,1);
Trns_func_NOFRF_LS2 = zeros(len_frq_rng,1);
Trns_func_NOFRF_LS2_3H = zeros(len_frq_rng,1);

% Output frequency response contribution from each nonlinear order.
Yn_fund = zeros(N,len_frq_rng);
Yn_3H = zeros(N,len_frq_rng);

%% Evaluate NOFRFs
for i = 1:len_frq_rng
    f1 = frq_rng(i);
    f2 = f1;

    u_full = sin(2*pi*f1*tspan_full);
    u = u_full(i_tz:end);
    u = u(1:len);

    % Generate the input-output data for NOFRF evaluation.
    ode_func = @(t,y,Amp_i) ODE_func(t,y,f1,Amp_i);
    Y_full = zeros(length(tspan_full),n_A);
    Y_full = NOFRF_data_sim(ode_func,tspan_full,n_A,A,Y_full,y0);
    Y = Y_full(i_tz:end,:);
    Y = Y(1:len,:);

    % Response used by SISO_NOFRF for NOFRF testing/validation.
    y_full = ode4(@(t,y) ode_func(t,y,A(1)),tspan_full,y0);
    y_test = y_full(i_tz:end,1);
    y_test = y_test(1:len);

    [~,~,~,~,~,~,~,G_LS_2,~,~] = ...
        SISO_NOFRF(Fs,Ts,tspan,fftn,f1,f2,u,A,A(1), ...
        nl_ord_set,gc,harm_inpt,lw,norm,displ,Y,y_test);

    G_LS2_cell{i} = G_LS_2;

    if mod(i,20) == 0 || i == len_frq_rng
        disp(['NOFRF evaluation: ',num2str(i),'/',num2str(len_frq_rng)]);
    end
end

%% Generate transmissibility using the evaluated NOFRFs
for i = 1:len_frq_rng
    f1 = frq_rng(i);
    f2 = f1;

    % Unit-amplitude reference input and the test input at Amp = 0.4.
    u_ref_full = sin(2*pi*f1*tspan_full);
    u_ref = u_ref_full(i_tz:end);
    u_ref = u_ref(1:len);
    u = Amp.*u_ref;

    % Actual nonlinear-system response.
    ode_func = @(t,y,Amp_i) ODE_func(t,y,f1,Amp_i);
    y_full = ode4(@(t,y) ode_func(t,y,Amp),tspan_full,y0);
    y = y_full(i_tz:end,1);
    y = y(1:len);

    y_fft = fft(y,len_adj).*Ts;
    u_fft = fft(u,len_adj).*Ts;
    Y_fft = y_fft(1:len_adj_hlf);
    U_fft = u_fft(1:len_adj_hlf);

    % Generate output response using the evaluated NOFRFs.
    G_LS_2 = G_LS2_cell{i};
    [Yn,Y_NOFRF,len_NOFRF] = ...
        nofrf_test(Amp,u_ref,G_LS_2,len_adj,len_adj_hlf, ...
        N,Fs,Ts,f1,f2);

    Y_NOFRF = [Y_NOFRF;zeros(len_adj_hlf-len_NOFRF,1)];
    Yn = [Yn;zeros(len_adj_hlf-len_NOFRF,N)];

    % Frequency indices for the excitation frequency and the 3rd harmonic.
    frq_ind = round((f1/Fs)*len_adj)+1;
    frq_ind3 = round((3*f1/Fs)*len_adj)+1;

    input_mag = abs(U_fft(frq_ind));

    %% Transmissibility at excitation frequency
    Trns_func(i) = abs(Y_fft(frq_ind))/input_mag;
    Trns_func_NOFRF_LS2(i) = abs(Y_NOFRF(frq_ind))/input_mag;

    %% Transmissibility at 3rd harmonic
    Trns_func_3H(i) = abs(Y_fft(frq_ind3))/input_mag;
    Trns_func_NOFRF_LS2_3H(i) = abs(Y_NOFRF(frq_ind3))/input_mag;

    %% Contribution from each nonlinear order
    for n = nl_ord_set
        Yn_fund(n,i) = abs(Yn(frq_ind,n));
        Yn_3H(n,i) = abs(Yn(frq_ind3,n));
    end
end

%% NMSE
NMSE_Trans = sum((Trns_func-Trns_func_NOFRF_LS2).^2)/var(Trns_func);
NMSE_Trans_3H = sum((Trns_func_3H-Trns_func_NOFRF_LS2_3H).^2)/var(Trns_func_3H);

disp(['NMSE_Trans    = ',num2str(NMSE_Trans)]);
disp(['NMSE_Trans_3H = ',num2str(NMSE_Trans_3H)]);

%% Actual and NOFRF generated transmissibility - Fig. 5.9
figure;

subplot(2,2,1);
plot(frq_rng,Trns_func,'.-');hold on;
plot(frq_rng,Trns_func_NOFRF_LS2,'o','MarkerSize',4);
grid on;
xlabel('Excitation Frequency (Hz)');
ylabel('|Trans|');
legend('Actual','NOFRFs estimate','Location','best');
xlim([0 35]);

subplot(2,2,2);
plot(3.*frq_rng,Trns_func_3H,'.-');hold on;
plot(3.*frq_rng,Trns_func_NOFRF_LS2_3H,'o','MarkerSize',4);
grid on;
xlabel('3 * Excitation Frequency (Hz)');
ylabel('|Trans_3|');
xlim([0 120]);

subplot(2,2,3);
semilogy(frq_rng,Trns_func,'.-');hold on;
semilogy(frq_rng,Trns_func_NOFRF_LS2,'o','MarkerSize',4);
grid on;
xlabel('Excitation Frequency (Hz)');
ylabel('|Trans| (log scale)');
xlim([0 35]);

subplot(2,2,4);
semilogy(3.*frq_rng,Trns_func_3H,'.-');hold on;
semilogy(3.*frq_rng,Trns_func_NOFRF_LS2_3H,'o','MarkerSize',4);
grid on;
xlabel('3 * Excitation Frequency (Hz)');
ylabel('|Trans_3| (log scale)');
xlim([0 120]);

%% Output frequency response contribution of each nonlinear order - Fig. 5.10
odd_fund = 1:2:N;
odd_3H = 3:2:N;

figure;

subplot(2,2,1);
plot(frq_rng,Yn_fund(odd_fund,:)','.-');
grid on;
xlabel('Excitation Frequency (Hz)');
ylabel('Magnitude of Y_n');
legend('1st order contribution','3rd order contribution', ...
    '5th order contribution','7th order contribution', ...
    '9th order contribution','Location','best');
xlim([0 35]);

subplot(2,2,2);
plot(3.*frq_rng,Yn_3H(odd_3H,:)','.-');
grid on;
xlabel('3 * Excitation Frequency (Hz)');
ylabel('Magnitude of Y_n');
legend('3rd order contribution','5th order contribution', ...
    '7th order contribution','9th order contribution', ...
    'Location','best');
xlim([0 120]);

subplot(2,2,3);
semilogy(frq_rng,Yn_fund(odd_fund,:)','.-');
grid on;
xlabel('Excitation Frequency (Hz)');
ylabel('Magnitude of Y_n (log scale)');
xlim([0 35]);

subplot(2,2,4);
semilogy(3.*frq_rng,Yn_3H(odd_3H,:)','.-');
grid on;
xlabel('3 * Excitation Frequency (Hz)');
ylabel('Magnitude of Y_n (log scale)');
xlim([0 120]);

%% Local Functions
function dy = ODE_func(t,Y,f1,Amp)

C = 3.84*pi;
K1 = (12*pi)^2;
K3 = 0.1*(12*pi)^6;

u = Amp.*sin(2*pi*f1*t);

dy = [Y(2);...
      u - C.*Y(2) - K1.*Y(1) - K3.*Y(1).^3];

end
