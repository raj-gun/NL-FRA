clear; clc; close all;

%% NOFRF local approximation of Duffing jump phenomena
% This example reproduces the analysis style of Figs. 5.9 and 5.10 in:
% S. R. A. S. Gunawardena, "A Control Systems Perspective to Condition
% Monitoring and Fault Diagnosis", MPhil thesis, University of Sheffield,
% 2018, Chapter 5.
%
% The nonlinear oscillator is
%
%   y''(t) + C y'(t) + K1 y(t) + K3 y^3(t) = A cos(2*pi*f_F*t)
%
% with
%   C  = 0.96*pi,
%   K1 = (12*pi)^2,
%   K3 = 0.1*(12*pi)^6.
%
% NOFRFs are evaluated locally over the low-amplitude interval
%   A1 = [1.3, 1.2]
% using N = 9 nonlinear orders. They are then tested outside this interval
% at A_test = 1.5, as in Fig. 5.9 of the thesis.
%
% Only two figures are produced:
%   Figure 1 - actual vs NOFRF-generated transmissibility (Fig. 5.9 style)
%   Figure 2 - order-wise output contributions (Fig. 5.10 style)

%% Repository paths
this_file = mfilename('fullpath');
examples_dir = fileparts(this_file);
repo_root = fileparts(examples_dir);
addpath(fullfile(repo_root,'NOFRF'));
addpath(examples_dir);  % ode4.m

%% Frequency sweep and numerical settings
% Fig. 5.9 spans approximately 0-35 Hz at the fundamental; the right-hand
% plots use the corresponding third harmonic, 3*f_F.
frq_rng = linspace(1,35,200);     % excitation frequency f_F (Hz)
len_frq_rng = numel(frq_rng);

% Sampling must cover the highest NOFRF harmonic considered: 9*35 Hz.
Fs = 1000;                        % sampling frequency (Hz)
Ts = 1/Fs;
fftn = 10000;                     % 10 s record -> 0.1 Hz FFT spacing

% Simulate a transient before t = 0, then retain the steady-state record.
t_pre = 10;
t_record = fftn/Fs;
tspan_full = -t_pre:Ts:t_record;
[~,i_tz] = min(abs(tspan_full));
tspan = tspan_full(i_tz:end);
tspan = tspan(1:fftn+1);
len = numel(tspan);

len_adj = fftn;
len_adj_hlf = floor(len_adj/2)+1;

%% NOFRF settings from the thesis low-amplitude local approximation
nl_ord_set = 1:9;
N = max(nl_ord_set);

% Nine amplitudes over A1 = [1.3,1.2], matching Table 2 / Fig. 5.8.
A_eval = linspace(1.3,1.2,N)';
n_A = numel(A_eval);

% Fig. 5.9 tests the locally evaluated NOFRFs outside A1 at A_p = 1.5.
A_test = 1.5;

% Current SISO_NOFRF options. All internal plotting is disabled because
% this example produces only the two thesis-style figures below.
gc = 'b';
harm_inpt = 1;
lw = 0.5;
norm_plot = 0;
displ = [0,0,0,0];

y0 = [0,0]';

%% Storage
G_cell = cell(1,len_frq_rng);
Yn_cell = cell(1,len_frq_rng);

Trans_actual = NaN(len_frq_rng,1);
Trans_NOFRF = NaN(len_frq_rng,1);
Trans3_actual = NaN(len_frq_rng,1);
Trans3_NOFRF = NaN(len_frq_rng,1);

% Output-frequency contributions |Y_n| used for Fig. 5.10.
Yn_fund = NaN(N,len_frq_rng);
Yn_3H = NaN(N,len_frq_rng);

%% 1. Evaluate local NOFRFs over A1 = [1.3,1.2]
for i = 1:len_frq_rng
    f1 = frq_rng(i);
    f2 = f1;                      % harmonic/single-tone input

    % Unit-amplitude reference input. Constant amplitude scaling is handled
    % through A_eval by SISO_NOFRF/NOFRF_LS2.
    u_full = cos(2*pi*f1*tspan_full);
    u_ref = u_full(i_tz:end);
    u_ref = u_ref(1:len);

    % Generate the nine input-output experiments used for LS evaluation.
    ode_func = @(t,y,Amp) duffing_jump_oscillator(t,y,f1,Amp);
    Y_full = zeros(numel(tspan_full),n_A);
    Y_full = NOFRF_data_sim(ode_func,tspan_full,n_A,A_eval,Y_full,y0);
    Y_train = Y_full(i_tz:end,:);
    Y_train = Y_train(1:len,:);

    % SISO_NOFRF requires a validation response. Use the upper edge of the
    % local evaluation interval; internal validation plots are disabled.
    y_eval_full = ode4(@(t,y) ode_func(t,y,A_eval(1)),tspan_full,y0);
    y_eval = y_eval_full(i_tz:end,1);
    y_eval = y_eval(1:len);

    [~,~,~,~,~,~,~,G_LS_2,~,~] = ...
        SISO_NOFRF(Fs,Ts,tspan,fftn,f1,f2,u_ref,A_eval,A_eval(1), ...
        nl_ord_set,gc,harm_inpt,lw,norm_plot,displ,Y_train,y_eval);

    G_cell{i} = G_LS_2;

    if mod(i,20) == 0 || i == len_frq_rng
        fprintf('NOFRF evaluation: %d/%d excitation frequencies complete.\n', ...
            i,len_frq_rng);
    end
end

%% 2. Test the local NOFRFs at A_test = 1.5
for i = 1:len_frq_rng
    f1 = frq_rng(i);
    f2 = f1;

    u_ref_full = cos(2*pi*f1*tspan_full);
    u_ref = u_ref_full(i_tz:end);
    u_ref = u_ref(1:len);
    u_test = A_test.*u_ref;

    % Direct nonlinear-system response at A_test.
    ode_func = @(t,y,Amp) duffing_jump_oscillator(t,y,f1,Amp);
    y_full = ode4(@(t,y) ode_func(t,y,A_test),tspan_full,y0);
    y_test = y_full(i_tz:end,1);
    y_test = y_test(1:len);

    Y_fft_full = fft(y_test,len_adj).*Ts;
    U_fft_full = fft(u_test,len_adj).*Ts;
    Y_fft = Y_fft_full(1:len_adj_hlf);
    U_fft = U_fft_full(1:len_adj_hlf);

    % Reconstruct the output using the locally evaluated NOFRFs.
    G_LS_2 = G_cell{i};
    [Yn,Y_NOFRF,len_NOFRF] = ...
        nofrf_test(A_test,u_ref,G_LS_2,len_adj,len_adj_hlf, ...
        N,Fs,Ts,f1,f2);

    Y_NOFRF = [Y_NOFRF; zeros(len_adj_hlf-len_NOFRF,1)];
    Yn = [Yn; zeros(len_adj_hlf-len_NOFRF,N)];
    Yn_cell{i} = Yn;

    % FFT bins for the excitation frequency and its third harmonic.
    ind1 = round((f1/Fs)*len_adj)+1;
    ind3 = round((3*f1/Fs)*len_adj)+1;

    input_mag = abs(U_fft(ind1));
    if input_mag == 0
        continue;
    end

    % Fig. 5.9 quantities: transmissibility at f_F and 3*f_F.
    Trans_actual(i) = abs(Y_fft(ind1))/input_mag;
    Trans_NOFRF(i) = abs(Y_NOFRF(ind1))/input_mag;

    if ind3 <= len_adj_hlf
        Trans3_actual(i) = abs(Y_fft(ind3))/input_mag;
        Trans3_NOFRF(i) = abs(Y_NOFRF(ind3))/input_mag;
    end

    % Fig. 5.10 quantities: individual nonlinear-order output responses.
    for n = nl_ord_set
        Yn_fund(n,i) = abs(Yn(ind1,n));
        if ind3 <= size(Yn,1)
            Yn_3H(n,i) = abs(Yn(ind3,n));
        end
    end
end

%% Error measures reported in the thesis discussion
NMSE_Trans = nmse_valid(Trans_actual,Trans_NOFRF);
NMSE_Trans3 = nmse_valid(Trans3_actual,Trans3_NOFRF);

fprintf('\nTest amplitude A_test = %.3f\n',A_test);
fprintf('NMSE fundamental transmissibility = %.6e\n',NMSE_Trans);
fprintf('NMSE third-harmonic transmissibility = %.6e\n',NMSE_Trans3);

%% Figure 1 - thesis Fig. 5.9 style
% Top row: linear magnitude.
% Bottom row: logarithmic magnitude.
% Left: fundamental transmissibility versus excitation frequency f_F.
% Right: third-harmonic transmissibility versus 3*f_F.
figure;

subplot(2,2,1);
plot(frq_rng,Trans_actual,'.-'); hold on;
plot(frq_rng,Trans_NOFRF,'o','MarkerSize',4);
hold off; grid on;
xlabel('Excitation frequency (Hz)');
ylabel('|Trans(f_F)|');
title('Fundamental transmissibility');
legend('Actual','NOFRF estimate','Location','best');
xlim([0 35]);

subplot(2,2,2);
plot(3*frq_rng,Trans3_actual,'.-'); hold on;
plot(3*frq_rng,Trans3_NOFRF,'o','MarkerSize',4);
hold off; grid on;
xlabel('3rd-harmonic frequency (Hz)');
ylabel('|Trans_3(3f_F)|');
title('Third-harmonic transmissibility');
xlim([0 120]);

subplot(2,2,3);
semilogy(frq_rng,Trans_actual,'.-'); hold on;
semilogy(frq_rng,Trans_NOFRF,'o','MarkerSize',4);
hold off; grid on;
xlabel('Excitation frequency (Hz)');
ylabel('|Trans(f_F)|');
xlim([0 35]);

subplot(2,2,4);
semilogy(3*frq_rng,Trans3_actual,'.-'); hold on;
semilogy(3*frq_rng,Trans3_NOFRF,'o','MarkerSize',4);
hold off; grid on;
xlabel('3rd-harmonic frequency (Hz)');
ylabel('|Trans_3(3f_F)|');
xlim([0 120]);

sgtitle(sprintf('Actual and NOFRF-generated transmissibility, A = %.1f',A_test));

%% Figure 2 - thesis Fig. 5.10 style
% For the cubic Duffing oscillator under a harmonic input, the relevant
% contributions are odd nonlinear orders. At the fundamental we display
% n = 1,3,5,7,9. At 3*f_F, the first possible contributing odd order is 3,
% so n = 3,5,7,9 are displayed as in the thesis discussion.
odd_fund = 1:2:N;
odd_3H = 3:2:N;

legend_fund = arrayfun(@(n) sprintf('%d%s order contribution', ...
    n,ordinal_suffix(n)),odd_fund,'UniformOutput',false);
legend_3H = arrayfun(@(n) sprintf('%d%s order contribution', ...
    n,ordinal_suffix(n)),odd_3H,'UniformOutput',false);

figure;

subplot(2,2,1);
plot(frq_rng,Yn_fund(odd_fund,:)','.-');
grid on;
xlabel('Excitation frequency (Hz)');
ylabel('|Y_n(f_F)|');
title('Contributions at f_F');
legend(legend_fund,'Location','best');
xlim([0 35]);

subplot(2,2,2);
plot(3*frq_rng,Yn_3H(odd_3H,:)','.-');
grid on;
xlabel('3rd-harmonic frequency (Hz)');
ylabel('|Y_n(3f_F)|');
title('Contributions at 3f_F');
legend(legend_3H,'Location','best');
xlim([0 120]);

subplot(2,2,3);
semilogy(frq_rng,Yn_fund(odd_fund,:)','.-');
grid on;
xlabel('Excitation frequency (Hz)');
ylabel('|Y_n(f_F)|');
xlim([0 35]);

subplot(2,2,4);
semilogy(3*frq_rng,Yn_3H(odd_3H,:)','.-');
grid on;
xlabel('3rd-harmonic frequency (Hz)');
ylabel('|Y_n(3f_F)|');
xlim([0 120]);

sgtitle(sprintf('Order-wise NOFRF output contributions, A = %.1f',A_test));

%% Local functions
function dy = duffing_jump_oscillator(t,Y,f_exc,Amp)
% Duffing oscillator from Eq. (5.21) of the thesis:
%   y'' + C*y' + K1*y + K3*y^3 = A*cos(2*pi*f_F*t)

C = 0.96*pi;
K1 = (12*pi)^2;
K3 = 0.1*(12*pi)^6;

u = Amp*cos(2*pi*f_exc*t);

dy = [Y(2); ...
      u - C*Y(2) - K1*Y(1) - K3*Y(1)^3];
end

function value = nmse_valid(actual,estimate)
% NMSE using finite points only.
valid = isfinite(actual) & isfinite(estimate);
a = actual(valid);
e = estimate(valid);

if isempty(a) || var(a) == 0
    value = NaN;
else
    value = sum((a-e).^2)/var(a);
end
end

function s = ordinal_suffix(n)
% English ordinal suffix for plot legends.
if mod(n,100) >= 11 && mod(n,100) <= 13
    s = 'th';
elseif mod(n,10) == 1
    s = 'st';
elseif mod(n,10) == 2
    s = 'nd';
elseif mod(n,10) == 3
    s = 'rd';
else
    s = 'th';
end
end
