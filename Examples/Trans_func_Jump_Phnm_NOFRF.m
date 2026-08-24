clear; clc; close all;

%% Transmissibility and jump-phenomenon analysis using current NL-FRA
% This example is a cleaned and updated version of the legacy
% Trans_func_Jump_Phnm_NOFRF.m script. It evaluates NOFRFs for a nonlinear
% Duffing-type oscillator over a swept single-tone excitation frequency and
% compares the directly simulated transmissibility with the response
% reconstructed from the identified NOFRFs.
%
% The current SISO_NOFRF interface is used for NOFRF identification, as in
% Examples/TEST.m. Internal SISO_NOFRF plotting is disabled during the
% sweep using the four-element display flag:
%   displ = [LS_info, prediction_plots, diagnostic_plots, NOFRF_plots].
%
% The optional plots near the end reproduce the useful plotting ideas that
% were commented out in the legacy script. Set the corresponding flags to
% true to enable them.

%% Repository paths
this_file = mfilename('fullpath');
examples_dir = fileparts(this_file);
repo_root = fileparts(examples_dir);
addpath(fullfile(repo_root,'NOFRF'));
addpath(examples_dir);  % ode4.m

%% Analysis settings
frq_rng = linspace(1,100,200);   % excitation-frequency sweep (Hz)
len_frq_rng = numel(frq_rng);

Fs = 2000;                       % sampling frequency (Hz)
Ts = 1/Fs;                       % current NL-FRA convention
fftn = 20000;                    % FFT length -> 0.1 Hz bin spacing

% Simulate a transient interval before t = 0 and analyse the steady-state
% interval after t = 0. Ten seconds is long relative to the decay time of
% the oscillator used below and also gives fftn samples at Fs = 2000 Hz.
t_pre = 10;
t_record = fftn/Fs;
tspan_full = -t_pre:Ts:t_record;
[~,i_tz] = min(abs(tspan_full));
tspan = tspan_full(i_tz:end);

% Keep exactly fftn+1 time samples, matching the convention used elsewhere
% in the repository where fftn = length(tspan)-1.
tspan = tspan(1:fftn+1);
len = numel(tspan);
len_adj = fftn;
len_adj_hlf = floor(len_adj/2)+1;
w = 0:Fs/len_adj:Fs/2;

%% NOFRF settings
nl_ord_set = 1:5;
N = max(nl_ord_set);

% Input scaling factors used to identify the NOFRFs. These are retained
% from the legacy example and provide N distinct amplitude experiments.
A = linspace(1.1,1,N)';
n_A = numel(A);
Amp_test = 1;

% SISO_NOFRF options
gc = 'b';
harm_inpt = 1;                   % single-tone/harmonic probing input
lw = 0.5;
norm_plot = 0;
displ = [0,0,0,0];              % suppress all internal sweep plots

% Initial state for all simulations
y0 = [0,0]';

%% Optional plotting controls
plot_transmissibility_3d = true;
plot_transmissibility_curves = false;
plot_log_transmissibility_curves = false;
plot_order_contributions = false;
plot_nofrf_maps = false;
plot_nofrf_scatter3 = false;
plot_nofrf_log_surface = false;
plot_nofrf_waterfall = false;
plot_nofrf_ribbon = false;

% NOFRF order used by the optional map/waterfall/ribbon plots.
map_order = 3;

%% Preallocate storage
G_LS2_cell = cell(1,len_frq_rng);

Trns_func_mat = zeros(len_frq_rng,n_A);
Trns_func_mat_3H = zeros(len_frq_rng,n_A);
Trns_func_NOFRF_LS2_mat = zeros(len_frq_rng,n_A);
Trns_func_NOFRF_LS2_mat_3H = zeros(len_frq_rng,n_A);

% Order-wise contribution to the fundamental transmissibility.
% Dimensions: nonlinear order x excitation frequency x input amplitude.
Trns_cont = zeros(N,len_frq_rng,n_A);

%% Identify NOFRFs at each excitation frequency
% For every frequency, generate the multi-amplitude training data required
% by the LS method, then call the current SISO_NOFRF wrapper. The wrapper
% performs the same identification pipeline used by Examples/TEST.m.
for i = 1:len_frq_rng
    f1 = frq_rng(i);
    f2 = f1;                     % single-tone input

    % Unit-amplitude reference input used for U_n(jw).
    u_full = cos(2*pi*f1*tspan_full);
    u_ref = u_full(i_tz:end);
    u_ref = u_ref(1:len);

    % Generate identification outputs for all amplitude scaling factors.
    ode_func = @(t,y,Amp) nonlinear_oscillator(t,y,f1,Amp);
    Y_full = zeros(numel(tspan_full),n_A);
    Y_full = NOFRF_data_sim(ode_func,tspan_full,n_A,A,Y_full,y0);
    Y_train = Y_full(i_tz:end,:);
    Y_train = Y_train(1:len,:);

    % Independent nominal-amplitude response used by SISO_NOFRF for
    % reconstruction/validation. Only the steady-state region is retained.
    y_full = ode4(@(t,y) ode_func(t,y,Amp_test),tspan_full,y0);
    y_test = y_full(i_tz:end,1);
    y_test = y_test(1:len);

    [~,~,~,~,~,~,~,G_LS_2,~,~] = ...
        SISO_NOFRF(Fs,Ts,tspan,fftn,f1,f2,u_ref,A,Amp_test, ...
        nl_ord_set,gc,harm_inpt,lw,norm_plot,displ,Y_train,y_test);

    G_LS2_cell{i} = G_LS_2;

    if mod(i,20) == 0 || i == len_frq_rng
        fprintf('NOFRF identification: %d/%d frequencies complete.\n', ...
            i,len_frq_rng);
    end
end

%% Test the identified NOFRFs over the training-amplitude range
% The NOFRFs are identified once per excitation frequency. For each input
% amplitude, the actual nonlinear response and the NOFRF reconstruction are
% then compared at the fundamental and third-harmonic frequencies.
for amp_i = 1:n_A
    Amp = A(amp_i);

    for i = 1:len_frq_rng
        f1 = frq_rng(i);
        f2 = f1;

        % Reference input for the NOFRF composition and actual scaled input
        % used for the direct nonlinear simulation.
        u_ref_full = cos(2*pi*f1*tspan_full);
        u_ref = u_ref_full(i_tz:end);
        u_ref = u_ref(1:len);
        u = Amp.*u_ref;

        ode_func = @(t,y,Amp_local) nonlinear_oscillator(t,y,f1,Amp_local);
        y_full = ode4(@(t,y) ode_func(t,y,Amp),tspan_full,y0);
        y = y_full(i_tz:end,1);
        y = y(1:len);

        % Direct FFTs of the input and output.
        y_fft = fft(y,len_adj).*Ts;
        u_fft = fft(u,len_adj).*Ts;
        Y_fft = y_fft(1:len_adj_hlf);
        U_fft = u_fft(1:len_adj_hlf);

        % Reconstruct the output from the NOFRFs at this excitation
        % frequency. nofrf_test is the same reconstruction routine used by
        % NOFRF_Y_pred inside SISO_NOFRF.
        G_LS_2 = G_LS2_cell{i};
        [Yn,Y_NOFRF,len_NOFRF] = ...
            nofrf_test(Amp,u_ref,G_LS_2,len_adj,len_adj_hlf, ...
            N,Fs,Ts,f1,f2);

        Y_NOFRF = [Y_NOFRF; zeros(len_adj_hlf-len_NOFRF,1)];
        Yn = [Yn; zeros(len_adj_hlf-len_NOFRF,N)];

        % FFT indices corresponding to f1 and 3*f1.
        frq_ind = round((f1/Fs)*len_adj)+1;
        frq_ind3 = round((3*f1/Fs)*len_adj)+1;

        input_mag = abs(U_fft(frq_ind));
        if input_mag == 0
            Trns_func_mat(i,amp_i) = NaN;
            Trns_func_NOFRF_LS2_mat(i,amp_i) = NaN;
            Trns_func_mat_3H(i,amp_i) = NaN;
            Trns_func_NOFRF_LS2_mat_3H(i,amp_i) = NaN;
            Trns_cont(:,i,amp_i) = NaN;
            continue;
        end

        % Fundamental transmissibility.
        Trns_func_mat(i,amp_i) = abs(Y_fft(frq_ind))/input_mag;
        Trns_func_NOFRF_LS2_mat(i,amp_i) = ...
            abs(Y_NOFRF(frq_ind))/input_mag;

        % Third-harmonic output normalized by the fundamental input.
        if frq_ind3 <= len_adj_hlf
            Trns_func_mat_3H(i,amp_i) = abs(Y_fft(frq_ind3))/input_mag;
            Trns_func_NOFRF_LS2_mat_3H(i,amp_i) = ...
                abs(Y_NOFRF(frq_ind3))/input_mag;
        else
            Trns_func_mat_3H(i,amp_i) = NaN;
            Trns_func_NOFRF_LS2_mat_3H(i,amp_i) = NaN;
        end

        % Order-wise contribution at the fundamental. This replaces the
        % legacy commented block and works for every requested order.
        for n = nl_ord_set
            Trns_cont(n,i,amp_i) = abs(Yn(frq_ind,n))/input_mag;
        end
    end

    fprintf('Validation: amplitude %d/%d complete (Amp = %.3g).\n', ...
        amp_i,n_A,Amp);
end

%% Error metrics
NMSE_Trans_mat = column_nmse(Trns_func_mat,Trns_func_NOFRF_LS2_mat);
NMSE_Trans_3H_mat = column_nmse(Trns_func_mat_3H, ...
    Trns_func_NOFRF_LS2_mat_3H);

disp(['NMSE_Trans    = ',num2str(NMSE_Trans_mat)]);
disp(['NMSE_Trans_3H = ',num2str(NMSE_Trans_3H_mat)]);

%% Main transmissibility comparison
if plot_transmissibility_3d
    F = repmat(frq_rng',1,n_A);
    A_plot = repmat(A',len_frq_rng,1);

    figure;
    subplot(2,1,1);
    plot3(F,A_plot,Trns_func_mat,'Marker','o'); hold on;
    plot3(F,A_plot,Trns_func_NOFRF_LS2_mat, ...
        'Marker','.','MarkerSize',10,'LineStyle','none');
    hold off; grid on;
    xlabel('Excitation frequency (Hz)');
    ylabel('Input amplitude');
    zlabel('Fundamental transmissibility');
    title('Direct simulation and NOFRF reconstruction');

    subplot(2,1,2);
    plot3(3.*F,A_plot,Trns_func_mat_3H,'Marker','o'); hold on;
    plot3(3.*F,A_plot,Trns_func_NOFRF_LS2_mat_3H, ...
        'Marker','.','MarkerSize',10,'LineStyle','none');
    hold off; grid on;
    xlabel('Third-harmonic frequency (Hz)');
    ylabel('Input amplitude');
    zlabel('Third-harmonic transmissibility');
end

%% Optional 2-D transmissibility curves
% These replace the two commented transmissibility figures in the legacy
% script. One figure is produced for each input amplitude.
if plot_transmissibility_curves || plot_log_transmissibility_curves
    for amp_i = 1:n_A
        figure;
        subplot(2,1,1);
        plot(frq_rng,Trns_func_mat(:,amp_i),'*-'); hold on;
        plot(frq_rng,Trns_func_NOFRF_LS2_mat(:,amp_i),'o-');
        hold off; grid on;
        xlabel('Excitation frequency (Hz)');
        ylabel('Transmissibility');
        title(['Fundamental, Amp = ',num2str(A(amp_i))]);
        legend('Direct','NOFRF','Location','best');

        subplot(2,1,2);
        plot(3.*frq_rng,Trns_func_mat_3H(:,amp_i),'*-'); hold on;
        plot(3.*frq_rng,Trns_func_NOFRF_LS2_mat_3H(:,amp_i),'o-');
        hold off; grid on;
        xlabel('Third-harmonic frequency (Hz)');
        ylabel('Transmissibility');
        title('Third harmonic');
        legend('Direct','NOFRF','Location','best');

        if plot_log_transmissibility_curves
            subplot(2,1,1); set(gca,'YScale','log');
            subplot(2,1,2); set(gca,'YScale','log');
        end
    end
end

%% Optional order-wise transmissibility contributions
% This generalizes the legacy Y1/Y3/Y5 contribution plotting code. Select
% any amplitude column by changing contribution_amp_index.
if plot_order_contributions
    contribution_amp_index = n_A;
    figure;
    plot(frq_rng,squeeze(Trns_cont(nl_ord_set,:,contribution_amp_index))', ...
        'Marker','.');
    grid on;
    xlabel('Excitation frequency (Hz)');
    ylabel('|Y_n(f_1)| / |U(f_1)|');
    title(['Order-wise contribution, Amp = ', ...
        num2str(A(contribution_amp_index))]);
    legend(arrayfun(@(n) ['n = ',num2str(n)],nl_ord_set, ...
        'UniformOutput',false),'Location','best');
end

%% Assemble NOFRF magnitude maps across the excitation sweep
% Each excitation frequency can have a different valid NOFRF frequency
% extent. NaN padding keeps unavailable regions out of the plots.
if plot_nofrf_maps || plot_nofrf_scatter3 || plot_nofrf_log_surface || ...
        plot_nofrf_waterfall || plot_nofrf_ribbon
    max_nofrf_len = max(cellfun(@(G) size(G,1),G_LS2_cell));
    w_nofrf = (0:max_nofrf_len-1).*(Fs/len_adj);
    NOFRF_mag_freq_rng = cell(N,1);

    for n = 1:N
        NOFRF_mag_freq_rng{n} = NaN(max_nofrf_len,len_frq_rng);
        for i = 1:len_frq_rng
            G = G_LS2_cell{i};
            if n <= size(G,2)
                this_len = size(G,1);
                NOFRF_mag_freq_rng{n}(1:this_len,i) = abs(G(:,n));
            end
        end
    end

    G_map = NOFRF_mag_freq_rng{map_order};
    [F_grid,W_grid] = meshgrid(frq_rng,w_nofrf);

    if plot_nofrf_maps
        figure;
        surf(F_grid,W_grid,G_map,'EdgeColor','none');
        xlabel('Excitation frequency (Hz)');
        ylabel('NOFRF frequency (Hz)');
        zlabel(['|G_',num2str(map_order),'|']);
        title(['NOFRF magnitude map, order ',num2str(map_order)]);
        colorbar;
    end

    if plot_nofrf_scatter3
        valid_map = isfinite(G_map) & (G_map > 0);
        figure;
        plot3(F_grid(valid_map),W_grid(valid_map),G_map(valid_map),'.');
        grid on;
        xlabel('Excitation frequency (Hz)');
        ylabel('NOFRF frequency (Hz)');
        zlabel(['|G_',num2str(map_order),'|']);
        title(['NOFRF scatter map, order ',num2str(map_order)]);
    end

    if plot_nofrf_log_surface
        figure;
        surf(F_grid,W_grid,G_map,'EdgeColor','none');
        set(gca,'ZScale','log');
        xlabel('Excitation frequency (Hz)');
        ylabel('NOFRF frequency (Hz)');
        zlabel(['|G_',num2str(map_order),'|']);
        title(['Log-scale NOFRF magnitude map, order ',num2str(map_order)]);
        colorbar;
    end

    if plot_nofrf_waterfall
        figure;
        waterfall(F_grid',W_grid',G_map');
        xlabel('Excitation frequency (Hz)');
        ylabel('NOFRF frequency (Hz)');
        zlabel(['|G_',num2str(map_order),'|']);
        title(['NOFRF waterfall, order ',num2str(map_order)]);
        colorbar;
    end

    if plot_nofrf_ribbon
        figure;
        ribbon(F_grid',G_map');
        xlabel('Excitation-frequency index');
        ylabel('NOFRF frequency bin');
        zlabel(['|G_',num2str(map_order),'|']);
        title(['NOFRF ribbon, order ',num2str(map_order)]);
        colorbar;
    end
end

%% Local functions
function dy = nonlinear_oscillator(t,Y,f1,Amp)
% Nonlinear oscillator used by the legacy Trans_func_Jump_Phnm_NOFRF.m.
%
%   y'' + C*y' + K1*y + K3*y^3 = Amp*cos(2*pi*f1*t)
%
% with natural frequency w0 = 12*pi rad/s and damping ratio 0.04.

w0 = 12*pi;
C = 2*0.04*w0;
K1 = w0^2;
K3 = 0.1*w0^6;

u = Amp.*cos(2*pi*f1*t);

dy = [Y(2); ...
      u - C.*Y(2) - K1.*Y(1) - K3.*Y(1).^3];
end

function nmse = column_nmse(y,y_hat)
% Normalized mean-square error evaluated independently for each column.
% NaN values (for example, harmonics above Nyquist) are ignored.

n_col = size(y,2);
nmse = NaN(1,n_col);
for k = 1:n_col
    valid = isfinite(y(:,k)) & isfinite(y_hat(:,k));
    yk = y(valid,k);
    yhatk = y_hat(valid,k);
    if numel(yk) > 1 && var(yk) > 0
        nmse(k) = sum((yk-yhatk).^2)/var(yk);
    end
end
end
