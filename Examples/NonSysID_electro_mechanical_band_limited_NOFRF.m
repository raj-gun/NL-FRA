clear all; clc; close all

% NonSysID (2025), Journal of Open Source Software, 10(114), 8028.
% Real-data NARX model identified from the electro-mechanical system example.
%
% The example identifies a polynomial NARX model using the real input-output
% data in NonSysID, excites the identified model with a general band-limited
% input, and evaluates the first four NOFRFs.

% Add NL-FRA and NonSysID functions using paths relative to this example file.
addpath('C:\Users\rajin\OneDrive - Coventry University\PhD project\GitHub\NonSysID\NonSysID\');
addpath('C:\Users\rajin\OneDrive - Coventry University\PhD project\GitHub\NL-FRA\NOFRF');

%% Identify NARX model from the electro-mechanical system data
u_ID_dat = readmatrix('C:\Users\rajin\OneDrive - Coventry University\PhD project\GitHub\NonSysID\Examples\Electro-mechanical_system\Data\x_cc.csv'); 
y_ID_dat = readmatrix('C:\Users\rajin\OneDrive - Coventry University\PhD project\GitHub\NonSysID\Examples\Electro-mechanical_system\Data\y_cc.csv');

% Down sample data as in the NonSysID electro-mechanical system example.
dwn_smpl = 100;
u_ID_dat = u_ID_dat(1:dwn_smpl:end);
y_ID_dat = y_ID_dat(1:dwn_smpl:end);

tt_splt = 100:350;
u_ID = u_ID_dat(tt_splt);
y_ID = y_ID_dat(tt_splt);

% Model type ARX/AR
mod_type = 'ARX';

% Maximum and minimum output lags
na1 = 1; na2 = 3;

% Maximum and minimum input lags
nb1 = 1; nb2 = 3;

% Maximum order of polynomial nonlinearity considered
nl_ord_max = 2;

% Run more than one iteration of iOFR for [linear model, nonlinear model]
x_iOFR = [true,true];

% Stopping criteria for [linear model, nonlinear model].
stp_cri = {'PRESS_thresh','PRESS_thresh'};
D1_thresh = [10^(-10),10^(0.9)];

% Specify if bias/DC offset is required, 0, or not, 1.
is_bias = 0;

% Specify number of inputs
n_inpts = 1;

% Specify the number of steps for k-steps ahead prediction
KSA_h = 20;

% Specify which RCT method to use, 1-4, 0 for no RCT.
RCT = 3;

% Do not generate the NonSysID simulation plots in this example.
sim = [1,1];

% Set to 1 to display all models generated from iOFRs, 0 otherwise.
displ = 0;

% Set 1 or 0 to use parallel processing for [linear model, nonlinear model].
parall = [0,0];

[model,~,iOFR_table_lin,iOFR_table_nl,best_mod_ind_lin,best_mod_ind_nl,~] = ...
    NonSysID(mod_type,u_ID,y_ID,na1,na2,nb1,nb2,nl_ord_max,is_bias, ...
    n_inpts,KSA_h,RCT,x_iOFR,stp_cri,D1_thresh,displ,sim,parall);

disp('ARX model:');
disp(iOFR_table_lin{best_mod_ind_lin,1});

if best_mod_ind_nl ~= 0
    disp('NARX model:');
    tbl_NARX = join(iOFR_table_nl{best_mod_ind_nl,10}, ...
                    iOFR_table_nl{best_mod_ind_nl,1});
    disp(tbl_NARX);
end

%% Sampling and input definition
Fs = 200;
Ts = 1/Fs;
tspan = (-511:512).*Ts;
fftn = length(tspan)-1;

%% Generate data using the identified NARX model
% General band-limited probing input:
% u(t) = Amp_2*(3/(2*pi))*[sin(2*pi*55*t)-sin(2*pi*30*t)]/t
% Amp_2 = 1/15 gives a finite value of 5 at t = 0, matching the upper
% input level present in the electro-mechanical identification data.
f1 = 3.5;
f2 = 2;
Amp_2 = 1;
u = Amp_2 .* sinc_difference_input(tspan,f1,f2);
u = u(:);
u = u .* (std(u_ID)/std(u));
len = length(u);

% Multiple constant input gains for LS-based NOFRF evaluation.
% The NOFRFs are invariant to a constant scaling of the probing input.
A = [1, 0.85, 0.70, 0.55, 0.40]';
n_A = length(A);
Y = zeros(len,n_A);

Y = NOFRF_data_sim_NARX(model,u,n_A,A,Y,KSA_h);

%% Independent response at the nominal input amplitude
Amp = 1;
[~,y_hat] = model_simulation(model,Amp.*u,zeros(len,1),KSA_h);
max_lag = max([model{2},model{4}]);
y_test = [zeros(max_lag,1); y_hat(:,1)];

len_adj = fftn;
len_adj_hlf = floor(len_adj/2)+1;
w = 0:Fs/len_adj:Fs/2;

disp(['Input amplitude factors = ',num2str(A')]);
disp('Simulations complete');

%% Evaluate NOFRFs up to fourth order
nl_ord_set = 1:4;
N = max(nl_ord_set);

gc = 'b';
harm_inpt = 0;
lw = 0.5;
displ = [1,1,1];
norm = 0;

u_nofrf = u;

[Norm_SSE_abs_LS_2,Norm_SSE_arg_LS_2,Norm_SSE_LS_2, ...
 Fe_n_Yn,Fe_p_Yn,Y_NOFRF_MLS_n,Y_NOFRF_LS_2,G_LS_2, ...
 sse,Y_model] = ...
    SISO_NOFRF(Fs,Ts,tspan,fftn,f1,f2,u_nofrf,A,Amp, ...
    nl_ord_set,gc,harm_inpt,lw,norm,displ,Y,y_test);

%% Local Functions
function Y = NOFRF_data_sim_NARX(model,u,n_A,A,Y,KSA_h)
% Generate input-output data for LS-based NOFRF evaluation by repeatedly
% simulating an identified NonSysID model with scaled probing inputs.

max_lag = max([model{2},model{4}]);
y0 = zeros(length(u),1);

for i = 1:n_A
    Amp_1 = A(i);
    [~,y_hat] = model_simulation(model,Amp_1.*u,y0,KSA_h);
    Y(:,i) = [zeros(max_lag,1); y_hat(:,1)];
end
end

function u = sinc_difference_input(t,f1,f2)
% General band-limited input, including its finite limit at t = 0.
u = (3/(2*pi)).*(sin(2*pi*f1.*t)-sin(2*pi*f2.*t))./t;
zero_idx = (t == 0);
u(zero_idx) = 3.*(f1-f2);
end
