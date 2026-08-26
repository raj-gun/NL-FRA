%{1
clear all; clc; close all

% Wigren and Schoukens (2017), Coupled Electric Drives Data Set and
% Reference Models, Technical Report 2017-024, Uppsala University.
%
% The example identifies a polynomial NARX model from the uniformly
% distributed-input CED data set and evaluates the first four NOFRFs of
% the identified model using a general band-limited probing input.
%
% DATAUNIF.MAT contains two SISO records:
%   u11 -> z11 : first uniformly distributed-input realization
%   u12 -> z12 : second uniformly distributed-input realization
% The sampling period is Ts = 0.02 s (Fs = 50 Hz), and each record
% contains 500 input-output samples.

% Add NL-FRA and NonSysID functions using paths relative to this example file.
addpath('C:\Users\rajin\OneDrive - Coventry University\PhD project\GitHub\NonSysID\NonSysID\');
addpath('C:\Users\rajin\OneDrive - Coventry University\PhD project\GitHub\NL-FRA\NOFRF');

%% Import CED data and identify NARX model
% DATAUNIF.MAT is distributed with the Coupled Electric Drives benchmark:
% https://uu.diva-portal.org/smash/get/diva2:1165531/FULLTEXT01.zip
%
% Place DATAUNIF.MAT in Examples/Data before running this example.
load('C:\Users\rajin\OneDrive - Coventry University\PhD project\GitHub\NL-FRA\Examples\Data\Coupled_electric_drive\DATAUNIF.mat');


% Use the first uniformly distributed-input realization for identification.
u_ID = u11(:);
y_ID = z11(:);

% Keep the second realization as independent validation data.
u_val = u12(:);
y_val = z12(:);


mod_type = 'ARX'; % Model type ARX/AR
na1 = 1; na2 = 9; % Maximum and minimum output lags
nb1 = 1; nb2 = 7; % Maximum and minimum input lags
nl_ord_max = 3; % Maximum order of polynomial nonlinearity considered
is_bias = 0; % Specify if bias/DC offset is required, 0, or not, 1.
n_inpts = 1; % Specify number of inputs

x_iOFR = [false,false]; % Run more than one iteration of iOFR for [linear model, nonlinear model]

stp_cri = {'PRESS_min','PRESS_thresh'}; % Stopping criteria for [linear model, nonlinear model].
D1_thresh = [ 0 ,10^(-4.5)]; % PRESS_min for initial linear ARX model. I.e. automatic stopping, check NonSysID documentation.

KSA_h = 20; % Specify the number of steps for k-steps ahead prediction
RCT = 4; % Specify which RCT method to use, 1-4, 0 for no RCT.

sim = [1,1]; % Specify whether to simulate model and display results respectively
displ = 0; % Set to 1 to display all models generated from iOFRs, 0 otherwise

% Set 1 or 0 to use parallel processing to accelerate iOFRs,
% for [linear model, nonlinear model]
parall = [0,1];

[model,~,iOFR_table_lin,iOFR_table_nl,best_mod_ind_lin,best_mod_ind_nl,~] = ...
    NonSysID(mod_type,u_ID,y_ID,na1,na2,nb1,nb2,nl_ord_max,is_bias, ...
    n_inpts,KSA_h,RCT,x_iOFR,stp_cri,D1_thresh,displ,sim,parall);

% Print ARX/NARX models
disp('ARX model:');
disp(iOFR_table_lin{best_mod_ind_lin,1});
if best_mod_ind_nl ~= 0
    disp('NARX model:');
    tbl_NARX = join(iOFR_table_nl{best_mod_ind_nl,10}, ...
                    iOFR_table_nl{best_mod_ind_nl,1});
    disp(tbl_NARX);
end

%% Independent validation of the identified NARX model
[~,y_hat_val,error_val] = model_simulation(model,u_val,y_val,KSA_h);
disp(['Validation RMSE = ',num2str(sqrt(mean(error_val(:,1).^2)))]);

max_lag = max([model{2},model{4}]);
y_test = [zeros(max_lag,1); y_hat_val(:,1)];

figure;plot(y_val()); hold on; plot(y_test(:,1));

%}
%%
%{
% close all;
%% Sampling and input definition
Fs = 50;                 % Ts = 0.02 s in the CED benchmark
Ts = 1/Fs;
tspan =  (-511+0.5*Ts:512+0.5*Ts).*Ts; %[-511+0.5*Ts:Ts:511+0.5*Ts]; %
fftn = length(tspan)-1;

%% Generate data using the identified NARX model
% General band-limited probing input:
% u(t) = Amp_2*(3/(2*pi))*[sin(2*pi*10*t)-sin(2*pi*1*t)]/t
% The 1-10 Hz band lies below the 25 Hz Nyquist frequency of the CED data.
% Amp_2 = 0.05 gives a peak input magnitude of 1.35 V at t = 0.
f1 = 6.5;                 % upper edge of input band (Hz); f1 > f2 [3.5*1.833,3.5]
f2 = 3.5;                  % lower edge of input band (Hz)
u = sinc_difference_input(tspan,f1,f2)';
u = u .* (std(u_ID)/std(u)) .* 0.4;
len = length(u);

% Multiple constant input gains for LS-based NOFRF evaluation.
% The NOFRFs are invariant to a constant scaling of the probing input.
N = 4; % Order of nonlinearity to consider
n_A = N;
A = linspace(0.9,0.5,n_A)';
Y = zeros(len,n_A);

Y = NOFRF_data_sim_NARX(model,u,n_A,A,Y,KSA_h);

%% Independent response at the nominal input amplitude
Amp = 1;
y_test = simulate_NARX(model,Amp.*u,KSA_h);

len_adj = fftn;
len_adj_hlf = floor(len_adj/2)+1;
w = 0:Fs/len_adj:Fs/2;

disp(['Input amplitude factors = ',num2str(A')]);
disp('Simulations complete');

%% Evaluate NOFRFs up to fourth order
nl_ord_set = 1:N;
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

%}
%% Local Functions
function Y = NOFRF_data_sim_NARX(model,u,n_A,A,Y,KSA_h)
% Generate input-output data for LS-based NOFRF evaluation using the
% identified NARX model and the NonSysID model_simulation() function.

for i = 1:n_A
    Amp_1 = A(i);
    Y(:,i) = simulate_NARX(model,Amp_1.*u,KSA_h);
end
end

function y = simulate_NARX(model,u,KSA_h)
% Simulate the identified NARX model with zero initial conditions.
% model_simulation() removes the first max(na,nb) samples when forming
% the lagged information matrix; prepend them here to preserve the input
% record length required by SISO_NOFRF().

max_lag = max([model{2},model{4}]);
y0 = zeros(length(u),1);

[~,y_hat] = model_simulation(model,u,y0,KSA_h);
y = [zeros(max_lag,1); y_hat(:,1)];
end

function u = sinc_difference_input(t,f1,f2)
% General band-limited input, including its finite limit at t = 0.
u = (3/(2*pi)).*(sin(2*pi*f1.*t)-sin(2*pi*f2.*t))./t;
zero_idx = (t == 0);
u(zero_idx) = 3.*(f1-f2);
end
