clear all; clc; close all

% Multi-tone NOFRF example using the Duffing oscillator.
% The NOFRFs are evaluated using three discrete input frequencies at
% 5,7 and 8 Hz.

% Add NL-FRA functions using a path relative to this example file.
addpath('C:\Users\rajin\OneDrive - Coventry University\PhD project\GitHub\NL-FRA\NOFRF\');

%% Sampling and input definition
Fs = 1000;
Ts = 1/Fs;
fftn = 10000;

% Simulate the transient before t = 0 and use the response after t = 0
% for the NOFRF calculations.
t_pre = 10;
t_record = fftn/Fs;
tspan_full = -t_pre:Ts:t_record;
[~,i_tz] = min(abs(tspan_full));
tspan = tspan_full(i_tz:end);
tspan = tspan(1:fftn+1);
len = length(tspan);

pos_freq_comp = [5,7,8];

u_full = cos(2*pi*pos_freq_comp(1)*tspan_full) + ...
         cos(2*pi*pos_freq_comp(2)*tspan_full) + ...
         cos(2*pi*pos_freq_comp(3)*tspan_full);
u = u_full(i_tz:end);
u = u(1:len);

%% Generate data using the nonlinear oscillator and RK4
nl_ord_set = 1:5;
N = max(nl_ord_set);

% NOFRFs are evaluated over A1 = [1.3,1.2] using nine amplitudes.
A = linspace(1.3,1.2,N)';
n_A = length(A);
Y_full = zeros(length(tspan_full),n_A);

ode_func = @(t,y,Amp) ODE_func(t,y,pos_freq_comp,Amp);
y0 = [0,0]';
Y_full = NOFRF_data_sim(ode_func,tspan_full,n_A,A,Y_full,y0);
Y = Y_full(i_tz:end,:);
Y = Y(1:len,:);

%% Independent response at the test input amplitude
Amp = 1.4;
y_full = ode4(@(t,y) ode_func(t,y,Amp),tspan_full,y0);
y = y_full(i_tz:end,1);
y_test = y(1:len);

len_adj = fftn;
len_adj_hlf = floor(len_adj/2)+1;
w = 0:Fs/len_adj:Fs/2;

disp(['Input frequencies = ',num2str(pos_freq_comp),' Hz']);
disp(['Input amplitude factors = ',num2str(A')]);
disp(['Test input amplitude = ',num2str(Amp)]);
disp('Simulations complete');

%% Evaluate NOFRFs up to ninth order
gc = 'b';
harm_inpt = 1;
lw = 0.5;
displ = [1,1,1,1];
norm = 0;

u_nofrf = u;

[Norm_SSE_abs_LS_2,Norm_SSE_arg_LS_2,Norm_SSE_LS_2, ...
 Fe_n_Yn,Fe_p_Yn,Y_NOFRF_MLS_n,Y_NOFRF_LS_2,G_LS_2, ...
 sse,Y_model] = ...
    SISO_NOFRF_comp(Fs,Ts,tspan,fftn,pos_freq_comp,u_nofrf,A,Amp, ...
    nl_ord_set,gc,harm_inpt,lw,norm,displ,Y,y_test);

%% Local Functions
function dy = ODE_func(t,Y,pos_freq_comp,Amp)

w0 = 12*pi;
C = 2*0.04*w0;
K1 = w0^2;
K3 = 0.1*w0^6;

u = Amp.*(cos(2*pi*pos_freq_comp(1)*t) + ...
           cos(2*pi*pos_freq_comp(2)*t) + ...
           cos(2*pi*pos_freq_comp(3)*t));

dy = [Y(2);...
      u - C.*Y(2) - K1.*Y(1) - K3.*Y(1).^3];

end
