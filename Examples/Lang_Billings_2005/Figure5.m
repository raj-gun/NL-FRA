clear all; clc; close all

% Lang and Billings (2005), International Journal of Control, 78(5),
% pp. 345-362. Reproduction example for Figure 5.
%
% The example uses the nonlinear oscillator in Eq. (35), the first input
% in Eq. (36), and evaluates the first four NOFRFs. Figure 5 in the paper
% displays |G_n(j2*pi*f)|, n = 1,...,4, over 10-20 Hz.

% Add NL-FRA functions using a path relative to this example file.
addpath('C:\Users\rajin\OneDrive - Coventry University\PhD project\GitHub\NL-FRA\NOFRF\');
addpath('C:\Users\rajin\OneDrive - Coventry University\PhD project\GitHub\NL-FRA\Examples');

%% Sampling and input definition from Eq. (36)
Fs = 200;               % Ts = 0.005 s in the paper
Ts = 1/Fs;
tspan = (-511:512).*Ts; % t = -511*0.005,...,512*0.005 s
fftn = length(tspan)-1;

%% Generate data using the nonlinear oscillator and RK4
% First input in Eq. (36):
% u(t) = (3/(2*pi)) * [sin(2*pi*55*t)-sin(2*pi*30*t)]/t
% sinc_difference_input() handles the removable singularity at t = 0.
f1 = 55;                % upper edge of input band (Hz); f1 > f2
f2 = 30;                % lower edge of input band (Hz)
Amp_2 = 1;
u = Amp_2 .* sinc_difference_input(tspan,f1,f2);
len = length(u);

% Multiple constant input gains for LS-based NOFRF evaluation.
% The NOFRFs are invariant to a constant scaling of the probing input.
A = [1, 0.85, 0.70, 0.55, 0.40]';
n_A = length(A);
Y = zeros(len,n_A);

ode_func = @(t,y,Amp) ODE_func(t,y,f1,f2,Amp,Amp_2);
y0 = [0,0]';
Y = NOFRF_data_sim(ode_func,tspan,n_A,A,Y,y0);

%% Independent response at the nominal input amplitude
Amp = 1;
y = ode4(@(t,y) ode_func(t,y,Amp),tspan,y0);
y_test = y(:,1);

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
function dy = ODE_func(t,Y,f1,f2,Amp,Amp_2)
% Nonlinear oscillator from Eq. (35) of Lang & Billings (2005):
% m*y'' + c*y' + k1*y + k2*y^2 + k3*y^3 = u(t)

u = Amp_2 .* sinc_difference_input(t,f1,f2);
u = Amp .* u;

m = 1;
c = 20;
k1 = 1e4;
k2 = 1e7;
k3 = 5e9;

dy = [Y(2); ...
      (1/m).*(u - c.*Y(2) - k1.*Y(1) - k2.*Y(1).^2 - k3.*Y(1).^3)];
end

function u = sinc_difference_input(t,f1,f2)
% Eq. (36), including its finite limit at t = 0.
u = (3/(2*pi)).*(sin(2*pi*f1.*t)-sin(2*pi*f2.*t))./t;
zero_idx = (t == 0);
u(zero_idx) = 3.*(f1-f2);
end
