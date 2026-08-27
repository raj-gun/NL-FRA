# Example Code Structure

This page gives a general MATLAB code structure for evaluating and validating SISO NOFRFs with `NL-FRA`. The workflow follows the examples included in the repository.

---

## 1. Add NL-FRA to the MATLAB Path

```matlab
repo_root = 'path_to_NL-FRA';
addpath(fullfile(repo_root,'NOFRF'));
```

If the input-output data are generated using another package or model implementation, add that code to the MATLAB path separately. For example, the NARX example uses the `NonSysID` package to identify and simulate a polynomial NARX model.

---

## 2. Define Sampling and FFT Settings

```matlab
Fs = 1000;
Ts = 1/Fs;
fftn = 10000;

tspan = (0:fftn).*Ts;
len_adj = fftn;
len_adj_hlf = floor(len_adj/2)+1;
```

The probing input, evaluation outputs and independent validation output should use consistent sampling and record lengths.

---

## 3. Define the Probing Input

### General band-limited input

For a positive-frequency band `[f2,f1]`, with `f1 >= f2`, define the probing input and use `SISO_NOFRF`:

```matlab
f1 = 6.5;
f2 = 3.5;

u = (3/(2*pi)).*(sin(2*pi*f1.*tspan)-sin(2*pi*f2.*tspan))./tspan;
u(tspan == 0) = 3*(f1-f2);
u = u(:);
```

The repository's data-driven NARX example uses this type of input after normalising its amplitude relative to the identification input.

### Single-tone input

For sinusoidal probing, set `f1 = f2` and define a unit/reference sinusoid:

```matlab
f1 = 7;
f2 = f1;
u = cos(2*pi*f1*tspan)';
```

This is the form used when evaluating local NOFRFs during the transmissibility frequency sweep.

### Discrete multi-tone input

For a finite set of positive input frequencies, use `SISO_NOFRF_comp`:

```matlab
pos_freq_comp = [5,7,8];

u = cos(2*pi*pos_freq_comp(1)*tspan) + ...
    cos(2*pi*pos_freq_comp(2)*tspan) + ...
    cos(2*pi*pos_freq_comp(3)*tspan);
u = u(:);
```

---

## 4. Choose the Nonlinear Orders and Evaluation Amplitudes

```matlab
nl_ord_set = 1:5;
N = max(nl_ord_set);

A = linspace(1.3,1.2,N)';
n_A = length(A);

Amp = 1.4;
```

`A` contains the constant input-amplitude scalings used to generate the data for LS evaluation. `Amp` is the amplitude used to generate an independent response for testing the evaluated NOFRFs.

The NOFRF method uses the same reference probing waveform at several constant amplitudes. For each scaling `A(i)`, the applied input is

```matlab
u_i = A(i).*u;
```

and the corresponding output becomes column `i` of `Y`.

---

## 5. Generate or Acquire the Evaluation Data

Preallocate one output column per amplitude scaling:

```matlab
Y = zeros(length(u),n_A);
```

Then obtain the system response for each scaled input:

```matlab
for i = 1:n_A
    u_i = A(i).*u;

    % Replace this line with the experimental system, ODE simulation,
    % identified NARX model, neural-network model, or other dynamic model.
    Y(:,i) = simulate_system(u_i);
end
```

The data-driven LS formulation only requires the input-output responses. The mechanism used to generate those responses is separate from the NOFRF estimation itself.

---

## 6. Generate an Independent Validation Response

Use the test amplitude `Amp` to generate a separate response:

```matlab
y_test = simulate_system(Amp.*u);
```

This response is compared with the spectrum reconstructed from the evaluated NOFRFs.

---

## 7. Configure Display Options

```matlab
gc = 'b';
harm_inpt = 0;
lw = 0.5;
norm = 0;
displ = [1,1,1,1];
```

The display vector controls:

```text
displ(1) - LS evaluation information
displ(2) - NOFRF prediction/validation plots
displ(3) - U_n, Y and valid-frequency plots
displ(4) - NOFRF magnitude/phase plots
```

For a repeated frequency-sweep calculation, use

```matlab
displ = [0,0,0,0];
```

to suppress plotting and printed evaluation output.

---

## 8. Evaluate General Band-Limited or Sinusoidal NOFRFs

```matlab
[Norm_SSE_abs_LS_2,Norm_SSE_arg_LS_2,Norm_SSE_LS_2, ...
 Fe_n_Yn,Fe_p_Yn,Y_NOFRF_MLS_n,Y_NOFRF_LS_2,G_LS_2, ...
 sse,Y_model] = ...
    SISO_NOFRF(Fs,Ts,tspan,fftn,f1,f2,u,A,Amp, ...
    nl_ord_set,gc,harm_inpt,lw,norm,displ,Y,y_test);
```

`G_LS_2` contains the evaluated NOFRFs. The validation stage also returns the NOFRF-generated output spectrum and individual nonlinear-order contributions.

---

## 9. Evaluate Discrete Multi-Tone NOFRFs

For a multi-tone probing input, replace `f1,f2` with the positive discrete frequency vector:

```matlab
harm_inpt = 1;

[Norm_SSE_abs_LS_2,Norm_SSE_arg_LS_2,Norm_SSE_LS_2, ...
 Fe_n_Yn,Fe_p_Yn,Y_NOFRF_MLS_n,Y_NOFRF_LS_2,G_LS_2, ...
 sse,Y_model] = ...
    SISO_NOFRF_comp(Fs,Ts,tspan,fftn,pos_freq_comp,u,A,Amp, ...
    nl_ord_set,gc,harm_inpt,lw,norm,displ,Y,y_test);
```

`SISO_NOFRF_comp` determines the nonlinear output-frequency components generated from combinations of the supplied input frequencies and only solves for NOFRFs where the corresponding nonlinear input composition is valid.

---

## 10. Reconstruct a Response Separately

After evaluating `G_LS_2`, a response can be reconstructed at another constant input amplitude without rerunning the complete validation pipeline.

### Band-limited or sinusoidal input

```matlab
[Yn,Y_NOFRF,len_NOFRF] = ...
    nofrf_test(Amp,u,G_LS_2,len_adj,len_adj_hlf, ...
    N,Fs,Ts,f1,f2);
```

### Multi-tone input

```matlab
[Yn,Y_NOFRF,len_NOFRF] = ...
    nofrf_test_comp(Amp,u,G_LS_2,len_adj,len_adj_hlf, ...
    N,Fs,Ts,pos_freq_comp);
```

`Yn(:,n)` is the reconstructed contribution of nonlinear order `n`, and `Y_NOFRF` is their sum.

---

## 11. Example: Nonlinear Transmissibility

For a sinusoidal frequency sweep, evaluate local NOFRFs at each excitation frequency and store `G_LS_2`. Then reconstruct the response at the desired test amplitude using `nofrf_test`.

For excitation frequency `f1`, the fundamental and third-harmonic FFT indices are

```matlab
frq_ind  = round((f1/Fs)*len_adj)+1;
frq_ind3 = round((3*f1/Fs)*len_adj)+1;
```

If `U_fft` is the test-input spectrum, the transmissibilities are

```matlab
input_mag = abs(U_fft(frq_ind));

Trans  = abs(Y_NOFRF(frq_ind))/input_mag;
Trans3 = abs(Y_NOFRF(frq_ind3))/input_mag;
```

The individual nonlinear-order contributions are available directly from `Yn`:

```matlab
Yn_fund = abs(Yn(frq_ind,:));
Yn_3H   = abs(Yn(frq_ind3,:));
```

See [`Examples/Transmisibility_NOFRF.m`](../Examples/Transmisibility_NOFRF.m) for the complete frequency-sweep implementation.

---

## Repository Examples

- [`NARX_band_limited_NOFRF.m`](../Examples/NARX_band_limited_NOFRF.m): identifies a NARX model from the coupled electric-drive data and evaluates the first four NOFRFs using a general band-limited probing input.
- [`Lang_Billings_2005_exmpl_5_1.m`](../Examples/Lang_Billings_2005_exmpl_5_1.m): general band-limited NOFRF example based on a nonlinear dynamic model.
- [`Multi_tone_NOFRF.m`](../Examples/Multi_tone_NOFRF.m): evaluates NOFRFs for a Duffing oscillator excited at 5, 7 and 8 Hz.
- [`Transmisibility_NOFRF.m`](../Examples/Transmisibility_NOFRF.m): evaluates local NOFRFs over a sinusoidal frequency sweep and reconstructs the fundamental and third-harmonic transmissibilities.

---

## See Also

- [`SISO_NOFRF`](./SISO_NOFRF.md)
- [`nofrf_test`](./nofrf_test.md)
