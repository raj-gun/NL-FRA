# `nofrf_test`

`nofrf_test` reconstructs the individual nonlinear-order Output Frequency Response Functions (OFRFs) and the total output spectrum from previously evaluated NOFRFs for a general band-limited or sinusoidal probing input.

For a discrete multi-tone probing input, use `nofrf_test_comp`.

---

## Signatures

### General band-limited or sinusoidal input

```matlab
[Yn_NOFRF,Y_NOFRF,len_adj_hlf] = ...
    nofrf_test(Amp,u,G_LS_2,len_adj,len_adj_hlf,N,Fs,Ts,f1,f2);
```

### Discrete multi-tone input

```matlab
[Yn_NOFRF,Y_NOFRF,len_adj_hlf] = ...
    nofrf_test_comp(Amp,u,G_LS_2,len_adj,len_adj_hlf,N,Fs,Ts,pos_freq_comp);
```

---

## Parameters

| Name | Type | Required | Description |
|---|---|---:|---|
| `Amp` | `double` | Yes | Constant amplitude scaling applied to the reference probing input for the reconstructed response. |
| `u` | `vector` | Yes | Reference probing input used when the NOFRFs were evaluated. |
| `G_LS_2` | `matrix` | Yes | Evaluated NOFRFs. Columns correspond to nonlinear orders. |
| `len_adj` | `int` | Yes | FFT length used for the reconstruction. |
| `len_adj_hlf` | `int` | Yes | Requested number of non-negative frequency bins, normally `floor(len_adj/2)+1`. The routine may reduce this value to the maximum frequency supported by the nonlinear-order frequency mask. |
| `N` | `int` | Yes | Maximum nonlinear order represented in `G_LS_2`. |
| `Fs` | `double` | Yes | Sampling frequency in Hz. |
| `Ts` | `double` | Yes | Sampling period in seconds. |
| `f1` | `double` | `nofrf_test` only | Upper edge of the positive-frequency input band in Hz. For sinusoidal probing set `f1 = f2`. |
| `f2` | `double` | `nofrf_test` only | Lower edge of the positive-frequency input band in Hz. |
| `pos_freq_comp` | `vector` | `nofrf_test_comp` only | Positive discrete frequency components of the multi-tone probing input in Hz. |

---

## Returns

| Output | Type | Description |
|---|---|---|
| `Yn_NOFRF` | `matrix` | Reconstructed nonlinear-order output contributions. Column `n` contains `Y_n(jw) = G_n(jw)U_n(jw)` at the retained frequency bins. |
| `Y_NOFRF` | `vector` | Total reconstructed output spectrum obtained by summing `Yn_NOFRF` over nonlinear order. |
| `len_adj_hlf` | `int` | Number of retained non-negative frequency bins after applying the valid-frequency support. |

---

## Algorithm (High-Level)

1. Determine the valid frequency support of each nonlinear order using `nonlinear_freq_range` or `nonlinear_freq_comp`.
2. Construct the nonlinear input compositions from the FFTs of `u.^n` using

   ```matlab
   fft(u.^n,len_adj).*Ts.*((1/sqrt(n))/((2*pi)^(n-1)))
   ```

3. Apply the test-amplitude factor `Amp^n` to the `n`-th order input composition.
4. Multiply each valid input composition by the corresponding evaluated NOFRF.
5. Sum all nonlinear-order contributions to obtain the total reconstructed output spectrum.

In compact form, the routine evaluates

\[
Y_n(j\omega)=G_n(j\omega)\,A^nU_n(j\omega),
\qquad
Y(j\omega)=\sum_{n=1}^{N}Y_n(j\omega),
\]

within the valid frequency support of each nonlinear order.

---

## Example Usage

```matlab
len_adj = fftn;
len_adj_hlf = floor(len_adj/2)+1;
N = size(G_LS_2,2);

[Yn,Y_NOFRF,len_NOFRF] = ...
    nofrf_test(Amp,u,G_LS_2,len_adj,len_adj_hlf, ...
    N,Fs,Ts,f1,f2);
```

If the reconstructed spectrum is required over the full non-negative FFT range, zero-pad the returned arrays after `len_NOFRF`:

```matlab
Y_NOFRF = [Y_NOFRF; zeros(len_adj_hlf-len_NOFRF,1)];
Yn = [Yn; zeros(len_adj_hlf-len_NOFRF,N)];
```

This is the approach used in the transmissibility example.

---

## Transmissibility Example

For a sinusoidal input with excitation frequency `f1`, the repository evaluates the actual and NOFRF-generated transmissibilities at the excitation frequency and at the third harmonic. The corresponding FFT indices are

```matlab
frq_ind  = round((f1/Fs)*len_adj)+1;
frq_ind3 = round((3*f1/Fs)*len_adj)+1;
```

and the NOFRF-generated values are obtained from `Y_NOFRF`:

```matlab
input_mag = abs(U_fft(frq_ind));

Trans_NOFRF  = abs(Y_NOFRF(frq_ind))/input_mag;
Trans3_NOFRF = abs(Y_NOFRF(frq_ind3))/input_mag;
```

The columns of `Yn` can then be used to inspect the contribution of each nonlinear order at the fundamental and generated harmonic frequencies.

---

## Notes

- `nofrf_test` is intended for NOFRFs evaluated using `SISO_NOFRF`/`NOFRF_LS2`.
- `nofrf_test_comp` is intended for NOFRFs evaluated using `SISO_NOFRF_comp`/`NOFRF_LS2comp`.
- The reference input `u`, sampling settings and frequency definition should be consistent with those used during NOFRF evaluation.
- Frequencies outside the order-dependent valid support are excluded by the internal mask.

---

## See Also

- [`SISO_NOFRF`](./SISO_NOFRF.md)
- [`Example Code Structure`](./Example_code_structure.md)
- [`Transmisibility_NOFRF.m`](../Examples/Transmisibility_NOFRF.m)
