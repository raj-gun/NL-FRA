# NL-FRA Package Documentation

This provides documentation for nonlinear frequency response analysis using the `NL-FRA` package.  
The workflow includes evaluating NOFRFs from input-output data, reconstructing the output response, validating the evaluated NOFRFs and analysing nonlinear-order contributions.

- [**SISO_NOFRF**](./SISO_NOFRF.md)  
  Main functions for evaluating and validating SISO NOFRFs using general band-limited/sinusoidal inputs or discrete multi-tone inputs.  
  Describes function signatures, parameters, outputs, display controls and example usage.

- [**nofrf_test**](./nofrf_test.md)  
  Functions for reconstructing the individual OFRF contributions and the total output spectrum from evaluated NOFRFs.  
  Covers both band-limited/sinusoidal and discrete multi-tone probing inputs.

- [**Example Code Structure**](./Example_code_structure.md)  
  A step-by-step guide to applying the `NL-FRA` package in MATLAB.  
  Shows how to prepare probing inputs, generate or acquire scaled input-output data, evaluate NOFRFs and validate the reconstructed response.

---

## Getting Started

1. Prepare a probing input `u` and choose the nonlinear orders to evaluate.  
2. Obtain output responses for several constant input-amplitude scalings `A`.  
3. Use [`SISO_NOFRF`](./SISO_NOFRF.md) or `SISO_NOFRF_comp` to evaluate and validate the NOFRFs.  
4. Use [`nofrf_test`](./nofrf_test.md) or `nofrf_test_comp` when a separate NOFRF-generated response or nonlinear-order decomposition is required.  
5. Follow the workflow in the [Example Code Structure](./Example_code_structure.md) for a complete nonlinear frequency-response analysis pipeline.  
