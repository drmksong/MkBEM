# MkBEM Benchmark Roadmap

This directory now contains a first smoke benchmark, `MkBEM_benchmark_smoke`, that runs `THM_TeST2.dat` and checks that the solver creates output for every time step without `nan` or `inf` values.

Recommended next benchmark cases for publication work:

1. Pure elastic cavity benchmark
   - Compare displacement and traction around a circular cavity with a known Lamé/Kirsch solution.
   - Run mesh refinements such as `N = 8, 16, 32, 64` and report error convergence.

2. Pure thermal diffusion benchmark
   - Use a heat source or thermal boundary condition with an analytical transient solution based on the exponential integral.
   - Compare temperature profiles at early, middle, and late times.

3. Pure hydraulic diffusion benchmark
   - Use a pressure or fluid-source case with a Theis-type analytical solution.
   - Compare pressure and flow-rate profiles over logarithmic time steps.

4. Coupled THM benchmark
   - Start from the uncoupled thermal and hydraulic cases, then enable coupling terms and compare the sign and magnitude of displacement, temperature, and pressure responses.
   - Record all material parameters and units in the input file and manuscript table.

5. Conservation checks
   - Integrate boundary heat/fluid fluxes and compare them with prescribed source totals.
   - Treat conservation residuals as regression-test metrics.

For each benchmark, keep the input data, reference values, output parser, tolerance, and expected convergence trend together so the result can be reproduced from a fresh build.
