# NPSAT URF Algorithm Notes and Review

This note describes the current `NPSATurf` implementation in `CPP/NPSAT_URF_main.h`, reviews potential implementation issues, and compares the current Eigen-only lognormal fitting routine against the older dlib-based version.

## What `NPSATurf` Does

`NPSATurf` receives a discretized streamline as a vector of `segInfo` records. Each segment stores a velocity `v` and length `l`. For a selected porosity or velocity multiplier, the function computes a one-dimensional unit response function (URF) along the streamline, then fits the resulting breakthrough curve increments to a lognormal distribution.

The high-level process is:

1. Compute physical parameters for the streamline.
   - `vel = segment_velocity / velMult`
   - `age += segment_length / vel`
   - `aL = alpha * total_streamline_length^beta`
   - `Del = aL * vel + Dm`

2. Assemble finite-element matrices over the streamline.
   - `Dglo` is the advection-dispersion operator.
   - `Mglo` is the mass matrix.
   - If `calcDecay` is enabled, `DgloDecay` adds first-order decay terms using `lambda = log(2)/(halfTime*365)`.
   - The streamline is traversed in reverse order, so the calculation works from the discharge end back through the streamline.

3. Build the time-stepping system.
   - The code uses a theta-method form:
     - `Aglo = Mglo + omega * TimeStep * Dglo`
     - `Bglo = Mglo - (1 - omega) * TimeStep * Dglo`
   - The first node is treated as a fixed concentration boundary with `Cprev(0) = 1.0`.
   - The reduced system solves only unknown nodes `1..nel`.

4. Factor the sparse matrix once.
   - `KK = Aglo.block(1,1,nel,nel)`
   - `GG = Aglo.block(1,0,nel,1)`
   - `SparseLU` factorization is reused at every time step because the matrix is constant.

5. March forward in time.
   - At each step, the solver advances concentration and stores `Cprev(nel)` as the cumulative breakthrough at the outlet.
   - The loop stops when outlet cumulative breakthrough exceeds `URFtol`.
   - If the loop exceeds `maxTotalTime`, fitted outputs are set to `-66`.

6. Convert cumulative breakthrough to an incremental URF.
   - The first value is `unitBTC[0]`.
   - Later values are differences: `unitBTC[i] - unitBTC[i-1]`.
   - The peak position is tracked and used as the initial location for lognormal fitting.

7. Fit lognormal curves.
   - The base URF is fitted directly.
   - If decay is enabled, the decayed response is normalized by `1/unitBTCDecay.back()`.
   - If diffusion/difference output is enabled, `unitBTC - unitBTCDecay` is normalized and fitted separately.

## Review Findings

### Likely Bug: Decay Matrix Assembly

In the middle-element branch of the matrix assembly, the decay matrix uses the non-decay diagonal term:

```cpp
DgloTriDecay.emplace_back(ii, ii, D11 + DdiagPrevDecay);
```

This probably should be:

```cpp
DgloTriDecay.emplace_back(ii, ii, D11dec + DdiagPrevDecay);
```

The first and last element branches use `D11dec`, so the middle branch is inconsistent.

### Time Unit Ambiguity

`halfTime()` returns a daily decay rate:

```cpp
log(2)/(halfTime*365)
```

`fp.Age` is also converted from days to years by dividing by `365`. However, `TimeStep` is commented as `// in years` in `my_structures.h`, while the sample option says `TimeStep 365 // years`. The numerical use in `NPSATurf` looks like `TimeStep` should be in days, not years.

This should be clarified because decay, age, and transport time stepping must use consistent units.

### `maxTotalTime` Is Used as an Iteration Count

The loop checks:

```cpp
cnt++;
if (cnt > opt.maxTotalTime)
```

This treats `maxTotalTime` as a number of time steps, not as physical time. If `maxTotalTime` is meant to be years or days, the check should probably involve `cnt * TimeStep`.

### No Guard Against Zero Length or Zero Velocity

Several calculations divide by segment length or velocity:

```cpp
age += Lel / vel;
Del / Lel;
```

If a segment has zero length or nonpositive velocity, the matrix and age can become invalid. The upstream segment builder tries to avoid tiny segments, but it would still be safer for `NPSATurf` to explicitly reject invalid `Lel` or `vel`.

### Potential Normalization Divide-by-Zero

Decay and difference scaling can divide by zero:

```cpp
scaleDecay = 1/unitBTCDecay.back();
scaleDiff = 1/(unitBTC.back() - unitBTCDecay.back());
```

If the terminal values are zero, equal, or numerically very close, fitted outputs can become infinite or unstable.

### `calcDiff` Depends on `calcDecay`

The option parser currently forces `calcDiff = false` when decay is disabled. That is reasonable because the difference curve uses `unitBTCDecay`, but the dependency is implicit in `NPSATurf`. It may be worth documenting in the option file.

### Skip-Age Sentinel Values Are Odd

For young streamlines, the function returns sentinel values `-d` where `d` is the first integer threshold that the age is below. If `Age <= 0`, this can set values to `-0.0`, which is indistinguishable from `0.0` in many outputs.

### Unused Variables and Stale Comments

Variables such as `ddx`, `ddy`, `ddz`, `rho_term`, `RHS1`, and `URFDiff` are unused. Comments mention terms such as `Kd` and `rho` that are currently disabled. Cleaning these would make the numerical intent easier to audit.

## Current `fitLgnrm` vs. dlib Version

The current implementation preserves the same model, residual, and analytic derivative as the dlib version. The main difference is only the source of `sqrt(2)`: the old version uses `dlib::sqrt_2`; the current version uses the local constant `sqrt2`.

The dlib version delegates optimization to:

```cpp
dlib::solve_least_squares_lm(
    dlib::objective_delta_stop_strategy(1e-7),
    residual,
    residual_derivative,
    DS,
    x);
```

The current implementation manually performs a small Levenberg-Marquardt-like loop:

```cpp
(J^T J + lambda I) step = -J^T r
candidate = x + step
```

Then it accepts the candidate if the objective improves, decreases `lambda`, or rejects it and increases `lambda`.

Important differences:

- The current version returns `false` for empty data or nonpositive `maxYpos`; the dlib version does not explicitly guard these before `log(maxYpos)`.
- The current version clamps the fitted lognormal standard deviation parameter to `1e-8`, which prevents negative or zero `b`. The dlib version could step into invalid values unless dlib's optimizer avoided them.
- The current version returns `true` after `maxIterations` even if convergence was not reached. This differs from dlib, where failure would usually surface as an exception or nonconverged behavior.
- The current stopping test is an absolute objective change, `abs(prevObjective - objective) < 1e-7`. dlib's `objective_delta_stop_strategy` is similar in spirit but not guaranteed to behave identically.
- The current implementation does not check whether `candidateObjective` is finite.
- The current implementation uses `LDLT` without checking decomposition status. For a near-singular `J^T J + lambda I`, this may silently produce poor steps.
- The current implementation performs one trial per outer iteration. A more typical LM implementation may retry immediately with a larger damping parameter until it finds an acceptable step.

Overall, the current replacement is a reasonable lightweight substitute for removing dlib, and it keeps the same fitting target. The highest-priority improvement would be to return `false` when the maximum iteration count is reached without satisfying the convergence criterion, and to check all candidate objectives and steps for finite values.

## Suggested Follow-Up Fixes

1. Change the middle decay assembly term from `D11` to `D11dec`.
2. Clarify and enforce time units for `TimeStep`, `maxTotalTime`, `halfTime`, and reported `Age`.
3. Add explicit guards for nonpositive segment length and velocity.
4. Add safe checks before `scaleDecay` and `scaleDiff` normalization.
5. Make `fitLgnrm` report nonconvergence instead of returning `true` after exhausting iterations.
6. Remove unused variables once the numerical form is confirmed.
