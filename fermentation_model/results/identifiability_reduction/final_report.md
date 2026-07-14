# Identifiable fermentation model after deep selection

Generated: 2026-06-06

## Final estimated subset

The recommended structurally/practically identifiable effective-parameter model estimates:

- `mu0`
- `qN`
- `betaG0`
- `betaF0`
- `qEG`
- `qEF`
- `iG`

The fit uses batches `25026`, `25086`, `25150`, and `25170`.
The weighted Gaussian SSE objective is `6001.982781`.

This replaces the earlier strict six-parameter reduced model by freeing `iG`.
The earlier six-parameter model remains a conservative reference with WSSE `6010.753439`.

## Fixed parameters

The following effective parameters remain fixed at reference values:

- `sN`
- `qXG`
- `qXF`
- `sG`
- `sF`
- `iE`
- `Kd0`
- `m0`

`qXG` and `qXF` remain fixed because the unregularized free fit drives both to their lower bounds and the weighted-relative FIM has a near-null direction dominated by `qXF`/`qXG`.

`iG` is no longer fixed. In the deep run it satisfied all decision filters:

- no active estimated-parameter bound
- full-rank weighted-relative FIM/Q
- two-sided 95% profile-likelihood crossing
- weak-prior Bayesian/Laplace posterior dominated by data

## FIM and sensitivity

For the final seven-parameter model:

- FIM condition number at 10% relative perturbation: `8.381e3`
- minimum relative eigenvalue: `1.193e-4`
- perturbation solves: all optimal
- weakest direction: dominated by `iG`, then `qEF`, `betaF0`, `betaG0`, `qEG`

This is weaker than the six-parameter model, but still full rank and profile-identifiable.

The rejected unregularized nine-parameter model has:

- `qXG` and `qXF` at lower bounds
- FIM condition number about `4.845e7`
- a near-null direction dominated by `qXF` and `qXG`

That model is not accepted as structurally identifiable.

## Profile likelihood

All seven estimated parameters cross the 95% chi-square profile-likelihood threshold on both sides:

- `mu0`
- `qN`
- `betaG0`
- `betaF0`
- `qEG`
- `qEF`
- `iG`

The released `iG` profile is therefore practically identifiable on the current data under the effective formulation.

## Bayesian check

A Laplace posterior check was run in log-parameter space using weak lognormal priors centered at reference values with `log_sd = 1.0`.

All seven estimated parameters are data-dominated under this weak prior. For `iG`, the posterior variance over weak-prior variance is about `0.043`, so the posterior is informed mainly by the likelihood, not only by the prior.

## Regularized curve-fit sensitivity

The best regularized curve-fit candidate is:

- estimated parameters: `mu0`, `qN`, `betaG0`, `betaF0`, `qEG`, `qEF`, `qXG`, `qXF`, `iG`
- log-L2 penalized parameters: `qXG`, `qXF`, `iG`
- log-L2 weight: `10`
- direct WSSE: `5992.323835`
- direct WSSE + log-L2 penalty: `5998.054501`
- active estimated-parameter bounds: none

This candidate is useful for curve-fit sensitivity and visual comparison, but it is not promoted to the structural model because `qXG` and `qXF` do not pass the unregularized structural/practical identifiability filters.

## DOE implications

The current model is suitable as the calibration baseline for DOE with seven estimated parameters.

DOE should target the fixed weak directions:

- `qXG` and `qXF`: improve biomass/sugar-growth separation, especially through richer viable biomass observations and conditions separating glucose and fructose growth uptake.
- `sG`, `sF`, and `sN`: excite saturation behavior directly instead of absorbing it into `betaG0`, `betaF0`, and `mu0`.
- `iE`, `Kd0`, and `m0`: require experiments that isolate ethanol inhibition, death/decay, and maintenance dynamics.

Primary result files:

- `final_identifiable_parameter_decision.csv`
- `final_fixed_parameter_decision.csv`
- `final_identifiability_summary.txt`
- `theta_plus_iG_identifiable.csv`
- `physical_plus_iG_identifiable.csv`
- `theta_final_identifiable.csv`
- `physical_final_identifiable.csv`
- `profile_summary_plus_iG_identifiable.csv`
- `bayesian_laplace_plus_iG_identifiable.csv`
- `../deep_model_selection/deep_model_selection_report.md`
- `../fit_strategy_analysis/fit_strategy_summary.csv`
