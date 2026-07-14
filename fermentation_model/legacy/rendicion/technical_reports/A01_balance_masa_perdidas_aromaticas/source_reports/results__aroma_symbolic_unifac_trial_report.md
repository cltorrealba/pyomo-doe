# Symbolic UNIFAC trial

This trial tested replacing the fitted log-linear partition surrogate with the full original-UNIFAC activity-coefficient equations embedded symbolically in Pyomo.

## Implemented option

- Fructose is treated as glucose-equivalent for the liquid composition.
- Liquid components are water, ethanol, glucose-equivalent sugar, and a trace aroma species.
- The activity coefficient is computed with original UNIFAC group contributions:
  - ethyl acetate
  - isoamyl acetate
  - ethyl octanoate
- Vapor pressure is computed directly from the literature correlations selected by `thermo`:
  - ethyl acetate: `WAGNER_MCGARRY`
  - isoamyl acetate: `ANTOINE_WEBBOOK`
  - ethyl octanoate: `LANDOLT`

The symbolic Pyomo mode is selected with:

```powershell
python fermentation_model\run_aroma_campaign_doe.py `
  --partition-implementation symbolic_unifac `
  --partition-mode water_ethanol_total_sugar_as_glucose
```

The previous fitted surrogate remains the default:

```powershell
python fermentation_model\run_aroma_campaign_doe.py `
  --partition-implementation surrogate
```

## Numerical validation

The symbolic UNIFAC expression was compared against the package-level `thermo.UNIFAC_gammas()` calculation at representative fermentation points.

Maximum relative difference in `K_lg`: approximately `1.7%`.

The nonzero difference comes mainly from using the raw vapor-pressure literature correlations directly in Pyomo, while `Chemical.Psat` may apply extrapolation behavior outside the tabulated validity range for low-temperature fermentation conditions.

## Pyomo performance test

Model construction with symbolic UNIFAC was successful:

- build time for one aroma candidate: approximately `1.3 s`
- variables: `1663`
- active constraints: `1550`
- `K_lg` expressions evaluated successfully

However, solve/FIM performance was not acceptable:

- one-candidate `compute_FIM(method="sequential")` with `aroma_only` did not finish within `10 min`
- one nominal simulation with symbolic UNIFAC did not finish within `5 min`
- the surrogate regression smoke test finished successfully in approximately `42 s`

## Decision

The symbolic UNIFAC implementation is retained in the code as an experimental option, but it is not recommended for the current DOE/FIM workflow.

The operational recommendation is to keep the fitted UNIFAC-derived surrogate for model-based DOE, and use the symbolic UNIFAC mode for offline validation or spot-checking partition predictions.
