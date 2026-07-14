# MPCC parameter homologation from calibrated fermentation model

Date: 2026-06-29

Purpose: hand-off document for the MPCC agent. This file translates the calibrated
pilot-2025 fermentation model parameters into the naming and unit contract used
by the current `DC_dFVB_2026` MPCC implementation.

This document was prepared after checking:

- `C:/Users/ctorrealba/Downloads/estructura_modelo_mpcc_actual_para_homologacion_20260629.md`
- `fermentation_model/pilot_2025/pilot_2025_co2_solubility_integrated_doe.executed.ipynb`
- `fermentation_model/pilot_2025/run_pilot_2025_co2_solubility_integrated_doe.py`
- `fermentation_model/pilot_2025/results/co2_solubility_integrated_doe/theta_selected_integrated.csv`
- `fermentation_model/pilot_2025/results/co2_solubility_integrated_doe/co2_params_selected.csv`
- `fermentation_model/pilot_2025/results/co2_solubility_integrated_doe/final_fit_metrics.csv`
- `fermentation_model/pilot_2025/results/co2_solubility_integrated_doe/co2_o2_parameter_estimability_current.csv`
- `fermentation_model/pilot_2025/results/co2_solubility_integrated_doe/co2_o2_selected_campaign_hybrid.csv`

## Executive summary

The MPCC model is not a direct ODE copy of the calibrated model. The MPCC uses a
reduced DC-dFBA structure where the `ZentenoKineticParameters` define dynamic
bounds for the metabolic LP. Therefore, several calibrated parameters must be
converted before they can be inserted into the MPCC.

The most important conversion is that the calibrated ODE uses an effective
reformulation:

```text
sN  = mu0 / Kn0
qN  = mu0 / YXN
qXG = mu0 / YXG
qXF = mu0 / YXF
sG  = betaG0 / Kg0
sF  = betaF0 / Kf0
qEG = betaG0 / YEG
qEF = betaF0 / YEF
iG  = 1 / Kig0
iE  = 1 / Kie0
```

The inverse mapping is the one that should be used to populate the MPCC.

Recommended policy:

1. Load the mapped values below as a `scenario_first` MPCC parameter set.
2. Run sequential validation before any full MPCC optimization.
3. Do not promote the set to the default MPCC config until trajectory-level
   checks confirm sugar, ethanol, nitrogen, biomass and aroma behavior.
4. Treat CO2, O2, glycerol, acetate, acetaldehyde, pyruvate, dead biomass and
   condensate states as extensions unless they are explicitly added to the MPCC
   state vector.

## Verified MPCC contract

The current MPCC state vector has 19 states:

```text
biomass, nitrogen, glucose, fructose, ethanol, oxygen,
protein, carbohydrate,
phenylalanine, leucine, valine, methionine, tyrosine,
phenylethanol, isoamyl_acetate, isobutanol, methionol,
tyrosol, ethyl_acetate
```

Important compatibility notes:

- `r_1862` should be interpreted as `isoamyl_acetate`, not isoamyl alcohol.
- `oxygen` currently has `dO/dt = 0` in the MPCC RHS unless the O2 extension is
  explicitly implemented.
- `glycerol`, `dead_biomass`, `pyruvate`, `acetaldehyde`, `acetate`, liquid CO2
  and aroma condensate are not part of the current MPCC state vector.
- Aroma volatilization exists as a hook, but the calibrated gas-liquid model is
  richer than the current MPCC implementation.

## Selected calibrated model structures

The calibrated pilot-2025 model selected/inherited the following structures:

| Block | Selected structure | Interpretation |
|---|---|---|
| Primary fermentation | reduced primary + glycerol + secondary v2 | Core sugar, nitrogen, biomass, ethanol, glycerol and secondary states |
| Secondary chemistry | `secondary_full_chem_o2fixed` | Pyruvate, acetaldehyde and acetate chemistry with O2-related terms fixed/limited |
| Aroma model | `ea_ethanol_nlimited` | Ethyl acetate depends on ethanol, biomass and nitrogen limitation |
| CO2 model | `solubility_o2_slow_transition` | CO2 release includes dissolution/saturation and early O2 modulation |
| Gas-liquid aroma loss | ethanol-water-sugar equilibrium layer + CO2 stripping factor | Loss scales with gas release and liquid-gas partitioning |

## Primary kinetic homologation table

Use this table as the main hand-off table for `ZentenoKineticParameters`.

| mpcc_parameter_name | advanced_model_parameter_name | advanced_model_value | advanced_model_unit | mpcc_value | mpcc_unit | conversion_formula | fixed_or_scenario | confidence | notes |
|---|---:|---:|---|---:|---|---|---|---|---|
| `MU0` | `mu0` | 0.068971078 | h^-1 | 0.068971078 | h^-1 | `MU0 = mu0` | scenario_first | medium | Direct map. Lower than MPCC nominal 0.1800. |
| `Kn0` | `mu0, sN` | `mu0=0.068971078; sN=8.751445871` | mixed effective | 0.007881107 | kg/m3 or g/L N | `Kn0 = mu0 / sN` | scenario_first | medium_low | Slightly below the old MPCC nominal 0.010. Keep as calibrated-effective value, but validate bounds. |
| `YXN` | `mu0, qN` | `mu0=0.068971078; qN=0.015771507` | mixed effective | 4.37314436 | MPCC yield convention | `YXN = mu0 / qN` | scenario_first | medium_low | Much lower than old MPCC nominal 19.69. This will increase N demand; validate N trajectories. |
| `YXG` | `mu0, qXG` | `mu0=0.068971078; qXG=0.081490986` | mixed effective | 0.846364506 | MPCC yield convention | `YXG = mu0 / qXG` | scenario_first | medium | Direct inverse reformulation. |
| `YXF` | `mu0, qXF` | `mu0=0.068971078; qXF=0.070933203` | mixed effective | 0.972345225 | MPCC yield convention | `YXF = mu0 / qXF` | scenario_first | medium | Direct inverse reformulation. |
| `betaG0` | `betaG0` | 1.186335615 | h^-1 equivalent | 1.186335615 | MPCC beta convention | `betaG0 = betaG0` | scenario_first | medium | Much higher than old nominal 0.225. Validate sugar drying speed. |
| `betaF0` | `betaF0` | 0.314509525 | h^-1 equivalent | 0.314509525 | MPCC beta convention | `betaF0 = betaF0` | scenario_first | medium | Fructose beta remains closer to old nominal than glucose beta. |
| `Kg0` | `betaG0, sG` | `betaG0=1.186335615; sG=0.097455098` | mixed effective | 12.1731508 | g/L glucose | `Kg0 = betaG0 / sG` | scenario_first | medium_low | Above old MPCC nominal 7.5 and above a common 10 g/L upper bound. Do not hard-code if old bounds still enforce 10. |
| `Kf0` | `betaF0, sF` | `betaF0=0.314509525; sF=0.1782952` | mixed effective | 1.76398195 | g/L fructose | `Kf0 = betaF0 / sF` | scenario_first | medium | Lower than old nominal 7.5; implies stronger fructose affinity in the calibrated model. |
| `YEG` | `betaG0, qEG` | `betaG0=1.186335615; qEG=1.112513547` | mixed effective | 1.06635611 | MPCC yield convention | `YEG = betaG0 / qEG` | scenario_first | low_medium | High relative to old nominal 0.49. Validate ethanol/sugar stoichiometry before promotion. |
| `YEF` | `betaF0, qEF` | `betaF0=0.314509525; qEF=1.278009637` | mixed effective | 0.246093234 | MPCC yield convention | `YEF = betaF0 / qEF` | scenario_first | low_medium | Lower than old nominal 0.49. Validate ethanol contribution from fructose. |
| `Kig0` | `iG` | 0.013970773 | 1/(g/L) | 71.5780007 | g/L glucose | `Kig0 = 1 / iG` | scenario_first | medium | Higher than old nominal 55.0. |
| `Kie0` | `iE` | 0.015660122 | 1/(g/L) | 63.8564620 | g/L ethanol | `Kie0 = 1 / iE` | scenario_first | medium | Higher than old nominal 40.0. |
| `Kd0` | `Kd0` | 0.000656095 | h^-1 equivalent | 0.000656095 | MPCC death convention | `Kd0 = Kd0` | scenario_first | low_medium | Biomass/death data are limited. Prefer scenario/prior rather than free optimization. |
| `MRATE_0` | `m0` | 0.018338579 | h^-1 equivalent | 0.018338579 | MPCC maintenance convention | `MRATE_0 = m0` | scenario_first | medium | Higher than old nominal 0.010. Validate residual sugar trajectories. |

## Julia configuration block for first sequential test

Use this only as a first validation scenario. Do not overwrite defaults before
checking trajectories.

```julia
ZentenoKineticParameters(;
    MU0 = 0.068971078,
    YXN = 4.37314436,
    YXG = 0.846364506,
    YXF = 0.972345225,
    YEG = 1.06635611,
    YEF = 0.246093234,
    Kn0 = 0.007881107,
    Kg0 = 12.1731508,
    Kf0 = 1.76398195,
    Kig0 = 71.5780007,
    Kie0 = 63.8564620,
    Kd0 = 0.000656095,
    betaG0 = 1.186335615,
    betaF0 = 0.314509525,
    MRATE_0 = 0.018338579,
)
```

If the MPCC code enforces old parameter bounds, at least these values may need
temporary scenario bounds:

```text
Kg0 = 12.1731508
YEG = 1.06635611
YXN = 4.37314436
```

Do not clip them silently. Clipping changes the calibrated model.

## Reformulation-only parameters

The following calibrated parameters should not be inserted directly into
`ZentenoKineticParameters`; they are effective reformulation variables:

| calibrated parameter | calibrated value | Use in MPCC |
|---|---:|---|
| `sN` | 8.751445871 | Use only through `Kn0 = mu0 / sN`. |
| `qN` | 0.015771507 | Use only through `YXN = mu0 / qN`. |
| `qXG` | 0.081490986 | Use only through `YXG = mu0 / qXG`. |
| `qXF` | 0.070933203 | Use only through `YXF = mu0 / qXF`. |
| `sG` | 0.097455098 | Use only through `Kg0 = betaG0 / sG`. |
| `sF` | 0.178295200 | Use only through `Kf0 = betaF0 / sF`. |
| `qEG` | 1.112513547 | Use only through `YEG = betaG0 / qEG`. |
| `qEF` | 1.278009637 | Use only through `YEF = betaF0 / qEF`. |
| `iG` | 0.013970773 | Use only through `Kig0 = 1 / iG`. |
| `iE` | 0.015660122 | Use only through `Kie0 = 1 / iE`. |

## Glycerol extension

The calibrated model includes glycerol production:

| parameter | value | Interpretation | MPCC status |
|---|---:|---|---|
| `gammaG0` | 0.070309868 | glucose-associated glycerol term | extension_required |
| `gammaF0` | 0.003394444 | fructose-associated glycerol term | extension_required |

The current MPCC state vector does not include `glycerol`. These values should
not be mapped into an existing MPCC parameter. If glycerol is added later, use
these as first scenario parameters and validate against pilot data.

## Secondary state extension

The calibrated model includes pyruvate, acetaldehyde and acetate states. The
current MPCC does not. These parameters should be delivered to the MPCC agent as
documentation, but marked `extension_required`.

| parameter | value | Intended role | MPCC status |
|---|---:|---|---|
| `kPyrS` | 0.800000 | legacy pyruvate synthesis term | extension_required |
| `kPyrD` | 0.025000 | legacy pyruvate decay/drain term | extension_required |
| `kAldPyr` | 0.007254 | acetaldehyde from pyruvate | extension_required |
| `kAldS` | 0.350000 | legacy acetaldehyde source term | extension_required |
| `kAldRed` | 0.000010 | acetaldehyde reduction | extension_required |
| `kAcAld` | 0.016679 | acetate from acetaldehyde | extension_required |
| `kAcStress` | 0.002616 | stress acetate term | extension_required |
| `kAcAssim` | 0.000001 | acetate assimilation | extension_required |
| `kPyrS_N` | 1.150082 | nitrogen-modulated pyruvate source | extension_required |
| `kPyrS_stat` | 0.012565 | stationary pyruvate source | extension_required |
| `kPyrO2` | 0.168623 | oxygen modulation for pyruvate | extension_required |
| `kPyrDrain` | 0.002397 | pyruvate drain | extension_required |
| `kAldS_N` | 3.004251 | nitrogen-modulated acetaldehyde source | extension_required |
| `kAldS_stat` | 0.674876 | stationary acetaldehyde source | extension_required |
| `kAldO2` | 1.999963 | oxygen modulation for acetaldehyde | extension_required |

## Aroma homologation

The MPCC currently tracks six aroma states:

```text
phenylethanol, isoamyl_acetate, isobutanol,
methionol, tyrosol, ethyl_acetate
```

The calibrated pilot model explicitly estimates parameters for:

```text
ethyl_acetate, isoamyl_acetate, ethyl_octanoate
```

Therefore only `isoamyl_acetate` and `ethyl_acetate` can be linked to current
MPCC aroma states without adding new states. `ethyl_octanoate` is not currently
in the MPCC state vector.

| mpcc_parameter_name | advanced_model_parameter_name | advanced_model_value | advanced_model_unit | mpcc_value | mpcc_unit | conversion_formula | fixed_or_scenario | confidence | notes |
|---|---:|---:|---|---:|---|---|---|---|---|
| `aroma_forcing_constants[isoamyl_acetate]` | `k_IAA_growth` | 0.053515 | effective production coefficient | not direct | MPCC forcing convention | no 1:1 conversion | extension_or_fit | medium | Current MPCC uses a single forcing constant; calibrated model has growth/stationary split. Fit an equivalent constant by trajectory matching if the MPCC structure is unchanged. |
| `aroma_forcing_constants[isoamyl_acetate]` | `k_IAA_stationary` | 0.040084 | effective production coefficient | not direct | MPCC forcing convention | no 1:1 conversion | extension_or_fit | medium | Prefer implementing the growth/stationary structure if possible. |
| `aroma_loss[isoamyl_acetate]` | `alpha_IAA_loss` | 0.660648 | effective stripping coefficient | 0.660648 | gas-liquid loss convention | `loss = alpha_i * K_i_LG(T,E,G,F) * q_gas * A_i_liq` | extension_required | medium | Requires gas-liquid loss hook using CO2 release and partition coefficient. |
| `aroma_forcing_constants[ethyl_acetate]` | `k_EA_growth` | 0.000100 | effective production coefficient | not direct | MPCC forcing convention | no 1:1 conversion | extension_or_fit | low_medium | Growth term was near lower bound. Do not use as sole EA forcing. |
| `aroma_forcing_constants[ethyl_acetate]` | `k_EA_stationary` | 0.005947 | effective production coefficient | not direct | MPCC forcing convention | no 1:1 conversion | extension_or_fit | medium | Stationary EA term is more relevant than growth term. |
| `aroma_forcing_constants[ethyl_acetate]` | `k_EA_XE_Nlim` | 0.041786 | effective EA N-limited coefficient | not direct | MPCC forcing convention | no 1:1 conversion | extension_or_fit | medium | Selected EA model is `ea_ethanol_nlimited`; this term is structurally important. |
| `aroma_loss[ethyl_acetate]` | `alpha_EA_loss` | 0.200200 | effective stripping coefficient | 0.200200 | gas-liquid loss convention | `loss = alpha_i * K_i_LG(T,E,G,F) * q_gas * A_i_liq` | extension_required | medium | Requires gas-liquid loss hook. |
| `ethyl_octanoate_state` | `k_EO_growth` | 0.000898 | effective production coefficient | not available | not available | not available | extension_required | medium_low | Add `ethyl_octanoate` state before using. |
| `ethyl_octanoate_state` | `k_EO_stationary` | 0.000150 | effective production coefficient | not available | not available | not available | extension_required | medium_low | Add `ethyl_octanoate` state before using. |
| `ethyl_octanoate_loss` | `alpha_EO_loss` | 1.478193 | effective stripping coefficient | not available | not available | not available | extension_required | medium_low | Add `ethyl_octanoate` state before using. |

Recommended MPCC approach for aromas:

1. Keep current `aroma_forcing_constants` unchanged for first kinetic
   validation if no aroma structural extension is implemented.
2. If using only the existing MPCC single-constant aroma forcing, estimate
   equivalent constants by trajectory matching against the calibrated model,
   not by direct numeric copy of `k_IAA_growth` or `k_EA_XE_Nlim`.
3. If modifying the MPCC, implement state-dependent terms for IAA and EA, then
   use the calibrated parameters above.
4. Add volatilization as a smooth sink term only after base production profiles
   are stable.

## CO2 and O2 homologation

The selected CO2 model is `solubility_o2_slow_transition`.

| parameter | value | unit | Status for MPCC |
|---|---:|---|---|
| `kCO2_release_h` | 0.266927696 | h^-1 | extension_required |
| `CO2sat_scale` | 0.406634334 | dimensionless | extension_required |
| `O2_qmax_mg_gdw_h` | 0.150000 | mg O2/gDW/h | extension_required |
| `O2_K_mg_l` | 0.250000 | mg/L | extension_required |
| `O2_ana_K_mg_l` | 0.750000 | mg/L | extension_required |
| `O2_ana_hill` | 2.000000 | dimensionless | extension_required |
| `O2_crabtree_floor` | 0.080000 | dimensionless | extension_required |
| `O2sat_scale` | 1.000000 | dimensionless | extension_required |
| `O2_kLa_h` | 0.000000 | h^-1 | extension_required |
| `O2_25171_fraction` | 0.050000 | dimensionless | data curation/scenario |

The focused CO2/O2 FIM was acceptable for the two estimated CO2 parameters:

| parameter | theta | approx log std | approx 95 percent multiplier | classification |
|---|---:|---:|---:|---|
| `kCO2_release_h` | 0.266927696 | 0.293301 | 1.776900 | well_estimated |
| `CO2sat_scale` | 0.406634334 | 0.448227 | 2.407347 | moderate |

Focused CO2/O2 FIM metrics:

```text
logdet = 4.213919699
min_eigenvalue = 4.553992132
max_eigenvalue = 14.848746622
condition_number = 3.260599973
trace_inv = 0.286933306
rank_1e-8 = 2
```

Important limitation: the full aroma-coupled FIM was not recomputed as a final
expensive run here because the aroma residual evaluation is slow. The current
document therefore supports CO2/O2 homologation more strongly than full aroma
structural promotion.

## Fit quality context

These metrics explain where the calibrated model is reliable enough for MPCC
homologation and where caution is needed.

| State/block | relative RMSE | correlation | Interpretation |
|---|---:|---:|---|
| Ethanol `E` | 0.178 | 0.909 | Good trend and acceptable scale. |
| Fructose `F` | 0.569 | 0.870 | Trend captured, scale errors remain. |
| Glucose `G` | 0.539 | 0.905 | Trend captured, scale errors remain. |
| Glycerol `Gly` | 0.108 | 0.910 | Good fit, but MPCC lacks glycerol state. |
| Nitrogen `N` | 0.944 | 0.674 | Weak. Do not over-trust N-related parameters without validation. |
| Biomass `X` | 0.338 | 0.645 | Moderate. Death/viability structure remains uncertain. |
| Acetaldehyde `AcAld` | 0.113 | 0.635 | Scale is reasonable, trend moderate. MPCC lacks this state. |
| Pyruvate `Pyr` | 0.331 | 0.431 | Moderate to weak. MPCC lacks this state. |
| Ethyl acetate total | 0.455 | 0.627 | Partial fit; should not be promoted blindly. |
| Isoamyl acetate total | 0.390 | 0.775 | Best aroma among the tracked aroma states. |
| CO2 25170 | 0.779 | 0.843 | Trend acceptable, scale still noisy. |
| CO2 25171 | 0.435 | 0.931 | Good trend after using the corrected reinoculation start. |

## Parameters to leave out of the current MPCC homologation

Leave these outside the current MPCC parameter replacement unless the model is
extended:

```text
gammaG0, gammaF0
kPyrS, kPyrD, kAldPyr, kAldS, kAldRed, kAcAld, kAcStress, kAcAssim
kPyrS_N, kPyrS_stat, kPyrO2, kPyrDrain
kAldS_N, kAldS_stat, kAldO2
k_EA_growth, k_EA_stationary, k_EA_X, k_EA_XE, k_EA_XE_Nlim, k_EA_AcAld, q10_EA
k_IAA_growth, k_IAA_stationary
k_EO_growth, k_EO_stationary
alpha_EA_loss, alpha_IAA_loss, alpha_EO_loss
kCO2_release_h, CO2sat_scale
O2_qmax_mg_gdw_h, O2_K_mg_l, O2_ana_K_mg_l, O2_ana_hill,
O2_crabtree_floor, O2sat_scale, O2_kLa_h, O2_25171_fraction
```

Reason: these parameters correspond to states, losses or structural terms that
are absent or only partially represented in the current MPCC.

## Recommended validation protocol for the MPCC agent

### Step 1: sequential kinetic check

Run the MPCC sequential simulator with the proposed `ZentenoKineticParameters`
and no dynamic optimization.

Required cases:

```text
case_01: nominal MPCC initial state, T = 18 C, no nutrient addition
case_02: nominal MPCC initial state, T = 20 C, no nutrient addition
case_03: pilot-like natural must initial state, measured T profile, measured N policy
case_04: pilot-like natural must initial state, warm-to-cool T profile, no N pulse
case_05: pilot-like natural must initial state, cold-to-warm T profile, early N pulse
```

Minimum checks:

```text
S(t) = G(t) + F(t)
t_S50, t_S30, t_S5 if reached
G(tf), F(tf), E(tf), N(tf), X(tf)
min_state_value
ethanol/sugar plausibility
glucose/fructose preference
solver retcode/status
```

Acceptance rules:

- No significant negative states.
- No unrealistic ethanol production relative to sugar consumed.
- Nitrogen should not be depleted instantly unless the measured case supports it.
- Sugar drying should be faster than the old nominal MPCC if the calibrated
  parameters are active.
- If `Kg0` or `YEG` are clipped by old bounds, report the clipped trajectory as
  a separate scenario, not as the calibrated scenario.

### Step 2: fixed-control MPCC check

Run a fixed-control MPCC trajectory using an existing feasible seed. Do not
release all controls yet.

Minimum checks:

```text
max_rhs_residual
max_collocation_integration_residual
max_continuity_residual
max_lower_complementarity_residual
selected_candidate_residual_feasible
selected_candidate_trajectory_ready
terminal_total_sugar_g_L
terminal_isoamyl_acetate_g_L
terminal_ethyl_acetate_g_L
```

### Step 3: aroma check

If the MPCC aroma structure is unchanged, compare only qualitative behavior.
If the aroma production/loss extension is implemented, compare:

```text
isoamyl_acetate retained
ethyl_acetate retained
condensate-equivalent ethyl_acetate
condensate-equivalent isoamyl_acetate
loss_fraction = condensate / (retained + condensate)
```

## Selected DOE campaign context

The current focused CO2/O2 MBDoE campaign selected these natural-must
experiments:

| Order | Design | Temperature segments | N pulse |
|---:|---|---|---|
| 1 | `natural_pilot_noN_dynamic_temperature` | 15, 22, 18, 22 C | none |
| 2 | `natural_pilot_CO2_warm_start_cool_retention` | 23, 22, 16, 15 C | none |
| 3 | `natural_pilot_warm_to_cool_noN` | 22, 22, 17, 16 C | none |
| 4 | `natural_pilot_EA_warm_early_noN` | 24, 23, 19, 17 C | none |
| 5 | `natural_pilot_CO2_Npulse_release_probe` | 18, 18, 23, 20 C | 50 h: 0.045 kg/m3 |
| 6 | `natural_pilot_cold_to_warm_earlyN` | 14, 16, 22, 20 C | 30 h: 0.045 kg/m3 |

This campaign is useful context for future data generation, but it is not a
replacement for validating the homologated MPCC trajectories.

## Final recommendation

For immediate MPCC work, only homologate the primary kinetic block:

```text
MU0, Kn0, YXN, YXG, YXF, betaG0, betaF0, Kg0, Kf0,
YEG, YEF, Kig0, Kie0, Kd0, MRATE_0
```

Keep these as a scenario parameter set. Run sequential validation first. Then
run fixed-control MPCC. Only after those pass should the agent modify aroma,
CO2/O2 or secondary-state structure.

The highest-risk direct replacements are:

```text
YXN = 4.37314436
Kg0 = 12.1731508
YEG = 1.06635611
YEF = 0.246093234
```

These are not necessarily wrong; they are calibrated-effective values from the
ODE reformulation. They must be validated at trajectory level before they are
treated as physically final MPCC defaults.
