# Extended fermentation campaign DOE

This report screens synthetic-must experiments using the Pyomo DoE sequential FIM workflow.

Extended design actions: temperature setpoint profile, initial synthetic-must composition, and fixed candidate pulses of YAN, glucose, fructose, ethanol, and viable biomass.

Extended model states: viable biomass `X`, dead/non-viable biomass `Xd`, assimilable nitrogen `N`, glucose `G`, fructose `F`, ethanol `E`, and cumulative `CO2`.

CO2 is represented here as cumulative mass-equivalent production with fixed stoichiometry `dCO2/dt = (44.01/46.07) * dE_prod/dt`. Zenteno (2010) includes a richer CO2 readout with growth-associated production and gas-liquid transfer/saturation; use that layer before final calibration if the online sensor reports gas flow instead of cumulative CO2.

Manual liquid/biomass sampling is assumed every 12 h. CO2 is assumed online and downsampled every 6 h in the FIM to avoid overwhelming the manual observations.

## Selected Campaign

|   campaign_order | candidate                           | family           |   horizon_h | temperature_c   |   campaign_logdet |   campaign_min_relative_eigenvalue |   campaign_condition_number |   correctable_mean_var_reduction |   correctable_worst_var_reduction | coverage_note           | rationale                                                                                                       |
|-----------------:|:------------------------------------|:-----------------|------------:|:----------------|------------------:|-----------------------------------:|----------------------------:|---------------------------------:|----------------------------------:|:------------------------|:----------------------------------------------------------------------------------------------------------------|
|                1 | combined_stress_long_horizon        | integrated       |         240 | 15, 22, 25, 18  |           86.5682 |                        3.60933e-06 |                    277060   |                         0.721955 |                          0.157033 |                         | Broad stress/input design for integrated all-parameter information.                                             |
|                2 | fructose_rich_iG_probe              | sugar_initial    |         192 | 18, 20, 24, 22  |           97.8867 |                        3.32445e-06 |                    300802   |                         0.83445  |                          0.368519 |                         | Fructose-rich must plus glucose pulse separates fructose capacity from glucose inhibition iG.                   |
|                3 | glucose_rich_growth_yield           | sugar_initial    |         168 | 18, 22, 24, 20  |          103.761  |                        6.52118e-06 |                    153346   |                         0.888583 |                          0.452026 |                         | High glucose must isolates glucose growth/fermentation terms and tests fructose response after glucose history. |
|                4 | high_biomass_low_N_maintenance      | maintenance      |         168 | 16, 18, 20, 18  |          106.728  |                        8.99366e-06 |                    111189   |                         0.896305 |                          0.45796  |                         | High biomass and low N create low-growth sugar consumption windows for m0.                                      |
|                5 | low_temp_growth_separation_plus_co2 | legacy_plus      |         168 | 15, 17, 20, 22  |          108.979  |                        1.08426e-05 |                     92228.9 |                         0.910785 |                          0.52468  |                         | Previous robust N/T design, now measured with Xd and online CO2.                                                |
|                6 | ethanol_pulse_death_probe           | death            |         216 | 20, 25, 25, 22  |          111.006  |                        1.00861e-05 |                     99146   |                         0.920432 |                          0.546496 |                         | Late ethanol stress with Xd observation targets iE/Kd0 separately from growth.                                  |
|                7 | ethanol_initial_challenge           | ethanol          |         192 | 18, 22, 25, 22  |          112.526  |                        1.00801e-05 |                     99205   |                         0.927214 |                          0.578202 |                         | Initial ethanol decouples ethanol inhibition iE from ethanol produced by fermentation.                          |
|                8 | glucose_pulse_after_N_depletion     | pulse_separation |         192 | 17, 20, 24, 22  |          113.771  |                        1.20621e-05 |                     82904.5 |                         0.933967 |                          0.60994  |                         | Glucose pulse after likely N limitation separates fermentation/maintenance from growth uptake.                  |
|                9 | viable_biomass_step                 | biomass_input    |         168 | 18, 22, 22, 18  |          114.721  |                        1.1421e-05  |                     87557.9 |                         0.936081 |                          0.618248 | forced X-pulse coverage | Known viable biomass addition tests whether rates scale with X and improves q/yield separation.                 |

## Top Single-Experiment Candidates

| candidate                           | family           |   horizon_h |   combined_logdet |   combined_min_relative_eigenvalue |   combined_condition_number |   correctable_mean_var_reduction |   correctable_worst_var_reduction | status   |
|:------------------------------------|:-----------------|------------:|------------------:|-----------------------------------:|----------------------------:|---------------------------------:|----------------------------------:|:---------|
| glucose_rich_growth_yield           | sugar_initial    |         168 |           85.7898 |                        2.91647e-06 |                      342880 |                         0.761066 |                         0.332089  | ok       |
| fructose_rich_iG_probe              | sugar_initial    |         192 |           86.1654 |                        3.62685e-06 |                      275721 |                         0.739823 |                         0.293045  | ok       |
| combined_stress_long_horizon        | integrated       |         240 |           86.5682 |                        3.60933e-06 |                      277060 |                         0.721955 |                         0.157033  | ok       |
| ethanol_initial_challenge           | ethanol          |         192 |           81.7083 |                        6.53786e-06 |                      152955 |                         0.698356 |                         0.122962  | ok       |
| low_temp_growth_separation_plus_co2 | legacy_plus      |         168 |           75.4044 |                        3.02218e-06 |                      330887 |                         0.570269 |                         0.40287   | ok       |
| ethanol_pulse_death_probe           | death            |         216 |           77.0489 |                        1.72733e-06 |                      578930 |                         0.522228 |                         0.0725297 | ok       |
| viable_biomass_step                 | biomass_input    |         168 |           74.8229 |                        1.62444e-06 |                      615598 |                         0.506829 |                         0.289874  | ok       |
| low_sugar_saturation_scan           | saturation       |         144 |           70.7452 |                        2.86251e-06 |                      349343 |                         0.456491 |                         0.111601  | ok       |
| glucose_pulse_after_N_depletion     | pulse_separation |         192 |           72.1559 |                        1.81534e-06 |                      550862 |                         0.444821 |                         0.0577146 | ok       |
| high_biomass_low_N_maintenance      | maintenance      |         168 |           74.6381 |                        2.35886e-06 |                      423934 |                         0.430751 |                         0.0144656 | ok       |
| fructose_pulse_after_N_depletion    | pulse_separation |         192 |           72.1929 |                        1.73048e-06 |                      577874 |                         0.352815 |                         0.0397875 | ok       |
| yan_saturation_scan                 | saturation       |         144 |           69.6918 |                        2.10711e-06 |                      474583 |                         0.305407 |                         0.0318964 | ok       |

## Final Parameter Variance Reductions

| parameter   | group                  |   campaign_var_ratio |   campaign_var_reduction |
|:------------|:-----------------------|---------------------:|-------------------------:|
| mu0         | accepted_current       |          0.0769238   |                 0.923076 |
| sN          | currently_fixed_target |          0.381752    |                 0.618248 |
| qN          | accepted_current       |          0.107036    |                 0.892964 |
| qXG         | currently_fixed_target |          0.0186115   |                 0.981388 |
| qXF         | currently_fixed_target |          0.018446    |                 0.981554 |
| betaG0      | accepted_current       |          0.0132449   |                 0.986755 |
| sG          | currently_fixed_target |          0.0180666   |                 0.981933 |
| betaF0      | accepted_current       |          0.00345121  |                 0.996549 |
| sF          | currently_fixed_target |          0.0370674   |                 0.962933 |
| qEG         | accepted_current       |          0.0195648   |                 0.980435 |
| qEF         | accepted_current       |          0.00356892  |                 0.996431 |
| iG          | accepted_current       |          0.00271512  |                 0.997285 |
| iE          | currently_fixed_target |          0.0160981   |                 0.983902 |
| Kd0         | currently_fixed_target |          0.000945029 |                 0.999055 |
| m0          | currently_fixed_target |          0.0203621   |                 0.979638 |

## Interpretation

- `qXG/qXF` are attacked through orthogonal glucose/fructose initial compositions and single-sugar pulses.
- `sG/sF/sN` are attacked through low/intermediate substrate and YAN ladders rather than only rich musts.
- `iE` is attacked by initial/pulsed ethanol so inhibition is not only generated endogenously by fermentation.
- `Kd0` is attacked by adding `Xd` as an observed state; without `Xd`, death remains confounded with slow growth or low viability.
- `m0` is attacked through high-biomass, low-N sugar consumption windows where growth is constrained.