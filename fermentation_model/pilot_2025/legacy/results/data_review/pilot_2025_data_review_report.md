# Pilot 2025 data review

## Scope

This first pass loads the pilot-scale natural-must fermentations, normalizes units where assumptions are clear, aligns online CO2 files to fermentation start times, and checks which states are informative for calibration and MBDoE.

## Unit assumptions

- `ETANOL` is treated as g/L. Values reach about 95-98, which is plausible as g/L and impossible as % v/v for wine.
- Viable cell counts are converted with 30 pg/cell, so 1 million cells/mL = 0.03 kg/m3.
- `PAN`, `AMMONIA`, `YAN`, and pyruvic/acetaldehyde observations are kept as mg/L.
- Aroma `*_total` columns are interpreted as wine-retained plus condenser-accumulated equivalent concentration in mg/L.
- Aroma `*_condensate` columns are interpreted as equivalent must concentration accumulated in the condenser in mg/L.

## Data volume

- Fermentation sheets loaded: 8.
- Total calibration rows: 312.
- Raw CO2 sensor batches loaded: 4.
- Curated CO2 sensor batches retained for calibration: 2.

## Density to sugar proxy

Linear fit: `S_GF_g_l = -2024.539 + 2.036 density`, `R2 = 0.987`, `RMSE = 5.82 g/L`, `n = 136`.

## Observation coverage

|   batch |   n_rows |   t_min_h |   t_max_h |   n_temperature_c |   n_density |   n_S_GF_g_l |   n_X_viable_kg_m3 |   n_X_dry_weight_g_l |   n_YAN_mg_l |   n_PAN_mg_l |   n_NH4_mg_l |   n_glycerol_g_l |   n_pyruvic_acid_mg_l |   n_acetaldehyde_mg_l |   n_acetic_acid_g_l |   n_E_g_l |   n_ethyl_acetate_total |   n_ethyl_acetate_condensate |   n_isoamyl_acetate_total |   n_isoamyl_acetate_condensate |   n_ethyl_octanoate_total |   n_ethyl_octanoate_condensate |
|--------:|---------:|----------:|----------:|------------------:|------------:|-------------:|-------------------:|---------------------:|-------------:|-------------:|-------------:|-----------------:|----------------------:|----------------------:|--------------------:|----------:|------------------------:|-----------------------------:|--------------------------:|-------------------------------:|--------------------------:|-------------------------------:|
|   25026 |       24 |         0 |       138 |                24 |          24 |           11 |                 11 |                   10 |           11 |           11 |           11 |               11 |                    11 |                     0 |                   0 |        11 |                       5 |                            4 |                         5 |                              4 |                         5 |                              4 |
|   25027 |       32 |         0 |       186 |                32 |          32 |           13 |                 13 |                   12 |           13 |           13 |           13 |               13 |                    13 |                     0 |                   0 |        13 |                       5 |                            4 |                         5 |                              4 |                         5 |                              4 |
|   25085 |       35 |         0 |       216 |                35 |          34 |           17 |                 17 |                   17 |           17 |           17 |           17 |               17 |                    17 |                    17 |                   0 |        17 |                       5 |                            4 |                         5 |                              4 |                         5 |                              4 |
|   25086 |       35 |         0 |       216 |                35 |          34 |           17 |                 17 |                   17 |           17 |           17 |           17 |               17 |                    17 |                    17 |                   0 |        17 |                       5 |                            4 |                         5 |                              4 |                         5 |                              4 |
|   25150 |       30 |         0 |       174 |                30 |          30 |           14 |                 14 |                   14 |           14 |           14 |           14 |               14 |                    13 |                    14 |                   0 |        14 |                       6 |                            5 |                         6 |                              5 |                         6 |                              5 |
|   25151 |       30 |         0 |       174 |                30 |          30 |           14 |                 14 |                   14 |           14 |           14 |           14 |               14 |                    13 |                    14 |                   0 |        14 |                       6 |                            5 |                         6 |                              5 |                         6 |                              5 |
|   25170 |       63 |         0 |       372 |                63 |          63 |           26 |                 23 |                   26 |           24 |           24 |           24 |               26 |                    26 |                    26 |                   0 |        26 |                      11 |                           10 |                        11 |                             10 |                        11 |                             10 |
|   25171 |       63 |         0 |       372 |                63 |          63 |           26 |                 22 |                   26 |           24 |           24 |           24 |               26 |                    26 |                    26 |                   0 |        26 |                       9 |                            8 |                         9 |                              8 |                         9 |                              8 |

## CO2 curation decisions

|   batch |   raw_rows | used_for_calibration   | reason                                                 |   activation_time_h |   curated_rows |
|--------:|-----------:|:-----------------------|:-------------------------------------------------------|--------------------:|---------------:|
|   25150 |         56 | False                  | excluded: CO2 file not usable for this batch           |                 nan |              0 |
|   25151 |         56 | False                  | excluded: CO2 file not usable for this batch           |                 nan |              0 |
|   25170 |      22877 | True                   | kept                                                   |                   0 |          22837 |
|   25171 |      20000 | True                   | trimmed to operational CO2/process activation override |                 100 |          13963 |

## CO2 sensor coverage

|   batch |   n_rows |   t_original_min_h |   t_original_max_h |   t_effective_min_h |   t_effective_max_h |   co2_min |   co2_median |   co2_max |
|--------:|---------:|-------------------:|-------------------:|--------------------:|--------------------:|----------:|-------------:|----------:|
|   25170 |    22837 |         0.00583333 |            380.711 |          0.00583333 |             380.711 |         0 |     0.305313 |   4.79852 |
|   25171 |    13963 |       100.007      |            332.777 |          0.00694444 |             232.777 |         0 |     0.406425 |   3.61356 |

## Initial interpretation

- The pilot data are mostly natural-must process trajectories, so design variables are less flexible than in synthetic-lab DOE.
- The strongest immediate calibration value is the richer time series for ethanol, glycerol, pyruvate, acetaldehyde, ethyl acetate, isoamyl acetate, and condenser aroma accumulation.
- Online CO2 is available for a subset of batches. Batches 25150 and 25151 are excluded from calibration. Batch 25171 is trimmed to the sustained CO2 activation point and carries both original and effective time columns.
- The next step is to map these normalized batches into the extended model residual structure and run calibration/estimability with the existing reduced secondary/aroma framework.
