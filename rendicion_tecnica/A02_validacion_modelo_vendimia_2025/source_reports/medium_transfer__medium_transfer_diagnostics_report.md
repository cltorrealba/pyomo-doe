# Medium transfer diagnostics

This report separates the historical fermentation database by must medium when metadata are available.
The goal is to decide whether synthetic-must experiments can be used directly for wine-must parameter estimation, or whether natural-must data require medium-specific correction terms.

## Data audit

|   batch |   n_rows |   t_max_h |   G_initial |   F_initial |   N_initial |   E_initial |   temperatura_mean |   n_positive_nutrient_pulses |   isoamyl_acetate_n_obs_raw |   ethyl_octanoate_n_obs_raw |   ethyl_acetate_n_obs_raw |
|--------:|---------:|----------:|------------:|------------:|------------:|------------:|-------------------:|-----------------------------:|----------------------------:|----------------------------:|--------------------------:|
|   25026 |       24 |       138 |       74.99 |       72.96 |      154.5  |     6.1542  |            20.7917 |                            2 |                           5 |                           5 |                         5 |
|   25027 |       32 |       186 |       77.81 |       77.54 |      162.5  |     6.1542  |            18.8125 |                            2 |                           5 |                           5 |                         5 |
|   25085 |       35 |       216 |       46.05 |       54.34 |       93.34 |     9.79179 |            18.3    |                            2 |                           5 |                           5 |                         5 |
|   25086 |       35 |       216 |       63.28 |       72.98 |       92.7  |    10.8884  |            18.4829 |                            2 |                           5 |                           5 |                         5 |
|   25150 |       30 |       174 |       65.07 |       67.57 |       94.24 |     8.27343 |            20.3867 |                            2 |                           6 |                           6 |                         6 |
|   25151 |       30 |       174 |       63.2  |       73.6  |       90.78 |     8.27314 |            20.5467 |                            2 |                           6 |                           6 |                         6 |
|   25170 |       63 |       372 |       65.64 |       74.52 |       98.42 |     6.56092 |            15.6317 |                            2 |                          10 |                          10 |                        10 |
|   25171 |       63 |       372 |       63.97 |       77.1  |       96.6  |     7.04976 |            14.8524 |                            2 |                           8 |                           8 |                         8 |

## Medium metadata status

No usable synthetic/natural medium mapping was provided or found in the Excel file.
Fill `medium_map_template.csv` with `medium_type = synthetic` or `natural`, then rerun this script with `--medium-map C:\Users\ctorrealba\OneDrive - Viña Concha y Toro S.A\Documentos\Doctorado\Artículos\Artículo - Estimación_dFBA\pyomo-doe\fermentation_model\results\medium_transfer_diagnostics\medium_map_template.csv`.

Until that mapping exists, the current DOE can use the historical data as a pooled prior, but it cannot distinguish model-structure limitations from medium-transfer limitations.
