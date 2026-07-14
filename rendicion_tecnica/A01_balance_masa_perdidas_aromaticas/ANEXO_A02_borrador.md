# Anexo A02 - Desarrollo y calibracion de modelos de balance de masa para cuantificar perdidas aromaticas

## Identificacion

- Objetivo especifico: OE2
- Actividad: 2.13
- Nombre actividad: Desarrollar y calibrar modelos de balance de masa para cuantificar perdidas aromaticas
- Periodo Gantt: 2025-05-01 a 2025-09-30
- Fecha de cierre documental del bundle: 2026-06-22
- Estado propuesto de evidencia: cobertura documental completa de la evidencia minima requerida.

## Resumen ejecutivo

Se desarrollo y calibro un modelo dinamico de balance de masa para fermentacion que integra biomasa, nitrogeno asimilable, glucosa, fructosa, etanol y glicerol. El modelo fue parametrizado y evaluado usando bases VL3, mosto sintetico y datos natural/sintetico homologados. Sobre esta base se implemento una capa mecanistica de perdidas aromaticas por stripping con CO2, con particion gas-liquido derivada de calculos UNIFAC y una prediccion de condensado terminal para cierre de masa.

El paquete A02 cubre las categorias solicitadas: ecuaciones y supuestos, datos de calibracion, parametros, cierre de masa, metricas de ajuste, version de codigo y analisis de sensibilidad/estimabilidad.

## Matriz de cumplimiento

| evidencia_minima_requerida         | estado_bundle   | archivos_bundle                                                                                    | detalle                                                                                                 |
|:-----------------------------------|:----------------|:---------------------------------------------------------------------------------------------------|:--------------------------------------------------------------------------------------------------------|
| Ecuaciones y supuestos del balance | cubierto        | 03_ecuaciones_y_supuestos_A02.md; code_sources/*.py                                                | Incluye balances base X/Xd/N/G/F/E/Gly, capa A_liq/A_loss/A_cond y supuestos UNIFAC/CO2.                |
| Datos de calibracion               | cubierto        | data_sources/*.xlsx; tables/02_resumen_datos_fuente.csv; tables/02_diccionario_variables_clave.csv | Incluye libros VL3 y mosto sintetico, mas resumen de hojas, columnas y variables.                       |
| Parametros                         | cubierto        | tables/05_parametros_cineticos.csv; tables/05_parametros_particion_aromas.csv                      | Incluye theta calibrado y parametros de particion aromaticos derivados de UNIFAC.                       |
| Cierre de masa                     | cubierto        | tables/06_cierre_masa_aromatico.csv; 03_ecuaciones_y_supuestos_A02.md                              | Cierre por especie/candidato: integral de sintesis, liquido final, perdida final y condensado esperado. |
| Metricas de ajuste                 | cubierto        | tables/04_metricas_calibracion.csv; source_reports/*fit*; source_reports/*glycerol*                | Incluye WSSE, residuales, grados de libertad, exito de ajuste y capacidad de ajuste secundaria.         |
| Version de codigo                  | cubierto        | 08_version_codigo.md; 00_manifest.csv                                                              | Incluye commit, rama, estado git, diff stat y hashes SHA256 de artefactos.                              |
| Analisis de sensibilidad           | cubierto        | tables/07_sensibilidad_estimabilidad.csv; source_reports/*eigen*; source_reports/*campaign*        | Incluye FIM, reduccion de varianza, ranking DOE, eigenvalores y direcciones debiles.                    |

## Datos de calibracion

Los datos fuente principales son:

- `Calibration_data_vl3.xlsx`: fermentaciones historicas VL3 con azucares, nitrogeno, glicerol, piruvato, acetaldehido, etanol y aromas totales.
- `mosto_sintetico_vl3.xlsx`: datos homologados de mosto sintetico, diseno CCD, datos originales, flags de calidad y diccionario de mapeo.
- `mosto_natural_xthiol.xlsx`: fuente complementaria natural usada por el flujo de carga cuando esta disponible.

Resumen de hojas:

| source_file                                       | sheet             |   data_rows_estimated |   columns | aroma_total_columns                                                                                                                |
|:--------------------------------------------------|:------------------|----------------------:|----------:|:-----------------------------------------------------------------------------------------------------------------------------------|
| fermentation_model/data/Calibration_data_vl3.xlsx | 25026             |                    64 |        25 | benzaldehido_total;hexil_acetate_total;isoamil_acetate_total;octanoate_de_etilo_total;phenylethylacetate_total;Ethyl_Acetate_total |
| fermentation_model/data/Calibration_data_vl3.xlsx | 25027             |                    64 |        26 | benzaldehido_total;hexil_acetate_total;isoamil_acetate_total;octanoate_de_etilo_total;phenylethylacetate_total;Ethyl_Acetate_total |
| fermentation_model/data/Calibration_data_vl3.xlsx | 25085             |                    64 |        25 | benzaldehido_total;hexil_acetate_total;isoamil_acetate_total;octanoate_de_etilo_total;phenylethylacetate_total;Ethyl_Acetate_total |
| fermentation_model/data/Calibration_data_vl3.xlsx | 25086             |                    64 |        25 | benzaldehido_total;hexil_acetate_total;isoamil_acetate_total;octanoate_de_etilo_total;phenylethylacetate_total;Ethyl_Acetate_total |
| fermentation_model/data/Calibration_data_vl3.xlsx | 25150             |                    64 |        25 | benzaldehido_total;hexil_acetate_total;isoamil_acetate_total;octanoate_de_etilo_total;phenylethylacetate_total;Ethyl_Acetate_total |
| fermentation_model/data/Calibration_data_vl3.xlsx | 25151             |                    64 |        26 | benzaldehido_total;hexil_acetate_total;isoamil_acetate_total;octanoate_de_etilo_total;phenylethylacetate_total;Ethyl_Acetate_total |
| fermentation_model/data/Calibration_data_vl3.xlsx | 25170             |                    63 |        25 | benzaldehido_total;hexil_acetate_total;isoamil_acetate_total;octanoate_de_etilo_total;phenylethylacetate_total;Ethyl_Acetate_total |
| fermentation_model/data/Calibration_data_vl3.xlsx | 25171             |                    63 |        25 | benzaldehido_total;hexil_acetate_total;isoamil_acetate_total;octanoate_de_etilo_total;phenylethylacetate_total;Ethyl_Acetate_total |
| fermentation_model/data/mosto_sintetico_vl3.xlsx  | 00_Resumen        |                    10 |         2 |                                                                                                                                    |
| fermentation_model/data/mosto_sintetico_vl3.xlsx  | Datos_homologados |                    86 |        92 | benzaldehido_total;hexil_acetate_total;isoamil_acetate_total;octanoate_de_etilo_total;phenylethylacetate_total;Ethyl_Acetate_total |
| fermentation_model/data/mosto_sintetico_vl3.xlsx  | MS007             |                    10 |        92 | benzaldehido_total;hexil_acetate_total;isoamil_acetate_total;octanoate_de_etilo_total;phenylethylacetate_total;Ethyl_Acetate_total |
| fermentation_model/data/mosto_sintetico_vl3.xlsx  | MS008             |                    10 |        92 | benzaldehido_total;hexil_acetate_total;isoamil_acetate_total;octanoate_de_etilo_total;phenylethylacetate_total;Ethyl_Acetate_total |
| fermentation_model/data/mosto_sintetico_vl3.xlsx  | MS009             |                    10 |        92 | benzaldehido_total;hexil_acetate_total;isoamil_acetate_total;octanoate_de_etilo_total;phenylethylacetate_total;Ethyl_Acetate_total |
| fermentation_model/data/mosto_sintetico_vl3.xlsx  | MS010             |                    10 |        92 | benzaldehido_total;hexil_acetate_total;isoamil_acetate_total;octanoate_de_etilo_total;phenylethylacetate_total;Ethyl_Acetate_total |
| fermentation_model/data/mosto_sintetico_vl3.xlsx  | MS011             |                    10 |        92 | benzaldehido_total;hexil_acetate_total;isoamil_acetate_total;octanoate_de_etilo_total;phenylethylacetate_total;Ethyl_Acetate_total |
| fermentation_model/data/mosto_sintetico_vl3.xlsx  | MS012             |                    10 |        92 | benzaldehido_total;hexil_acetate_total;isoamil_acetate_total;octanoate_de_etilo_total;phenylethylacetate_total;Ethyl_Acetate_total |
| fermentation_model/data/mosto_sintetico_vl3.xlsx  | MS013             |                     9 |        92 | benzaldehido_total;hexil_acetate_total;isoamil_acetate_total;octanoate_de_etilo_total;phenylethylacetate_total;Ethyl_Acetate_total |
| fermentation_model/data/mosto_sintetico_vl3.xlsx  | MS014             |                     9 |        92 | benzaldehido_total;hexil_acetate_total;isoamil_acetate_total;octanoate_de_etilo_total;phenylethylacetate_total;Ethyl_Acetate_total |

## Modelo y supuestos

El detalle formal esta en `03_ecuaciones_y_supuestos_A02.md`. En sintesis, el modelo base resuelve balances de masa para `X, Xd, N, G, F, E, Gly`; la extension aromatica agrega estados `A_liq`, `A_loss` y `A_cond` para cuantificar aroma retenido, aroma perdido y condensado terminal.

La perdida aromatica se calcula como:

```text
r_loss(i,t) = alpha_i K_lg(i,t) Q_CO2(t) A_liq(i,t)
```

y el cierre acumulado queda:

```text
dA_liq/dt  = r_syn - r_loss
dA_loss/dt = r_loss
A_cond     = eta A_loss(t_final)
```

## Parametros y calibracion

Metricas principales de calibracion:

| workflow                           | fit                 | fit_label                                     |   n_batches | mediums           |   n_parameters | success   |   initial_wsse |   final_wsse |   n_residuals |   wsse_per_residual |   wsse_per_dof |
|:-----------------------------------|:--------------------|:----------------------------------------------|------------:|:------------------|---------------:|:----------|---------------:|-------------:|--------------:|--------------------:|---------------:|
| new_must_glycerol_estimability_doe | natural_reduced11   | nan                                           |           9 | natural           |             11 | True      |      54080.6   |  13409.5     |           608 |           22.0551   |        22.4615 |
| new_must_glycerol_estimability_doe | synthetic_reduced11 | nan                                           |          10 | synthetic         |             11 | True      |      54936.9   |  10723.8     |           542 |           19.7856   |        20.1955 |
| new_must_glycerol_estimability_doe | mixed_reduced11     | nan                                           |          19 | natural,synthetic |             11 | True      |     109017     |  27138.9     |          1150 |           23.5991   |        23.827  |
| new_must_glycerol_estimability_doe | mixed_full17_l2     | nan                                           |          19 | natural,synthetic |             17 | True      |      27167.6   |  26132.3     |          1160 |           22.5278   |        22.8629 |
| secondary_fit_capacity             | nan                 | global_v2_reduced_o2fixed                     |          27 | nan               |            nan | True      |        nan     |   5754.09    |           588 |            9.78587  |       nan      |
| secondary_fit_capacity             | nan                 | medium_v2_reduced_o2fixed_historical_vl3      |           8 | nan               |            nan | True      |       2503.15  |   2418.63    |           251 |            9.63597  |       nan      |
| secondary_fit_capacity             | nan                 | medium_v2_reduced_o2fixed_natural             |           9 | nan               |            nan | True      |       2647.61  |   1949.82    |           253 |            7.70681  |       nan      |
| secondary_fit_capacity             | nan                 | medium_v2_reduced_o2fixed_synthetic           |          10 | nan               |            nan | True      |        603.339 |    365.77    |            87 |            4.20426  |       nan      |
| secondary_fit_capacity             | nan                 | batch_v2_reduced_o2fixed_historical_vl3_25026 |           1 | nan               |            nan | True      |        123.343 |      6.33316 |            11 |            0.575742 |       nan      |
| secondary_fit_capacity             | nan                 | batch_v2_all_o2free_historical_vl3_25026      |           1 | nan               |            nan | True      |        123.343 |      4.75219 |            11 |            0.432017 |       nan      |
| secondary_fit_capacity             | nan                 | batch_v2_reduced_o2fixed_historical_vl3_25027 |           1 | nan               |            nan | True      |        260.455 |     61.4148  |            13 |            4.72422  |       nan      |
| secondary_fit_capacity             | nan                 | batch_v2_all_o2free_historical_vl3_25027      |           1 | nan               |            nan | False     |        260.455 |     32.8603  |            13 |            2.52771  |       nan      |

Parametros cineticos:

| source_file                                                                       | fit                 |      mu0 |      sN |         qN |       qXG |      qXF |   betaG0 |        sG |   betaF0 |        sF |      qEG |      qEF |         iG |         iE |         Kd0 |         m0 |   gammaG0 |     gammaF0 |
|:----------------------------------------------------------------------------------|:--------------------|---------:|--------:|-----------:|----------:|---------:|---------:|----------:|---------:|----------:|---------:|---------:|-----------:|-----------:|------------:|-----------:|----------:|------------:|
| fermentation_model/results/new_must_glycerol_estimability_doe/theta_fit_table.csv | natural_reduced11   | 0.10091  | 18      | 0.0289682  | 0.1125    | 0.1125   | 0.502511 | 0.03      | 0.487841 | 0.03      | 1.17707  | 1.04431  | 0.0068611  | 0.00789219 | 0.0037451   | 0.01       | 0.0751174 | 0.00505987  |
| fermentation_model/results/new_must_glycerol_estimability_doe/theta_fit_table.csv | synthetic_reduced11 | 0.050014 | 18      | 0.00863748 | 0.1125    | 0.1125   | 0.554456 | 0.03      | 0.304715 | 0.03      | 1.12466  | 0.735091 | 0.00309291 | 0.0141443  | 0.000908616 | 0.01       | 0.0711298 | 0.000222185 |
| fermentation_model/results/new_must_glycerol_estimability_doe/theta_fit_table.csv | mixed_reduced11     | 0.073746 | 18      | 0.0173737  | 0.1125    | 0.1125   | 0.40601  | 0.03      | 0.288099 | 0.03      | 0.941025 | 0.636993 | 0.00305606 | 0.00589815 | 0.000888893 | 0.01       | 0.0607084 | 0.000158468 |
| fermentation_model/results/new_must_glycerol_estimability_doe/theta_fit_table.csv | mixed_full17_l2     | 0.061443 | 22.8223 | 0.0140776  | 0.0830335 | 0.147088 | 0.429577 | 0.0382373 | 0.286346 | 0.0509012 | 0.998681 | 0.600556 | 0.00330299 | 0.00536305 | 0.000772526 | 0.00451924 | 0.0584596 | 0.00312368  |

Parametros de particion aromatica:

| workflow                 | species         |        K20 |   temp_slope |   ethanol_slope |   sugar_slope |   fit_rmse_log |   fit_max_abs_log_error |   trap_efficiency |
|:-------------------------|:----------------|-----------:|-------------:|----------------:|--------------:|---------------:|------------------------:|------------------:|
| aroma_campaign_doe       | ethyl_acetate   | 0.00443687 |    0.0419447 |     -0.00443587 |   -0.00129469 |      0.0161414 |               0.0507759 |              0.9  |
| aroma_campaign_doe       | isoamyl_acetate | 0.00468664 |    0.0505482 |     -0.00839826 |   -0.00283999 |      0.0237955 |               0.077343  |              0.88 |
| aroma_campaign_doe       | ethyl_octanoate | 0.00509631 |    0.0584619 |     -0.0127983  |   -0.00454746 |      0.0318688 |               0.108173  |              0.85 |
| aroma_joint_campaign_doe | ethyl_acetate   | 0.00443687 |    0.0419447 |     -0.00443587 |   -0.00129469 |      0.0161414 |               0.0507759 |              0.9  |
| aroma_joint_campaign_doe | isoamyl_acetate | 0.00468664 |    0.0505482 |     -0.00839826 |   -0.00283999 |      0.0237955 |               0.077343  |              0.88 |
| aroma_joint_campaign_doe | ethyl_octanoate | 0.00509631 |    0.0584619 |     -0.0127983  |   -0.00454746 |      0.0318688 |               0.108173  |              0.85 |

## Cierre de masa aromatico

La tabla `tables/06_cierre_masa_aromatico.csv` resume el cierre por candidato y especie. La columna `mass_closure_rel_error` compara la integral numerica de sintesis con el incremento de `A_liq + A_loss`. La columna `loss_balance_rel_error` compara la integral numerica de perdida con el incremento de `A_loss`. Estos errores son diagnosticos de postproceso numerico; las ecuaciones del modelo imponen el cierre en la discretizacion.

Vista preliminar:

| workflow           | candidate                    | species         |   liquid_final |   loss_final |   condensate_final |   mass_closure_rel_error |   loss_balance_rel_error |
|:-------------------|:-----------------------------|:----------------|---------------:|-------------:|-------------------:|-------------------------:|-------------------------:|
| aroma_campaign_doe | cold_synthesis_hot_stripping | ethyl_acetate   |      42.4852   |    5.19289   |          4.6736    |             -0.00210507  |              -0.00595078 |
| aroma_campaign_doe | cold_synthesis_hot_stripping | ethyl_octanoate |       0.452972 |    0.0513013 |          0.0436061 |             -0.00235136  |              -0.00546047 |
| aroma_campaign_doe | cold_synthesis_hot_stripping | isoamyl_acetate |       1.8036   |    0.208659  |          0.18362   |             -0.00239743  |              -0.00578265 |
| aroma_campaign_doe | combined_stress_long_horizon | ethyl_acetate   |      37.087    |    3.38409   |          3.04568   |              0.000401911 |              -0.0013891  |
| aroma_campaign_doe | combined_stress_long_horizon | ethyl_octanoate |       0.403678 |    0.0259983 |          0.0220986 |              0.000160121 |              -0.00101903 |
| aroma_campaign_doe | combined_stress_long_horizon | isoamyl_acetate |       1.59496  |    0.120846  |          0.106345  |              0.000115083 |              -0.00120945 |
| aroma_campaign_doe | ethanol_partition_late_pulse | ethyl_acetate   |      35.2744   |    3.2707    |          2.94363   |              0.000211987 |              -0.00256079 |
| aroma_campaign_doe | ethanol_partition_late_pulse | ethyl_octanoate |       0.378157 |    0.0308918 |          0.026258  |             -9.25014e-05 |              -0.00205903 |
| aroma_campaign_doe | ethanol_partition_late_pulse | isoamyl_acetate |       1.50467  |    0.128632  |          0.113196  |             -0.00014924  |              -0.00233739 |
| aroma_campaign_doe | fructose_rich_iG_probe       | ethyl_acetate   |      44.193    |    5.21986   |          4.69788   |             -0.00161164  |              -0.00432456 |
| aroma_campaign_doe | fructose_rich_iG_probe       | ethyl_octanoate |       0.466067 |    0.0567782 |          0.0482614 |             -0.00185646  |              -0.00389901 |
| aroma_campaign_doe | fructose_rich_iG_probe       | isoamyl_acetate |       1.86711  |    0.219427  |          0.193096  |             -0.00190223  |              -0.00416161 |

## Sensibilidad y estimabilidad

El paquete incluye analisis FIM, reduccion de varianza, ranking DOE y diagnostico de eigen-direcciones. La tabla consolidada esta en `tables/07_sensibilidad_estimabilidad.csv`.

Vista preliminar:

| workflow              | source_file                                                                                      | analysis        | parameter   |      theta |   std_log_approx |   approx_95_multiplier |    fim_diag | active_bound   | classification     |   group |   campaign_var_ratio |   campaign_var_reduction |   design |   n_directions |   min_relative_eigenvalue |   condition_number |   numeric_rank |   effective_rank |
|:----------------------|:-------------------------------------------------------------------------------------------------|:----------------|:------------|-----------:|-----------------:|-----------------------:|------------:|:---------------|:-------------------|--------:|---------------------:|-------------------------:|---------:|---------------:|--------------------------:|-------------------:|---------------:|-----------------:|
| new_must_estimability | fermentation_model/results/new_must_glycerol_estimability_doe/parameter_estimability_summary.csv | natural_current | mu0         |  0.10091   |        0.0458004 |                1.09392 | 25536       | False          | well_estimated     |     nan |                  nan |                      nan |      nan |            nan |                       nan |                nan |            nan |              nan |
| new_must_estimability | fermentation_model/results/new_must_glycerol_estimability_doe/parameter_estimability_summary.csv | natural_current | sN          | 18         |        1.29125   |               12.5642  |     7.81023 | False          | weak_or_confounded |     nan |                  nan |                      nan |      nan |            nan |                       nan |                nan |            nan |              nan |
| new_must_estimability | fermentation_model/results/new_must_glycerol_estimability_doe/parameter_estimability_summary.csv | natural_current | qN          |  0.0289682 |        0.0457167 |                1.09374 | 18127.4     | False          | well_estimated     |     nan |                  nan |                      nan |      nan |            nan |                       nan |                nan |            nan |              nan |
| new_must_estimability | fermentation_model/results/new_must_glycerol_estimability_doe/parameter_estimability_summary.csv | natural_current | qXG         |  0.1125    |        1.51765   |               19.5815  |    11.6115  | False          | weak_or_confounded |     nan |                  nan |                      nan |      nan |            nan |                       nan |                nan |            nan |              nan |
| new_must_estimability | fermentation_model/results/new_must_glycerol_estimability_doe/parameter_estimability_summary.csv | natural_current | qXF         |  0.1125    |        1.37711   |               14.8668  |    12.2974  | False          | weak_or_confounded |     nan |                  nan |                      nan |      nan |            nan |                       nan |                nan |            nan |              nan |
| new_must_estimability | fermentation_model/results/new_must_glycerol_estimability_doe/parameter_estimability_summary.csv | natural_current | betaG0      |  0.502511  |        0.128106  |                1.28542 | 22812.6     | False          | well_estimated     |     nan |                  nan |                      nan |      nan |            nan |                       nan |                nan |            nan |              nan |
| new_must_estimability | fermentation_model/results/new_must_glycerol_estimability_doe/parameter_estimability_summary.csv | natural_current | sG          |  0.03      |        0.229029  |                1.56658 |   997.28    | False          | well_estimated     |     nan |                  nan |                      nan |      nan |            nan |                       nan |                nan |            nan |              nan |
| new_must_estimability | fermentation_model/results/new_must_glycerol_estimability_doe/parameter_estimability_summary.csv | natural_current | betaF0      |  0.487841  |        0.180151  |                1.42347 | 20813.3     | False          | well_estimated     |     nan |                  nan |                      nan |      nan |            nan |                       nan |                nan |            nan |              nan |
| new_must_estimability | fermentation_model/results/new_must_glycerol_estimability_doe/parameter_estimability_summary.csv | natural_current | sF          |  0.03      |        0.195403  |                1.46666 |  1060.99    | False          | well_estimated     |     nan |                  nan |                      nan |      nan |            nan |                       nan |                nan |            nan |              nan |
| new_must_estimability | fermentation_model/results/new_must_glycerol_estimability_doe/parameter_estimability_summary.csv | natural_current | qEG         |  1.17707   |        0.0888974 |                1.19034 | 32123.3     | False          | well_estimated     |     nan |                  nan |                      nan |      nan |            nan |                       nan |                nan |            nan |              nan |
| new_must_estimability | fermentation_model/results/new_must_glycerol_estimability_doe/parameter_estimability_summary.csv | natural_current | qEF         |  1.04431   |        0.135985  |                1.30543 | 22249.3     | False          | well_estimated     |     nan |                  nan |                      nan |      nan |            nan |                       nan |                nan |            nan |              nan |
| new_must_estimability | fermentation_model/results/new_must_glycerol_estimability_doe/parameter_estimability_summary.csv | natural_current | iG          |  0.0068611 |        0.336524  |                1.93399 |   594.362   | False          | well_estimated     |     nan |                  nan |                      nan |      nan |            nan |                       nan |                nan |            nan |              nan |

## Archivos de respaldo

- Datos fuente: `data_sources/`
- Codigo fuente clave: `code_sources/`
- Notebooks de trabajo: `notebooks/`
- Reportes originales generados por scripts: `source_reports/`
- Tablas originales de resultados: `source_tables/`
- Figuras de campanas y diagnosticos: `figures/`
- Inventario con hashes: `00_manifest.csv`

## Texto sugerido para informe

Se desarrollo y calibro un modelo dinamico de balance de masa para fermentacion que integra biomasa, nitrogeno, glucosa, fructosa, etanol y glicerol usando datos VL3 y mosto sintetico/natural. Sobre esta base se implemento una extension para perdidas aromaticas por stripping con CO2, particion gas-liquido derivada de UNIFAC, acumulacion de perdidas y prediccion de condensado terminal. El bundle A02 documenta ecuaciones, supuestos, datos fuente, parametros, metricas de ajuste, cierre de masa aromatico, version de codigo y analisis de sensibilidad/estimabilidad, cubriendo la evidencia minima requerida para la actividad OE2-2.13.

## Limitacion declarable

La evidencia cubre el desarrollo, calibracion base, simulacion y diseno experimental de perdidas aromaticas. Si la rendicion exige validacion experimental directa de gas o condensado, debe anexarse la medicion externa correspondiente como complemento; el presente bundle deja el punto trazado y preparado mediante la tabla de cierre de masa y el balance `A_cond = eta A_loss(t_final)`.
