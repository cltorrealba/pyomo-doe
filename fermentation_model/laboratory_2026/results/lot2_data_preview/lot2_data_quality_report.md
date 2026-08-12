# Lot 2 data ingestion and QC

## Scope

This artifact integrates Y15 chemistry, manual ethanol entries, Oculyze biomass, online temperature, and online CO2 for the three synthetic-must MBDoE fermentations in Lot 2.

## Main conventions

- Time zero is the first Oculyze sample in each reactor.
- Total Oculyze concentration is interpreted as million cells/mL. Biomass uses 30 pg/cell.
- Negative low-range analytical readings are preserved in raw columns and censored at zero only in explicitly named model-safe columns.
- Ethanol is recomputed from the entered replicates as `% v/v * 7.89 = g/L`; cached Excel formula cells are ignored.
- Sampling removes 50 mL when ethanol was measured and 5 mL for another sample with analytical evidence. Planned rows without any measurement evidence remove 0 mL. The CO2 flow is normalized back to the initial 2 L volume.
- Planned input events are not silently treated as executed. The F1 fructose pulse and F3 biomass step have trajectory evidence; nitrogen pulse times remain protocol-based and unverified.

## QC summary

| process   |   n_sample_rows |   n_samples_with_observation_evidence |   n_y15_samples |   n_ethanol |   n_ethanol_used_after_qc |   n_biomass |   n_time_inferred |   n_negative_censored |   final_sampling_only_volume_ml |   temperature_first_h |   temperature_last_h |   co2_first_h |   co2_last_h |   chemistry_last_h |
|:----------|----------------:|--------------------------------------:|----------------:|------------:|--------------------------:|------------:|------------------:|----------------------:|--------------------------------:|----------------------:|---------------------:|--------------:|-------------:|-------------------:|
| F1        |              18 |                                    18 |              15 |          13 |                        13 |          17 |                 2 |                     4 |                            1325 |              0.666667 |              99.8333 |      0.666667 |      99.8333 |                190 |
| F2        |              18 |                                    18 |              15 |          13 |                        12 |          17 |                 2 |                     4 |                            1325 |              0.666667 |              99.8333 |      0.833333 |      99.8333 |                190 |
| F3        |              18 |                                    16 |              14 |          11 |                        11 |          15 |                 1 |                     5 |                            1425 |              0.833333 |              99.8333 |      0.833333 |      99.8333 |                168 |

## Repaired Oculyze records

| process   | sample_id_raw   | sample_id     | raw_datetime        | used_datetime       | repair_reason                       | record_changed   |
|:----------|:----------------|:--------------|:--------------------|:--------------------|:------------------------------------|:-----------------|
| F1        | DOE-LAB004-16   | DOE-LAB004-16 | 2026-07-13 12:00:00 | 2026-07-13 12:00:00 | duplicate_time_context              | False            |
| F1        | DOE-LAB004-17   | DOE-LAB004-17 | 2026-07-13 12:00:00 | 2026-07-13 16:30:00 | duplicate_description_time_repaired | True             |
| F2        | DOE-LAB005-1    | DOE-LAB005-9  | 2026-07-08 14:30:00 | 2026-07-08 14:30:00 | sequence_id_repair                  | True             |
| F2        | DOE-LAB005-16   | DOE-LAB005-16 | 2026-07-13 12:00:00 | 2026-07-13 12:00:00 | duplicate_time_context              | False            |
| F2        | DOE-LAB005-17   | DOE-LAB005-17 | 2026-07-13 12:00:00 | 2026-07-13 16:30:00 | duplicate_description_time_repaired | True             |
| F3        | DOE-LAB006-5    | DOE-LAB006-15 | 2026-07-10 12:30:00 | 2026-07-10 12:30:00 | sequence_id_repair                  | True             |

## Ethanol points held for review

The following points fail a deliberately conservative sugar-to-ethanol screening rule. They remain visible in plots and processed data but are excluded from the default calibration table until sample timing or measurement identity is confirmed.

| process   | sample_id     |   t_h |   sugar_total_g_l |   ethanol_g_l_for_model | ethanol_mass_balance_screen        |
|:----------|:--------------|------:|------------------:|------------------------:|:-----------------------------------|
| F2        | DOE-LAB005-16 |   168 |             70.01 |                 94.8378 | review_time_mapping_or_measurement |

## Interpretation boundary

The online process files end near 99 h, whereas chemistry continues to approximately 190 h. Cumulative CO2 after the last online record is therefore a lower bound, not a full-batch mass balance.
