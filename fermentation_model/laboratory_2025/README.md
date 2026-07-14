# Laboratory 2025 campaign

## Scope

Synthetic-must, two-factor CCD campaign. The controlled factors are fixed
isothermal temperature and nutrient-addition timing. Ten fermentations are
recorded (`MS007`-`MS016`).

## Source of truth

- Raw/homologated workbook: `../data/Laboratorio 2025/mosto_sintetico_vl3.xlsx`
- Experiment registry: `../campaigns/experiments.csv`

The workbook is treated as an immutable source. Its `Diseño_CCD`,
`Datos_homologados`, `Datos_originales`, `Flags_calidad` and mapping sheets
must remain together.

## Scientific use

This campaign supplies the synthetic-must calibration and estimability prior
used by the shared fermentation model. It is distinct from the nine-protocol
model-based DOE executed in 2026.

Shared active modules are documented in `../shared/README.md`. Historical
calibration notebooks are retained under
`../legacy/development_2026/notebooks/initial_calibration/`.

The active shared calibration and estimability outputs are under
`../shared/results/new_must_glycerol_estimability_doe/`.
