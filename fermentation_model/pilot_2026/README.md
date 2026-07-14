# Pilot 2026 campaign

## Scope

Three pilot reactors were operated in parallel with different dynamic
temperature profiles and constant nutrition. Aroma characterization data have
recently arrived and the repository dataset is still under construction.

## Current source

- GC report: `../data/Piloto 2026/Lotes1a3/Resultados_GC.csv`

The file is an exported laboratory report rather than a tidy analytical table.
It must remain unchanged. A processed representation should be written to a
future campaign-specific results directory, never over this source file.

## Missing registry information

Before calibration or comparison, reconstruct one row per pilot fermentation
with:

- run/reactor identifier;
- lot or must batch;
- inoculation timestamp;
- temperature-profile identifier;
- confirmation that nutrition was constant;
- sampling timestamps and condensate fraction boundaries;
- mapping between GC sample identifiers and reactor/time;
- deviations or incidents.

Until that mapping exists, `PILOT26-LOTES1A3` remains a composite placeholder
in `../campaigns/experiments.csv`.
