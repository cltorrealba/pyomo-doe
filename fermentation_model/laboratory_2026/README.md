# Laboratory 2026 campaign

This workspace separates two scientifically different streams that share the
same year and scale.

## Transferability/replication stream

The natural-must workbook `../data/Laboratorio 2026/mosto_natural_xthiol.xlsx`
contains laboratory homologues of pilot profiles. The master sheet currently
identifies `ING26-LAB004`, `ING26-LAB005`, `ING26-LAB006`, `ING26-LAB010`,
`ING26-LAB011` and `ING26-LAB012`. The relationship of the first three to the
2026 pilot campaign remains marked as unconfirmed in the registry.

## Optimal-design stream

The authoritative design is the volume-constrained nine-protocol campaign:

- Design handoff: `../results/design_execution_bundle_2026-06-09/`
- Rich execution copy: `../../db_20260609/` (`rec/` is current; `ref/` is historical)
- Raw Lot 1: `../data/Laboratorio 2026/DOE_Lote_1/`
- Raw Lot 2: `../data/Laboratorio 2026/DOE_Lote_2/`
- Raw Lot 3: expected at `../data/Laboratorio 2026/DOE_Lote_3/`
- Current notebooks: `notebooks/`

As of 2026-07-14, raw folders are present for `DOE26-F01` through `DOE26-F06`.
The analytical exports for Lot 2 and the entire Lot 3 source folder are still
pending in Git.

## Recommended workflow order

1. Inspect and homogenize new raw inputs without modifying them.
2. Run the local IPOPT calibration with an explicit run configuration.
3. Calculate weighted FIM, eigenvalues and weak eigendirections.
4. Run profile likelihood only from a verified local optimum.
5. Compare predicted versus realized information gain before selecting the next
   sequential experiment.

Profile likelihood is a frequentist likelihood-based method. Any Laplace or
Bayesian approximation must be reported separately rather than being grouped
under the profile-likelihood label.
