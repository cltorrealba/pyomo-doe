# Computational reproducibility

The objective is that a result can be traced to one Git commit, one immutable
input snapshot and one explicit model/solver configuration.

## Windows/IDAES baseline

From an Anaconda or Miniforge prompt:

```powershell
conda env create -f fermentation_model\environment.yml
conda activate fermentation-doe
idaes get-extensions
python -c "import pyomo.environ as pyo; print(pyo.SolverFactory('ipopt').available())"
```

`idaes get-extensions` installs the IDAES binary extensions, including the
supported IPOPT distribution. Record the working Windows environment after
validation:

```powershell
conda env export --from-history > fermentation_model\environment.from-history.yml
conda list --explicit > fermentation_model\environment.windows-explicit.txt
```

The two exported files should only be committed after confirming they come from
the machine that reproduces the reference run.

## Before every scientific run

1. Select campaign and experiment IDs from `campaigns/experiments.csv`.
2. Copy `config/local_ipopt_reference.json` to a run-specific configuration and
   fill every required field listed there.
3. Confirm raw-data integrity:

   ```powershell
   python fermentation_model\tools\campaign_audit.py --check-hashes
   ```

4. Capture the execution context:

   ```powershell
   python fermentation_model\tools\capture_run_context.py `
     --config path\to\run_config.json `
     --output path\to\results\run_manifest.json
   ```

5. Store solver logs, parameter estimates, objective contributions, residuals,
   FIM diagnostics and figures beside that manifest.

## Scientific acceptance gates

- The selected local optimum must be reproducible from the recorded initial
  guess and solver options.
- No profile-likelihood point may improve the baseline objective beyond the
  declared numerical tolerance. If it does, rebase the optimum first.
- Parameters at bounds must be reported explicitly.
- FIM and eigenvalue comparisons must use the same parameter scaling and
  measurement-error model.
- Calibration and validation experiment sets must be listed separately.
- Executed notebooks are narrative evidence; they are not the source of model
  functions and must never be imported by cell index.

## What is not authoritative

Particle-swarm and broad multistart files are retained for traceability but are
not part of the preferred scientific workflow. Administrative `AXX` bundles
are renditions/deliverables and are excluded from model provenance.
