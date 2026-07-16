# Pilot 2026 notebooks

Regenerate and execute the data-loading review from the repository root:

```bash
python fermentation_model/pilot_2026/notebooks/build_and_execute_data_loading_qc_notebook.py
```

The builder uses only the repository's existing Python dependencies and embeds
all tables and figures in `pilot_2026_data_loading_qc.ipynb`. The notebook is
also compatible with a standard Python 3 Jupyter kernel for interactive review.
