"""Pilot 2026 adaptive model-based design workflow.

Campaign-specific data adapters and gate-controlled runners live here.  The
package deliberately does not expose an executable fermentation schedule:
calibration, safety bounds, and owner approval are separate hard gates.
"""

from .build_model_dataset import AdapterConfig, ModelDataset, load_model_dataset

__all__ = ["AdapterConfig", "ModelDataset", "load_model_dataset"]
