"""Descriptor-driven figure directory and range setup for CatMAP tutorials.

Only sets descriptor_ranges when the .mkm setup leaves it as None.
DFT/NEB energies, reaction network and scaling constraints are not modified.
"""
from math import isfinite
from pathlib import Path
import re


def configure_descriptor_output(model, script_file, margin=0.05):
    """Create a folder named after descriptor_names; infer missing ranges from input data."""
    names = list(model.descriptor_names)
    if not names:
        raise ValueError("CatMAP descriptor_names is empty")

    if getattr(model, "descriptor_ranges", None) is None:
        if not 0 <= margin < 1:
            raise ValueError("margin must be in [0, 1)")
        bounds = []
        for name in names:
            definition = model.species_definitions.get(name)
            if definition is None:
                raise ValueError(f"No CatMAP species definition for descriptor {name}")
            energies = definition.get("formation_energy")
            if not isinstance(energies, (list, tuple)) or len(energies) != len(model.surface_names):
                raise ValueError(f"{name}: formation_energy must contain one value per selected surface")

            values = []
            for surface, energy in zip(model.surface_names, energies):
                try:
                    value = float(energy)
                except (TypeError, ValueError) as exc:
                    raise ValueError(f"{surface}: missing or invalid {name} descriptor energy") from exc
                if not isfinite(value):
                    raise ValueError(f"{surface}: non-finite {name} descriptor energy")
                values.append(value)

            low, high = min(values), max(values)
            if low >= high:
                raise ValueError(f"{name}: data range has zero span; refusing to invent a plotting range")
            pad = margin * (high - low)
            bounds.append([low - pad, high + pad])
        model.descriptor_ranges = bounds

    # Only filesystem-safe characters; do not change CatMAP's species identifiers.
    tags = [re.sub(r"[^A-Za-z0-9_-]", "_", name.removesuffix("_s")) for name in names]
    if any(not tag or tag in (".", "..") for tag in tags):
        raise ValueError("A descriptor cannot be converted to a safe directory name")
    output_dir = Path(script_file).resolve().parent / ("descriptor_" + "_".join(tags))
    output_dir.mkdir(parents=True, exist_ok=True)

    print(f"Descriptors: {names}")
    print(f"Descriptor ranges: {model.descriptor_ranges}")
    print(f"Figure directory: {output_dir}")
    return output_dir
