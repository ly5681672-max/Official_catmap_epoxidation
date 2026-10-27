"""Export product production rates from official CatMAP results without interpolation.

The existing model, reaction network, energies and numerical parameters are untouched.
The catalyst table solves at exact descriptor coordinates for every selected surface.
"""
import csv
from pathlib import Path
from math import isfinite


def export_tof_tables(model, output_dir, product="C6H12O_g"):
    """Export grid values and exact per-catalyst product production rates.

    Uses CatMAP production_rate (nonnegative part of net gas TOF), not a
    separately inferred rate law. Fails explicitly if a catalyst cannot solve.
    """
    output_dir = Path(output_dir)
    if len(model.descriptor_names) != 2:
        raise ValueError("This export expects exactly two descriptors")
    names = list(model.descriptor_names)
    labels = list(model.output_labels.get("production_rate", []))
    if product not in labels:
        raise ValueError(f"Product {product} not found in CatMAP production_rate labels: {labels}")
    idx = labels.index(product)
    rate_map = getattr(model, "production_rate_map", None)
    if not rate_map:
        raise ValueError("CatMAP production_rate_map is empty or missing")

    # All values directly from CatMAP's solved descriptor grid.
    grid_rows = []
    for coords, rates in rate_map:
        if len(coords) != 2 or idx >= len(rates):
            raise ValueError("Invalid production_rate_map entry")
        grid_rows.append([str(coords[0]), str(coords[1]), str(rates[idx])])

    # Exact CatMAP point evaluations, not grid interpolation or nearest-point lookup.
    selected_surfaces = list(model.surface_names)
    point_rows = []
    old_outputs = list(model.output_variables)
    try:
        # Keep coverage, which CatMAP needs to determine steady-state rates.
        model.output_variables = ["coverage", "production_rate"]
        for surf_i, surf in enumerate(selected_surfaces):
            coords = []
            for name in names:
                vals = model.species_definitions[name].get("formation_energy")
                if not isinstance(vals, (tuple, list)) or len(vals) != len(selected_surfaces):
                    raise ValueError(f"Invalid descriptor formation energies: {name}")
                val = float(vals[surf_i])
                if not isfinite(val):
                    raise ValueError(f"Non-finite {name} descriptor for {surf}")
                coords.append(val)
            try:
                model.mapper.get_point_output(coords)
            except Exception as exc:
                raise RuntimeError(f"CatMAP exact-point production rate FAILED for {surf}: {coords}") from exc
            current_labels = list(model.output_labels["production_rate"])
            current_idx = current_labels.index(product)
            point_rows.append([surf, *(str(c) for c in coords),
                               str(model._production_rate[current_idx])])
    finally:
        model.output_variables = old_outputs

    # Write files only after every direct point has succeeded; no invented zeros.
    grid_file = output_dir / f"TOF_{product}_grid.csv"
    catalyst_file = output_dir / f"TOF_{product}_catalysts.csv"
    for path, header, rows in [
        (grid_file, [*names, f"{product}_production_rate_per_site_s-1"], grid_rows),
        (catalyst_file, ["surface_name", *names,
                         f"{product}_production_rate_per_site_s-1"], point_rows),
    ]:
        with path.open("w", newline="", encoding="utf-8-sig") as f:
            writer = csv.writer(f)
            writer.writerow(header)
            writer.writerows(rows)
    print(f"CatMAP TOF exports: {grid_file}, {catalyst_file}")
    return grid_file, catalyst_file
