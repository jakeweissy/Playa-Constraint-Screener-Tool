#!/usr/bin/env python3
"""
playa_screen.py
Playa-Aware Environmental Constraint Screening Tool

A repeatable GIS workflow that screens project boundaries (and candidate /
prospecting areas) against playa-centered environmental constraints and produces,
per site:
  (1) a composite constraint score (0-100) and rating,
  (2) environmental permitting risk flags,
  (3) ESG / reputational exposure,
  (4) operational + hydrologic risk.

The composite is built from the FIVE constraint dimensions in the project plan:
  C1  ESG / Reputation
  C2  Hydrology & Drainage
  C3  Migratory Bird Exposure
  C4  CWA Section 404 (WOTUS, scenario-adjusted)
  C5  PLJV Siting Guidance

Workflow: load -> reproject -> fix geometry -> clip -> intersect/buffer -> score -> export.

IMPORTANT: This is a SCREENING tool, not a jurisdictional determination or a wetland
delineation. Only the U.S. Army Corps of Engineers issues jurisdictional determinations,
and a field delineation by a qualified wetland scientist is still required for permitting.

Usage example:
  python playa_screen.py \
      --boundary data/project_boundaries.geojson \
      --playas   data/pljv_probable_playas_v5.shp \
      --clusters data/pljv_playa_clusters.shp \
      --nwi      data/nwi_wetlands.shp \
      --floodplain data/fema_nfhl.shp \
      --critical-habitat data/usfws_critical_habitat.shp \
      --id-field site_name \
      --scenario conforming \
      --outdir outputs

Dependencies: geopandas, pandas, shapely (>=2.0 recommended). See requirements.txt.
"""

from __future__ import annotations

import argparse
import os
import sys
import json
import shutil
import warnings
from dataclasses import dataclass, field
from datetime import date
from string import Template

import pandas as pd
import geopandas as gpd

# Square meters -> hectares / acres
M2_TO_HA = 1.0 / 10_000.0
M2_TO_ACRE = 0.000247105


# --------------------------------------------------------------------------- #
# Configuration
# --------------------------------------------------------------------------- #
@dataclass
class Config:
    # Equal-area projected CRS for valid area math. EPSG:5070 = NAD83 CONUS Albers (m).
    # Alternative for Texas-only work: EPSG:3083 (NAD83 Texas Centric Albers Equal Area).
    target_crs: str = "EPSG:5070"
    assume_crs: str = "EPSG:4326"  # used only if an input layer is missing a CRS

    scenario: str = "conforming"   # "conforming" (conservative) | "proposed" (Nov-2025 rule)
    buffer_m: float = 100.0        # avoidance buffer applied around playas
    nrhp_buffer_m: float = 1000.0  # proximity radius for NRHP points
    large_playa_ha: float = 4.0    # threshold for a "large" playa (PLJV avoidance guidance)
    large_playa_cap: int = 5       # # of large playas that maxes the PLJV sub-score

    # Derive cluster overlap from a column in the playa layer (e.g. PLJV "cluster")
    # instead of a separate clusters file. A playa counts as clustered when
    # its value in this column is greater than cluster_min.
    cluster_field: str = None
    cluster_min: float = 0.0
    # Geometry simplification (degrees) for the embedded web map. Smaller = more
    # detail/bigger file; 0.003 ~ 300 m, plenty for a statewide overview.
    map_simplify_tol: float = 0.003

    # Caps used to normalize percentage-based indicators to 0-100.
    playa_pct_cap: float = 10.0
    buffer_pct_cap: float = 25.0
    floodplain_pct_cap: float = 25.0
    nwi_pct_cap: float = 10.0
    # Share of the site occupied by clustered playas that maxes cluster intensity.
    # Kept small because playas are small features; tune up to be less sensitive.
    cluster_pct_cap: float = 5.0

    # Composite weights across the five constraint dimensions. Renormalized over
    # whichever dimensions can actually be assessed given the supplied layers.
    weights: dict = field(default_factory=lambda: {
        "esg": 0.25,
        "hydrology": 0.20,
        "migratory": 0.20,
        "cwa404": 0.15,
        "pljv": 0.20,
    })

    # WOTUS scenario factor applied to the §404 dimension. Under the proposed
    # rule most playas are less likely jurisdictional, but residual risk remains.
    scenario_404_factor: dict = field(default_factory=lambda: {
        "conforming": 1.0,
        "proposed": 0.5,
    })

    flag_threshold: float = 50.0   # a dimension at/above this raises a flag
    rating_bands: tuple = (
        (0, 25, "Low"),
        (25, 50, "Moderate"),
        (50, 75, "High"),
        (75, 101, "Severe"),
    )
    rating_colors: dict = field(default_factory=lambda: {
        "Low": "#2e7d32", "Moderate": "#f9a825",
        "High": "#ef6c00", "Severe": "#c62828",
    })


# --------------------------------------------------------------------------- #
# Geometry helpers
# --------------------------------------------------------------------------- #
def _union(geoseries):
    """Union a GeoSeries across geopandas versions."""
    try:
        return geoseries.union_all()
    except AttributeError:
        return geoseries.unary_union


def fix_geoms(gdf):
    """Repair invalid geometries; drop empties/nulls."""
    if gdf is None or gdf.empty:
        return gdf
    geom = gdf.geometry
    try:
        gdf = gdf.set_geometry(geom.make_valid())
    except Exception:
        gdf = gdf.set_geometry(geom.buffer(0))
    gdf = gdf[gdf.geometry.notna() & ~gdf.geometry.is_empty]
    return gdf


def load_layer(path, cfg, name, layer=None):
    """Read a vector layer, set/repair CRS, reproject to the analysis CRS.

    `layer` selects a specific layer inside a multi-layer source (e.g. a
    GeoPackage that contains more than one table)."""
    if not path:
        return None
    if not os.path.exists(path):
        warnings.warn(f"[{name}] file not found, skipping: {path}")
        return None
    gdf = gpd.read_file(path, layer=layer) if layer else gpd.read_file(path)
    if gdf.empty:
        warnings.warn(f"[{name}] layer is empty: {path}")
        return None
    if gdf.crs is None:
        warnings.warn(f"[{name}] missing CRS; assuming {cfg.assume_crs}. Verify this!")
        gdf = gdf.set_crs(cfg.assume_crs)
    gdf = gdf.to_crs(cfg.target_crs)
    return fix_geoms(gdf)


def clip_layer(gdf, mask_geom):
    """Clip a layer to a single mask geometry; return None if nothing remains."""
    if gdf is None or gdf.empty:
        return None
    try:
        clipped = gpd.clip(gdf, mask_geom)
    except Exception:
        sub = gdf[gdf.intersects(mask_geom)].copy()
        if sub.empty:
            return None
        sub["geometry"] = sub.geometry.intersection(mask_geom)
        clipped = sub
    clipped = clipped[clipped.geometry.notna() & ~clipped.geometry.is_empty]
    return clipped if not clipped.empty else None


def _pct(part_ha, whole_ha):
    return 0.0 if whole_ha <= 0 else 100.0 * part_ha / whole_ha


def cap_ratio(value, cap):
    """Normalize a value to a 0-100 sub-score, capped."""
    if cap <= 0:
        return 0.0
    return max(0.0, min(value / cap, 1.0)) * 100.0


# --------------------------------------------------------------------------- #
# Per-site metrics
# --------------------------------------------------------------------------- #
def compute_metrics(site_geom, layers, cfg):
    """Compute all spatial metrics for one site geometry."""
    m = {"boundary_area_ha": site_geom.area * M2_TO_HA,
         "boundary_area_ac": site_geom.area * M2_TO_ACRE}
    b_ha = m["boundary_area_ha"]

    # --- Playas ---------------------------------------------------------- #
    m["playas_provided"] = layers.get("playas") is not None
    if m["playas_provided"]:
        inside = clip_layer(layers["playas"], site_geom)
        if inside is not None:
            areas_ha = inside.geometry.area * M2_TO_HA
            m["playa_count"] = int(len(inside))
            m["playa_area_ha"] = float(areas_ha.sum())
            m["large_playa_count"] = int((areas_ha >= cfg.large_playa_ha).sum())
        else:
            m["playa_count"] = 0
            m["playa_area_ha"] = 0.0
            m["large_playa_count"] = 0
        m["playa_pct"] = _pct(m["playa_area_ha"], b_ha)

        # Buffered avoidance footprint (playas near the site, buffered, clipped to site)
        near = clip_layer(layers["playas"], site_geom.buffer(cfg.buffer_m))
        if near is not None:
            avoid = _union(near.buffer(cfg.buffer_m)).intersection(site_geom)
            m["playa_buffer_area_ha"] = float(avoid.area * M2_TO_HA)
        else:
            m["playa_buffer_area_ha"] = 0.0
        m["playa_buffer_pct"] = _pct(m["playa_buffer_area_ha"], b_ha)

    # --- Playa clusters -------------------------------------------------- #
    # Either a dedicated clusters layer, OR derived from a column in the playa
    # layer (e.g. PLJV's "cluster" field) when --cluster-field is set.
    # We record the clustered AREA and its share of the site so intensity can be
    # graded (scale-invariant) rather than a saturating yes/no.
    has_cluster_layer = layers.get("playa_clusters") is not None
    derive_from_field = bool(cfg.cluster_field) and layers.get("playas") is not None
    m["clusters_provided"] = has_cluster_layer or derive_from_field
    if has_cluster_layer:
        cl = clip_layer(layers["playa_clusters"], site_geom)
        area_ha = float((cl.geometry.area * M2_TO_HA).sum()) if cl is not None else 0.0
        m["cluster_overlap"] = cl is not None
        m["cluster_area_ha"] = area_ha
        m["cluster_pct"] = _pct(area_ha, b_ha)
        m["cluster_member_count"] = None  # not meaningful for a polygon layer
    elif derive_from_field:
        inside = clip_layer(layers["playas"], site_geom)
        if inside is not None and cfg.cluster_field in inside.columns:
            vals = pd.to_numeric(inside[cfg.cluster_field], errors="coerce")
            clustered = inside[vals > cfg.cluster_min]
            area_ha = float((clustered.geometry.area * M2_TO_HA).sum())
            m["cluster_overlap"] = len(clustered) > 0
            m["cluster_member_count"] = int(len(clustered))
            m["cluster_area_ha"] = area_ha
            m["cluster_pct"] = _pct(area_ha, b_ha)
        else:
            # Field missing on this layer -> can't assess after all.
            m["clusters_provided"] = False
            if inside is not None:
                warnings.warn(f"cluster field '{cfg.cluster_field}' not found in playa layer")

    # --- NWI wetlands ---------------------------------------------------- #
    m["nwi_provided"] = layers.get("nwi") is not None
    if m["nwi_provided"]:
        nwi = clip_layer(layers["nwi"], site_geom)
        m["nwi_area_ha"] = float((nwi.geometry.area * M2_TO_HA).sum()) if nwi is not None else 0.0
        m["nwi_pct"] = _pct(m["nwi_area_ha"], b_ha)

    # --- FEMA floodplain ------------------------------------------------- #
    m["floodplain_provided"] = layers.get("floodplain") is not None
    if m["floodplain_provided"]:
        fp = clip_layer(layers["floodplain"], site_geom)
        m["floodplain_area_ha"] = float((fp.geometry.area * M2_TO_HA).sum()) if fp is not None else 0.0
        m["floodplain_pct"] = _pct(m["floodplain_area_ha"], b_ha)

    # --- USFWS critical habitat ----------------------------------------- #
    m["habitat_provided"] = layers.get("critical_habitat") is not None
    if m["habitat_provided"]:
        ch = clip_layer(layers["critical_habitat"], site_geom)
        m["critical_habitat_overlap"] = ch is not None
        m["critical_habitat_area_ha"] = float((ch.geometry.area * M2_TO_HA).sum()) if ch is not None else 0.0

    # --- NRHP proximity (points) ---------------------------------------- #
    m["nrhp_provided"] = layers.get("nrhp") is not None
    if m["nrhp_provided"]:
        near = clip_layer(layers["nrhp"], site_geom.buffer(cfg.nrhp_buffer_m))
        m["nrhp_nearby_count"] = int(len(near)) if near is not None else 0

    return m


# --------------------------------------------------------------------------- #
# Dimension metadata: single source of truth (order sets the displayed number)
# key -> (number, short label, plain-language description, recommended action)
# --------------------------------------------------------------------------- #
DIM_META = {
    "esg": (1, "ESG / Reputation",
            "How visible the playa impact is for biodiversity/ESG commitments — driven by "
            "playa coverage on the site and clustered-playa intensity.",
            "Document playa avoidance for ESG reporting; align with post-acquisition "
            "biodiversity commitments."),
    "hydrology": (2, "Hydrology & Drainage",
            "Flooding and operational risk from building in playa basins and their "
            "avoidance buffers (affects panels, BESS, and access roads).",
            "Run a site drainage/catchment review; keep infrastructure out of playa basins "
            "and buffers."),
    "migratory": (3, "Migratory Bird Exposure",
            "Central Flyway exposure graded by clustered-playa intensity and any "
            "listed-species critical habitat.",
            "Engage USFWS/IPaC early; consider seasonal timing and layout to reduce Central "
            "Flyway impacts."),
    "cwa404": (4, "CWA Section 404",
            "Residual federal wetland-permitting risk, adjusted for the chosen WOTUS "
            "scenario (conforming vs. proposed).",
            "Request a USACE preliminary/approved jurisdictional determination; budget for a "
            "field wetland delineation."),
    "pljv": (5, "PLJV Siting Guidance",
            "Conflict with PLJV guidance to avoid large isolated playas and playa clusters.",
            "Follow PLJV guidance: redesign to avoid large isolated playas and playa "
            "clusters."),
}


# --------------------------------------------------------------------------- #
# Scoring: the five constraint dimensions
# --------------------------------------------------------------------------- #
def score_constraints(m, cfg):
    """Return dimension sub-scores, which were assessable, and the composite."""
    playa_present = m.get("playa_area_ha", 0.0) > 0 or m.get("playa_count", 0) > 0
    habitat = bool(m.get("critical_habitat_overlap", False))

    # Graded cluster intensity (0-100): share of the site occupied by clustered
    # playas, normalized to a cap. Scale-invariant, so it discriminates whether the
    # unit is a small parcel or a whole county, instead of saturating to 100 on any
    # cluster touch. Stored back on m for the report.
    cluster_intensity = 0.0
    if m.get("clusters_provided"):
        cluster_intensity = cap_ratio(m.get("cluster_pct", 0.0), cfg.cluster_pct_cap)
    m["cluster_intensity"] = round(cluster_intensity, 1)
    ci = cluster_intensity

    dims, assessed = {}, {}

    # 1 ESG / Reputation: visible playa presence + clustered-playa intensity
    assessed["esg"] = m.get("playas_provided", False)
    if assessed["esg"]:
        dims["esg"] = max(
            cap_ratio(m.get("playa_pct", 0.0), cfg.playa_pct_cap),
            ci,
        )

    # 2 Hydrology & Drainage: footprint in playas / catchment + floodplain
    assessed["hydrology"] = m.get("playas_provided", False) or m.get("floodplain_provided", False)
    if assessed["hydrology"]:
        dims["hydrology"] = max(
            cap_ratio(m.get("playa_pct", 0.0), cfg.playa_pct_cap),
            cap_ratio(m.get("playa_buffer_pct", 0.0), cfg.buffer_pct_cap),
            cap_ratio(m.get("floodplain_pct", 0.0), cfg.floodplain_pct_cap),
        )

    # 3 Migratory Bird Exposure: Central Flyway cluster intensity + listed-species habitat
    assessed["migratory"] = m.get("clusters_provided", False) or m.get("habitat_provided", False)
    if assessed["migratory"]:
        dims["migratory"] = max(
            ci,
            100.0 if habitat else 0.0,
            cap_ratio(m.get("playa_pct", 0.0), cfg.playa_pct_cap) * 0.5,
        )

    # 4 CWA Section 404: residual jurisdictional risk, WOTUS-scenario-adjusted
    assessed["cwa404"] = m.get("nwi_provided", False) or m.get("playas_provided", False)
    if assessed["cwa404"]:
        base = max(cap_ratio(m.get("nwi_pct", 0.0), cfg.nwi_pct_cap),
                   50.0 if playa_present else 0.0)
        factor = cfg.scenario_404_factor.get(cfg.scenario, 1.0)
        dims["cwa404"] = base * factor

    # 5 PLJV Guidance: avoid large isolated playas AND clusters
    assessed["pljv"] = m.get("playas_provided", False)
    if assessed["pljv"]:
        dims["pljv"] = max(
            ci,
            cap_ratio(m.get("large_playa_count", 0), cfg.large_playa_cap),
            cap_ratio(m.get("playa_pct", 0.0), cfg.playa_pct_cap) * 0.6,
        )

    # Composite: weighted mean over assessed dimensions (weights renormalized)
    w = {k: cfg.weights[k] for k in dims if assessed.get(k)}
    wsum = sum(w.values())
    composite = sum(dims[k] * w[k] for k in w) / wsum if wsum > 0 else 0.0

    rating = next((lbl for lo, hi, lbl in cfg.rating_bands if lo <= composite < hi), "Low")
    not_assessed = [k for k, ok in assessed.items() if not ok]
    return {
        "dimensions": dims,
        "not_assessed": not_assessed,
        "composite": round(composite, 1),
        "rating": rating,
    }


def build_flags(m, scoring, cfg):
    """Human-readable flags + recommended actions for dimensions over threshold."""
    flags = []
    for k, v in scoring["dimensions"].items():
        if v >= cfg.flag_threshold:
            num, label, _desc, action = DIM_META[k]
            flags.append({"dimension": f"{num}. {label}", "score": round(v, 1),
                          "action": action})
    return sorted(flags, key=lambda f: -f["score"])


# --------------------------------------------------------------------------- #
# Reporting
# --------------------------------------------------------------------------- #
def _fmt(v, nd=1):
    return f"{v:,.{nd}f}" if isinstance(v, (int, float)) else str(v)


def render_report(site_id, m, scoring, flags, cfg, template_path, out_path):
    """Fill the one-page HTML report template for a single site."""
    # Cluster value: graded intensity rather than a yes/no.
    if not m.get("clusters_provided"):
        cluster_val = "not assessed"
    elif m.get("cluster_overlap"):
        intensity = m.get("cluster_intensity", 0)
        cnt = m.get("cluster_member_count")
        extra = (f"{cnt} clustered playas, " if cnt is not None else "")
        cluster_val = f"intensity {intensity:.0f}/100 ({extra}{_fmt(m.get('cluster_pct', 0))}% of site)"
    else:
        cluster_val = "none in site"

    metric_rows = []
    rows = [
        ("Boundary area", f"{_fmt(m.get('boundary_area_ha', 0))} ha "
                          f"({_fmt(m.get('boundary_area_ac', 0))} ac)", ""),
        ("Playas in boundary", f"{m.get('playa_count', '-')} "
                               f"({_fmt(m.get('playa_area_ha', 0))} ha, "
                               f"{_fmt(m.get('playa_pct', 0))}% of site)",
         "Count and total area of probable playas whose footprint falls inside the site."),
        (f"Large playas (&ge; {cfg.large_playa_ha:g} ha)", m.get("large_playa_count", "-"),
         f"Playas at least {cfg.large_playa_ha:g} hectares in area. PLJV specifically flags "
         f"large isolated playas for avoidance, so these carry extra weight."),
        ("Playa avoidance footprint", f"{_fmt(m.get('playa_buffer_area_ha', 0))} ha "
                                      f"({_fmt(m.get('playa_buffer_pct', 0))}% of site)",
         f"Site area within {cfg.buffer_m:g} m of any playa &mdash; a proxy for the land you "
         f"would likely keep clear of panels, BESS, and access roads. Larger = more design "
         f"constraint."),
        ("Playa cluster intensity", cluster_val,
         f"Graded by the share of the site occupied by clustered playas (maxes at "
         f"{cfg.cluster_pct_cap:g}% of site), not a yes/no. Clusters are PLJV's highest "
         f"waterfowl-avoidance priority."),
    ]
    # Optional layers: only shown when actually supplied (keeps the report uncluttered).
    if m.get("nwi_provided"):
        rows.append(("NWI wetland area", f"{_fmt(m.get('nwi_area_ha', 0))} ha",
                     "Mapped wetlands (National Wetlands Inventory) within the site."))
    if m.get("floodplain_provided"):
        rows.append(("FEMA floodplain", f"{_fmt(m.get('floodplain_pct', 0))}% of site",
                     "Share of the site in a FEMA mapped flood hazard area."))
    if m.get("habitat_provided"):
        rows.append(("USFWS critical habitat",
                     "Overlap" if m.get("critical_habitat_overlap") else "No overlap",
                     "Designated critical habitat for ESA-listed species intersecting the site."))
    if m.get("nrhp_provided"):
        rows.append(("NRHP sites nearby", m.get("nrhp_nearby_count", 0),
                     "Public National Register historic sites within the proximity radius."))
    for label, val, hint in rows:
        h = f"<div class='hint'>{hint}</div>" if hint else ""
        metric_rows.append(f"<tr><td>{label}{h}</td><td>{val}</td></tr>")

    dim_rows = []
    for k, (num, label, desc, _action) in sorted(DIM_META.items(), key=lambda kv: kv[1][0]):
        name = f"{num}. {label}"
        if k in scoring["dimensions"]:
            val, cls = f"{scoring['dimensions'][k]:.1f}", ""
        else:
            val, cls = "not assessed", " class='na'"
        dim_rows.append(
            f"<tr><td><b>{name}</b><div class='desc'>{desc}</div></td><td{cls}>{val}</td></tr>")

    if flags:
        flag_html = "".join(
            f"<li><b>{f['dimension']}</b> (score {f['score']}): {f['action']}</li>"
            for f in flags)
        flag_html = f"<ul>{flag_html}</ul>"
    else:
        flag_html = "<p>No constraint dimension exceeded the flag threshold.</p>"

    with open(template_path, "r", encoding="utf-8") as fh:
        tmpl = Template(fh.read())

    html = tmpl.safe_substitute(
        site_name=str(site_id),
        run_date=date.today().isoformat(),
        scenario=cfg.scenario,
        crs=cfg.target_crs,
        score=f"{scoring['composite']:.1f}",
        rating=scoring["rating"],
        rating_color=cfg.rating_colors.get(scoring["rating"], "#555"),
        metric_rows="".join(metric_rows),
        dim_rows="".join(dim_rows),
        flags=flag_html,
    )
    with open(out_path, "w", encoding="utf-8") as fh:
        fh.write(html)
    return out_path


# --------------------------------------------------------------------------- #
# Driver
# --------------------------------------------------------------------------- #
def run(args):
    cfg = Config(target_crs=args.crs, scenario=args.scenario, buffer_m=args.buffer,
                 cluster_field=args.cluster_field, cluster_min=args.cluster_min,
                 cluster_pct_cap=args.cluster_pct_cap)
    os.makedirs(args.outdir, exist_ok=True)

    # Find the report template even if the user runs from another folder.
    script_dir = os.path.dirname(os.path.abspath(__file__))
    template_path = args.template
    if not os.path.exists(template_path):
        alt = os.path.join(script_dir, os.path.basename(args.template))
        if os.path.exists(alt):
            template_path = alt
        else:
            sys.exit(f"ERROR: report template not found: {args.template}\n"
                     f"Keep report_template.html next to playa_screen.py, or pass "
                     f"--template /full/path/to/report_template.html")

    boundary = load_layer(args.boundary, cfg, "boundary", layer=args.boundary_layer)
    if boundary is None or boundary.empty:
        sys.exit("ERROR: boundary layer could not be loaded.")

    layers = {
        "playas": load_layer(args.playas, cfg, "playas"),
        "playa_clusters": load_layer(args.clusters, cfg, "playa_clusters"),
        "nwi": load_layer(args.nwi, cfg, "nwi"),
        "floodplain": load_layer(args.floodplain, cfg, "floodplain"),
        "critical_habitat": load_layer(args.critical_habitat, cfg, "critical_habitat"),
        "nrhp": load_layer(args.nrhp, cfg, "nrhp"),
    }

    id_field = args.id_field if (args.id_field and args.id_field in boundary.columns) else None
    results, scored_geoms = [], []

    for idx, row in boundary.iterrows():
        site_id = row[id_field] if id_field else f"site_{idx}"
        geom = row.geometry
        if geom is None or geom.is_empty:
            warnings.warn(f"Skipping {site_id}: empty geometry")
            continue

        m = compute_metrics(geom, layers, cfg)
        scoring = score_constraints(m, cfg)
        flags = build_flags(m, scoring, cfg)

        safe_id = "".join(c if c.isalnum() else "_" for c in str(site_id))
        report_path = os.path.join(args.outdir, f"report_{safe_id}.html")
        render_report(site_id, m, scoring, flags, cfg, template_path, report_path)

        flat = {"site_id": site_id, "composite_score": scoring["composite"],
                "rating": scoring["rating"], "scenario": cfg.scenario}
        flat.update({k: round(v, 1) for k, v in scoring["dimensions"].items()})
        for key in ("playa_count", "playa_area_ha", "playa_pct", "large_playa_count",
                    "playa_buffer_pct", "cluster_overlap", "cluster_intensity",
                    "cluster_member_count", "cluster_pct", "nwi_pct",
                    "floodplain_pct", "critical_habitat_overlap", "boundary_area_ha"):
            if key in m:
                flat[key] = m[key]
        flat["flags"] = "; ".join(f["dimension"] for f in flags)
        flat["not_assessed"] = ";".join(scoring["not_assessed"])
        results.append(flat)

        # Build the output feature from ONLY clean fields we control. Copying the
        # boundary's original attributes can drag in NaN/odd types that produce
        # invalid JSON (e.g. literal NaN tokens) the web map can't parse.
        props = {
            "site_name": str(site_id),
            "composite_score": scoring["composite"],
            "rating": scoring["rating"],
            "report_file": f"report_{safe_id}.html",
            "cluster_intensity": float(m.get("cluster_intensity", 0) or 0),
            "flags": flat["flags"],
        }
        for dk, dv in scoring["dimensions"].items():
            props[f"dim_{dk}"] = round(dv, 1)
        g = gpd.GeoDataFrame([props], geometry=[geom], crs=cfg.target_crs)
        scored_geoms.append(g)

        print(f"  {site_id}: {scoring['composite']:.1f} ({scoring['rating']})")

    df = pd.DataFrame(results)
    csv_path = os.path.join(args.outdir, "screening_results.csv")
    df.to_csv(csv_path, index=False)

    if scored_geoms:
        out_gdf = gpd.GeoDataFrame(pd.concat(scored_geoms, ignore_index=True), crs=cfg.target_crs)
        out_gdf = out_gdf.to_crs("EPSG:4326")   # GeoJSON standard / web maps use WGS84
        out_gdf.to_file(os.path.join(args.outdir, "screening_results.geojson"), driver="GeoJSON")

        # Self-contained map: simplify geometry and embed it directly in the HTML, so the
        # map simply opens — no separate file to load, no browser size/parse limits, no
        # drag-and-drop. This is the file to hand off.
        try:
            map_gdf = out_gdf.copy()
            map_gdf["geometry"] = map_gdf.geometry.simplify(cfg.map_simplify_tol,
                                                            preserve_topology=True)
            html = build_map_html(map_gdf.to_json(), len(map_gdf))
            with open(os.path.join(args.outdir, "playa_map.html"), "w", encoding="utf-8") as fh:
                fh.write(html)
        except Exception as e:
            warnings.warn(f"map generation failed: {e}")

    with open(os.path.join(args.outdir, "run_config.json"), "w") as fh:
        json.dump({"scenario": cfg.scenario, "crs": cfg.target_crs,
                   "buffer_m": cfg.buffer_m, "weights": cfg.weights}, fh, indent=2)

    print(f"\nDone. {len(results)} site(s) screened. Outputs in: {args.outdir}")


_MAP_TEMPLATE = r'''<!DOCTYPE html>
<html lang="en"><head><meta charset="utf-8">
<meta name="viewport" content="width=device-width, initial-scale=1">
<title>Playa Constraint Screening &mdash; Map</title>
<link rel="stylesheet" href="https://unpkg.com/leaflet@1.9.4/dist/leaflet.css">
<style>
 :root{--paper:#f3eee3;--ink:#20241f;--basin:#1d3b2a;--clay:#b4612f;--muted:#6f6657;--line:#d9d0be;--panel:#fbf8f1;}
 html,body{height:100%;margin:0;font-family:-apple-system,Segoe UI,Roboto,Helvetica,Arial,sans-serif;color:var(--ink);}
 #map{position:absolute;inset:0;background:var(--paper);}
 .panel{position:absolute;top:16px;left:16px;z-index:1000;width:300px;max-width:calc(100% - 32px);
  background:var(--panel);border:1px solid var(--line);border-radius:10px;
  box-shadow:0 14px 40px -18px rgba(29,59,42,.55);padding:14px 16px;}
 .kicker{font-size:10px;letter-spacing:.2em;text-transform:uppercase;color:var(--clay);}
 h1{font-size:19px;margin:4px 0 3px;color:var(--basin);}
 .sub{font-size:11.5px;color:var(--muted);}
 .count{font-size:12px;color:var(--muted);margin:10px 0;}
 .count b{color:var(--basin);font-size:14px;}
 .lt{font-size:10px;letter-spacing:.16em;text-transform:uppercase;color:var(--muted);margin-bottom:6px;}
 .lg{display:flex;align-items:center;gap:8px;font-size:12px;margin-bottom:5px;}
 .sw{width:15px;height:15px;border-radius:4px;flex:none;border:1px solid rgba(0,0,0,.15);}
 .pp-name{font-size:15px;color:var(--basin);font-weight:700;}
 .pp-score{display:flex;align-items:baseline;gap:6px;margin:6px 0 8px;}
 .pp-num{font-size:26px;font-weight:700;}
 .pp-badge{color:#fff;font-size:10px;font-weight:700;padding:2px 7px;border-radius:4px;}
 table.dims{width:100%;border-collapse:collapse;font-size:11.5px;}
 table.dims td{padding:2px 0;border-bottom:1px solid var(--line);}
 table.dims td:last-child{text-align:right;font-weight:700;}
 .na{color:#aaa;font-style:italic;font-weight:400;}
 .ppr{display:inline-block;margin-top:9px;font-size:11px;font-weight:700;text-decoration:none;
  color:var(--paper);background:var(--basin);padding:6px 11px;border-radius:6px;}
 .msg{position:absolute;top:50%;left:50%;transform:translate(-50%,-50%);z-index:1000;background:#fff;
  border:1px solid var(--line);border-radius:10px;padding:18px 22px;font-size:13px;max-width:340px;display:none;}
</style></head><body>
<div id="map"></div>
<div class="panel">
 <div class="kicker">Environmental Permitting &middot; GIS</div>
 <h1>Playa Constraint Screening</h1>
 <div class="sub">Per-site playa risk &middot; Southern High Plains</div>
 <div class="count"><b>__NFEATURES__</b> sites</div>
 <div class="lt">Constraint rating</div>
 <div class="lg"><span class="sw" style="background:#2e7d32"></span>Low (0&ndash;25)</div>
 <div class="lg"><span class="sw" style="background:#e0a106"></span>Moderate (25&ndash;50)</div>
 <div class="lg"><span class="sw" style="background:#ef6c00"></span>High (50&ndash;75)</div>
 <div class="lg"><span class="sw" style="background:#c62828"></span>Severe (75&ndash;100)</div>
</div>
<div id="msg" class="msg"></div>
<script src="https://unpkg.com/leaflet@1.9.4/dist/leaflet.js"></script>
<script>
var DATA = __GEOJSON__;
var RATING = {Low:'#2e7d32',Moderate:'#e0a106',High:'#ef6c00',Severe:'#c62828'};
var DIMS = [['dim_esg','1. ESG / Reputation'],['dim_hydrology','2. Hydrology & Drainage'],
 ['dim_migratory','3. Migratory Bird Exposure'],['dim_cwa404','4. CWA Section 404'],
 ['dim_pljv','5. PLJV Siting Guidance']];
function showMsg(t){var m=document.getElementById('msg');m.innerHTML=t;m.style.display='block';}
if(typeof L==='undefined'){
 showMsg('The map library (Leaflet) could not load. Connect to the internet and reopen this file.');
}else{
 var imagery=L.tileLayer('https://server.arcgisonline.com/ArcGIS/rest/services/World_Imagery/MapServer/tile/{z}/{y}/{x}',{attribution:'Esri',maxZoom:19});
 var light=L.tileLayer('https://{s}.basemaps.cartocdn.com/light_all/{z}/{x}/{y}{r}.png',{attribution:'OSM, CARTO',maxZoom:20});
 var map=L.map('map',{layers:[imagery],preferCanvas:true}).setView([33.9,-101.0],7);
 L.control.layers({'Satellite':imagery,'Light map':light},{},{position:'topright'}).addTo(map);
 function ratingOf(p){if(p.rating&&RATING[p.rating])return p.rating;var s=+p.composite_score||0;
  return s>=75?'Severe':s>=50?'High':s>=25?'Moderate':'Low';}
 function popup(p){var r=ratingOf(p),c=RATING[r];
  var rows=DIMS.map(function(d){var v=p[d[0]];var cell=(v==null||v==='')?'<td class="na">n/a</td>':'<td>'+(+v).toFixed(0)+'</td>';return '<tr><td>'+d[1]+'</td>'+cell+'</tr>';}).join('');
  var fl=p.flags?'<div style="margin-top:7px;font-size:10.5px;color:#b4612f"><b>Flags:</b> '+p.flags+'</div>':'';
  var rf=p.report_file?'<a class="ppr" href="'+p.report_file+'" target="_blank" rel="noopener">Open full report &#8599;</a>':'';
  return '<div class="pp-name">'+(p.site_name||'Site')+'</div><div class="pp-score"><span class="pp-num" style="color:'+c+'">'+(+p.composite_score||0).toFixed(1)+'</span><span class="pp-badge" style="background:'+c+'">'+r+'</span><span style="font-size:11px;color:#888">/100</span></div><table class="dims">'+rows+'</table>'+fl+rf;}
 try{
  var layer=L.geoJSON(DATA,{style:function(f){var c=RATING[ratingOf(f.properties)];return {color:c,weight:1,fillColor:c,fillOpacity:.45};},
   onEachFeature:function(f,l){l.bindPopup(popup(f.properties||{}),{maxWidth:300});}}).addTo(map);
  map.fitBounds(layer.getBounds().pad(0.1));
 }catch(e){showMsg('Could not draw the data: '+e.message);}
}
</script></body></html>'''


def build_map_html(geojson_str, n_features):
    """Return a self-contained map HTML with the GeoJSON embedded (no external data file)."""
    return _MAP_TEMPLATE.replace("__NFEATURES__", str(n_features)).replace("__GEOJSON__", geojson_str)


def build_parser():
    p = argparse.ArgumentParser(description="Playa-aware environmental constraint screening.")
    p.add_argument("--boundary", required=True, help="Project/candidate boundary vector file")
    p.add_argument("--boundary-layer", dest="boundary_layer", default=None,
                   help="Layer name inside the boundary file (for multi-layer GeoPackages)")
    p.add_argument("--playas", help="PLJV probable playas")
    p.add_argument("--clusters", help="PLJV playa clusters")
    p.add_argument("--nwi", help="NWI wetlands")
    p.add_argument("--floodplain", help="FEMA NFHL floodplain")
    p.add_argument("--critical-habitat", dest="critical_habitat", help="USFWS critical habitat")
    p.add_argument("--nrhp", help="National Register of Historic Places (points)")
    p.add_argument("--id-field", help="Boundary attribute to use as the site name")
    p.add_argument("--crs", default="EPSG:5070", help="Equal-area analysis CRS (default EPSG:5070)")
    p.add_argument("--scenario", choices=["conforming", "proposed"], default="conforming",
                   help="WOTUS scenario for the §404 dimension")
    p.add_argument("--buffer", type=float, default=100.0, help="Playa avoidance buffer (meters)")
    p.add_argument("--cluster-field", dest="cluster_field", default=None,
                   help="Column in the playa layer to derive cluster overlap (e.g. 'cluster')")
    p.add_argument("--cluster-min", dest="cluster_min", type=float, default=0.0,
                   help="A playa counts as clustered when its cluster-field value exceeds this")
    p.add_argument("--cluster-pct-cap", dest="cluster_pct_cap", type=float, default=5.0,
                   help="%% of site covered by clustered playas that maxes cluster intensity "
                        "(lower = more sensitive; default 5.0)")
    p.add_argument("--outdir", default="outputs", help="Output directory")
    p.add_argument("--template", default="report_template.html", help="HTML report template")
    return p


if __name__ == "__main__":
    run(build_parser().parse_args())
