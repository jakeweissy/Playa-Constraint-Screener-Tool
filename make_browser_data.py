#!/usr/bin/env python3
"""
make_browser_data.py  —  RUN ONCE to build the data file the browser screener uses.

It reads your playa shapefile, clips it to Texas (if you pass a Texas boundary),
keeps only the two fields the app needs (cluster flag + area), simplifies the shapes,
reprojects to WGS84, and writes `playa_app_data.js`. The browser app loads that file
directly (via a <script> tag), so no web server is needed afterward.

Usage (from your project root, venv active):
    python make_browser_data.py --playas data/probable_playas_v5_shapefiles.shp --texas data/tx_counties_clean.gpkg

Then put `playa_app_data.js` in the SAME folder as `playa_screener.html` and open the html.
"""
import argparse
import os
import geopandas as gpd
import pandas as pd


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--playas", required=True, help="PLJV probable playas shapefile")
    ap.add_argument("--texas", help="Texas boundary (e.g. tx_counties_clean.gpkg) to clip + outline")
    ap.add_argument("--texas-layer", default=None, help="layer name if --texas is multi-layer")
    ap.add_argument("--cluster-field", default="cluster")
    ap.add_argument("--simplify", type=float, default=0.0015, help="degrees (~150 m)")
    ap.add_argument("--out", default="playa_app_data.js")
    a = ap.parse_args()

    print("Reading playas ...")
    playas = gpd.read_file(a.playas).to_crs("EPSG:4326")

    tx_outline_json = None
    if a.texas:
        print("Reading Texas boundary ...")
        tx = (gpd.read_file(a.texas, layer=a.texas_layer) if a.texas_layer
              else gpd.read_file(a.texas)).to_crs("EPSG:4326")
        tx_union = tx.geometry.union_all() if hasattr(tx.geometry, "union_all") else tx.geometry.unary_union
        print("Clipping playas to Texas ...")
        playas = gpd.clip(playas, tx_union).reset_index(drop=True)
        tx_outline_json = gpd.GeoSeries([tx_union], crs="EPSG:4326") \
            .simplify(0.01, preserve_topology=True).to_json()

    print("Computing areas + cluster flags ...")
    area_ha = (playas.to_crs("EPSG:5070").geometry.area / 10000.0).round(3)
    if a.cluster_field in playas.columns:
        cflag = (pd.to_numeric(playas[a.cluster_field], errors="coerce").fillna(0) > 0).astype(int)
    else:
        print(f"  (cluster field '{a.cluster_field}' not found; marking all 0)")
        cflag = pd.Series([0] * len(playas))

    slim = gpd.GeoDataFrame(
        {"c": cflag.values, "a": area_ha.values},
        geometry=playas.geometry.simplify(a.simplify, preserve_topology=True),
        crs="EPSG:4326",
    )
    slim = slim[slim.geometry.notna() & ~slim.geometry.is_empty].reset_index(drop=True)

    print(f"Writing {len(slim)} playas to {a.out} ...")
    with open(a.out, "w", encoding="utf-8") as f:
        f.write("window.PLAYA_DATA=" + slim.to_json() + ";\n")
        if tx_outline_json:
            f.write("window.TEXAS_BOUNDARY=" + tx_outline_json + ";\n")

    mb = os.path.getsize(a.out) / 1e6
    print(f"Done. {a.out} is {mb:.1f} MB. Put it next to playa_screener.html and open the html.")
    if mb > 25:
        print("  (Large file — if the app is slow to open, re-run with a bigger --simplify, e.g. 0.003)")


if __name__ == "__main__":
    main()
