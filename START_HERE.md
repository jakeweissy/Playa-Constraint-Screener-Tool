# START HERE — Playa Constraint Screening (clean setup)

Follow this top to bottom once and you'll have a folder you can hand to your boss:
an interactive statewide map, a report for every county, and a ranked spreadsheet.

**What's different now:** the map comes out *self-contained* — the data is built right
into the file. You just double-click it. No web server, no dragging files, no internet
needed except to show the background satellite imagery.

---

## 1. The folder layout

Make one folder (anywhere fast and local — **not** inside OneDrive), e.g.
`C:\GIS\Playa Screening`. It should contain exactly this:

```
Playa Screening\
├── playa_screen.py          ← the tool
├── report_template.html     ← used by the tool to build reports
├── check_setup.py           ← preflight checker
├── requirements.txt
├── START_HERE.md            ← this file
└── data\
    ├── probable_playas_v5_shapefiles.shp   (+ .dbf .prj .shx .cpg — keep them together)
    ├── tx_counties_clean.gpkg              ← Texas counties (single layer, ~254 features)
    ├── chab_tx.gpkg                        ← critical habitat, clipped to Texas
    └── nrhp_tx.gpkg                        ← NRHP points, clipped to Texas
```

**Delete these leftover files if you have them** — they were debugging aids and aren't
needed anymore: `fix_geojson.py`, `simplify_geojson.py`, any `playa_map.html` sitting in
the project root, `demo_boundary.geojson`, and any `screening_results_clean*.geojson` /
`*_simplified.geojson` inside output folders. The tool now handles all of that internally.

(`nwi_tx.gpkg` is optional — see the note in step 4.)

---

## 2. Set up Python

Open the folder in VS Code (**File ▸ Open Folder**), then **Terminal ▸ New Terminal**, and:

```
python -m venv .venv
.venv\Scripts\Activate.ps1
pip install -r requirements.txt
```

You should see `(.venv)` at the start of the prompt and the install finish cleanly.
(If activation is blocked, run `Set-ExecutionPolicy -Scope Process RemoteSigned` once, then
re-activate.)

---

## 3. Preflight check

```
python check_setup.py
```

Everything required should say `[ OK ]`. Fix any `[MISS]`. Optional layers you don't have
show `[WARN]`, which is fine. If a GeoPackage shows more than one layer, it lists them —
use the one with ~254 features.

---

## 4. Run the scan

Copy this as **one line** (PowerShell doesn't like the `\` breaks):

```
python playa_screen.py --boundary data/tx_counties_clean.gpkg --playas data/probable_playas_v5_shapefiles.shp --cluster-field cluster --critical-habitat data/chab_tx.gpkg --nrhp data/nrhp_tx.gpkg --id-field NAME --scenario conforming --outdir outputs_tx
```

It prints each county and its score as it goes, and finishes in a few minutes. When done,
you'll have an `outputs_tx\` folder containing everything.

> **About NWI:** it's left out above because it's a very large layer and slow to process,
> and it only feeds part of one dimension (CWA §404). If you want it included, add
> `--nwi data/nwi_tx.gpkg` to the command and expect a much longer run. Your scores are
> valid either way. FEMA floodplain is intentionally omitted and shows as "not assessed."

---

## 5. Look at your results

Open the `outputs_tx` folder and **double-click `playa_map.html`.** It opens in your
browser and shows all 254 counties shaded by risk, with no extra steps. Click any county
for its score breakdown, then click **"Open full report ↗"** to open that county's report.
To save a report as a file, in the report tab press **Ctrl+P ▸ Save as PDF**.

Also in that folder:
- **`screening_results.csv`** — open in Excel and sort by `composite_score` to rank the
  riskiest counties.
- **`report_<County>.html`** — the 254 individual reports (also reachable from the map).
- `screening_results.geojson`, `run_config.json` — the raw data and the settings used.

> The map needs internet only to draw the satellite/▸map background. If you're offline the
> county shapes still appear on a blank background, and if the map library itself can't
> load you'll see a short message saying so instead of a blank screen.

---

## 6. You’re Done!

The deliverable is the **`outputs_tx`** folder. Right-click it ▸ **Send to ▸ Compressed
(zipped) folder**, and share the zip. Your boss double-clicks `playa_map.html` inside it —
nothing to install. Include `screening_results.csv` for the ranked table.

If you also want them to see the methodology, include `playa_screen.py` and this guide.

---

## Doing specific project sites instead of counties

Same command — point `--boundary` at a file of your project polygons and set `--id-field`
to the column that names them, e.g.:

```
python playa_screen.py --boundary data/my_sites.gpkg --playas data/probable_playas_v5_shapefiles.shp --cluster-field cluster --critical-habitat data/chab_tx.gpkg --nrhp data/nrhp_tx.gpkg --id-field site_name --scenario conforming --outdir outputs_sites
```

---

## If something goes wrong

- **"file not found"** → check the path is spelled exactly and the file is in `data\`.
  Run `python check_setup.py` to see what it finds.
- **Map opens but is blank / shows a message** → you're likely offline; reconnect and
  reopen. The shapes need the map library (loaded from the internet the first time).
- **A run errors partway** → copy the red error text and keep it; it names the problem
  (usually a path or a layer issue).

