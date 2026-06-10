#!/usr/bin/env python3
"""Generate a self-contained HTML monitor page for a VASP NEB directory.

Usage:
    python3 neb_html_monitor.py /path/to/NEB_working_directory
    python3 neb_html_monitor.py /path/to/NEB_working_directory -o neb_report.html
"""

from __future__ import annotations

import argparse
from datetime import datetime
import html
import json
import math
import os
import re
from pathlib import Path
from typing import Any


FLOAT_RE = r"[-+]?(?:\d+\.\d*|\.\d+|\d+)(?:[Ee][-+]?\d+)?"


def read_text(path: Path, limit_bytes: int | None = None) -> str:
    if not path.exists() or not path.is_file():
        return ""
    if limit_bytes is None:
        return path.read_text(errors="replace")
    with path.open("rb") as handle:
        return handle.read(limit_bytes).decode(errors="replace")


def parse_incar(text: str) -> list[dict[str, str]]:
    settings: list[dict[str, str]] = []
    for raw in text.splitlines():
        line = raw.strip()
        if not line or line.startswith("#") or line.startswith("!"):
            continue
        match = re.match(r"^([A-Za-z][A-Za-z0-9_]*)\s*=\s*(.*)$", line)
        if not match:
            continue
        key, value = match.groups()
        comment = ""
        for token in ("!", "#", "("):
            if token in value:
                value, comment = value.split(token, 1)
                comment = comment.rstrip(")")
                break
        settings.append(
            {
                "key": key.strip(),
                "value": value.strip(),
                "comment": comment.strip(),
                "raw": raw.rstrip(),
            }
        )
    return settings


def numeric_image_dirs(workdir: Path) -> list[Path]:
    dirs = [p for p in workdir.iterdir() if p.is_dir() and re.fullmatch(r"\d+", p.name)]
    return sorted(dirs, key=lambda p: int(p.name))


def parse_oszicar(path: Path) -> list[dict[str, float]]:
    text = read_text(path)
    rows: list[dict[str, float]] = []
    pattern = re.compile(
        rf"^\s*(\d+)\s+F=\s*({FLOAT_RE})\s+E0=\s*({FLOAT_RE})\s+d\s*E\s*=\s*({FLOAT_RE})",
        re.MULTILINE,
    )
    for step, free_energy, e0, delta_e in pattern.findall(text):
        rows.append(
            {
                "step": int(step),
                "free_energy": float(free_energy),
                "e0": float(e0),
                "delta_e": float(delta_e),
            }
        )
    return rows


def parse_outcar_toten(path: Path) -> float | None:
    if not path.exists():
        return None
    last: float | None = None
    pattern = re.compile(rf"free\s+energy\s+TOTEN\s+=\s*({FLOAT_RE})")
    with path.open(errors="replace") as handle:
        for line in handle:
            match = pattern.search(line)
            if match:
                last = float(match.group(1))
    return last


def parse_force_conv(path: Path) -> list[dict[str, float]]:
    rows: list[dict[str, float]] = []
    for line in read_text(path).splitlines():
        parts = line.split()
        if len(parts) >= 2:
            try:
                rows.append({"step": int(float(parts[0])), "force": float(parts[1])})
            except ValueError:
                pass
    return rows


def parse_outcar_forces(path: Path) -> list[dict[str, float]]:
    rows: list[dict[str, float]] = []
    if not path.exists():
        return rows
    pattern = re.compile(rf"FORCES:\s+max atom,\s+RMS\s+({FLOAT_RE})\s+({FLOAT_RE})")
    with path.open(errors="replace") as handle:
        for line in handle:
            match = pattern.search(line)
            if match:
                rows.append(
                    {
                        "step": len(rows) + 1,
                        "max_force": float(match.group(1)),
                        "rms_force": float(match.group(2)),
                    }
                )
    return rows


def parse_poscar(path: Path) -> dict[str, Any] | None:
    lines = [line.strip() for line in read_text(path).splitlines() if line.strip()]
    if len(lines) < 8:
        return None
    try:
        scale = float(lines[1].split()[0])
        lattice = [[float(x) * scale for x in lines[i].split()[:3]] for i in range(2, 5)]
    except (ValueError, IndexError):
        return None

    names_line = lines[5].split()
    count_line_index = 6
    try:
        counts = [int(x) for x in lines[count_line_index].split()]
    except ValueError:
        names_line = [f"X{i + 1}" for i in range(len(lines[5].split()))]
        count_line_index = 5
        try:
            counts = [int(x) for x in lines[count_line_index].split()]
        except ValueError:
            return None

    coord_mode_index = count_line_index + 1
    if lines[coord_mode_index].lower().startswith("s"):
        coord_mode_index += 1
    direct = lines[coord_mode_index].lower().startswith(("d", "fractional"))
    coord_start = coord_mode_index + 1
    elements = []
    for element, count in zip(names_line, counts):
        elements.extend([element] * count)
    natoms = sum(counts)
    atoms = []
    for element, line in zip(elements, lines[coord_start : coord_start + natoms]):
        parts = line.split()
        if len(parts) < 3:
            continue
        try:
            a, b, c = [float(x) for x in parts[:3]]
        except ValueError:
            continue
        if direct:
            x = a * lattice[0][0] + b * lattice[1][0] + c * lattice[2][0]
            y = a * lattice[0][1] + b * lattice[1][1] + c * lattice[2][1]
            z = a * lattice[0][2] + b * lattice[1][2] + c * lattice[2][2]
        else:
            x, y, z = a * scale, b * scale, c * scale
        atoms.append({"e": element, "x": x, "y": y, "z": z})
    if not atoms:
        return None
    return {"source": path.name, "atoms": atoms}


def build_movie_frames(image_dirs: list[Path]) -> list[dict[str, Any]]:
    frames = []
    for image_dir in image_dirs:
        structure = parse_poscar(image_dir / "CONTCAR") or parse_poscar(image_dir / "POSCAR")
        if structure:
            structure["image"] = image_dir.name
            structure["comment"] = f"Image {image_dir.name} from {structure['source']}"
            frames.append(structure)
    return frames


def build_image_data(image_dirs: list[Path]) -> list[dict[str, Any]]:
    images = []
    for image_dir in image_dirs:
        oszicar_rows = parse_oszicar(image_dir / "OSZICAR")
        force_rows = parse_outcar_forces(image_dir / "OUTCAR")
        if not force_rows:
            force_rows = parse_force_conv(image_dir / "force.conv")
        final_energy = oszicar_rows[-1]["free_energy"] if oszicar_rows else parse_outcar_toten(image_dir / "OUTCAR")
        images.append(
            {
                "image": image_dir.name,
                "path": str(image_dir),
                "ionic_steps": len(oszicar_rows),
                "final_energy": final_energy,
                "energy_series": oszicar_rows,
                "force_series": force_rows,
                "files": sorted(p.name for p in image_dir.iterdir() if p.is_file()),
            }
        )
    return images


def make_relative_profile(images: list[dict[str, Any]]) -> list[dict[str, Any]]:
    energies = [img["final_energy"] for img in images if img["final_energy"] is not None]
    if not energies:
        return []
    reference = energies[0]
    return [
        {
            "image": img["image"],
            "energy": img["final_energy"],
            "relative_energy": None if img["final_energy"] is None else img["final_energy"] - reference,
        }
        for img in images
    ]


def status_from_data(settings: list[dict[str, str]], images: list[dict[str, Any]]) -> dict[str, Any]:
    lookup = {row["key"].upper(): row["value"] for row in settings}
    ediffg = None
    ediff = None
    try:
        ediffg = abs(float(lookup.get("EDIFFG", "")))
    except ValueError:
        pass
    try:
        ediff = abs(float(lookup.get("EDIFF", "")))
    except ValueError:
        pass
    latest_forces = []
    moving_images = images[1:-1] if len(images) > 2 else images
    for img in moving_images:
        forces = img["force_series"]
        if forces:
            row = forces[-1]
            latest_forces.append(row.get("force", row.get("max_force")))
    max_force = max([f for f in latest_forces if f is not None], default=None)
    converged = bool(ediffg and max_force is not None and max_force <= ediffg)
    return {"ediff_abs": ediff, "ediffg_abs": ediffg, "latest_max_force": max_force, "force_converged": converged}


def summarize_potcar(path: Path) -> str:
    text = read_text(path, 2_000_000)
    if not text:
        return "POTCAR not found."
    keep_prefixes = ("TITEL", "VRHFIN", "POMASS", "ENMAX", "LEXCH")
    lines = []
    current_block: list[str] = []
    for raw in text.splitlines():
        stripped = raw.strip()
        if any(stripped.startswith(prefix) for prefix in keep_prefixes):
            current_block.append(stripped)
        if stripped.startswith("End of Dataset"):
            if current_block:
                lines.append("\n".join(current_block))
                current_block = []
    if current_block:
        lines.append("\n".join(current_block))
    if not lines:
        return "POTCAR found, but no standard TITEL/VRHFIN/POMASS/ENMAX metadata lines were detected."
    return "\n\n".join(f"Dataset {i + 1}\n{block}" for i, block in enumerate(lines))


def html_page(data: dict[str, Any]) -> str:
    payload = json.dumps(data, ensure_ascii=False).replace("</", "<\\/")
    raw_incar = html.escape(data["raw_files"].get("INCAR", ""))
    raw_kpoints = html.escape(data["raw_files"].get("KPOINTS", ""))
    raw_slurm = html.escape(data["raw_files"].get("slurm.sh", ""))
    raw_potcar = html.escape(data["raw_files"].get("POTCAR", ""))
    title = html.escape(data["title"])
    return f"""<!doctype html>
<html lang="en">
<head>
  <meta charset="utf-8">
  <meta name="viewport" content="width=device-width, initial-scale=1">
  <title>{title}</title>
  <style>
    :root {{ color-scheme: light; --ink:#17202a; --muted:#65717f; --line:#d8dee8; --panel:#f7f9fc; --accent:#1b7f83; --warn:#b54708; }}
    body {{ margin:0; font:14px/1.45 system-ui, -apple-system, Segoe UI, sans-serif; color:var(--ink); background:#fff; }}
    header {{ padding:22px 28px; border-bottom:1px solid var(--line); background:#f9fbfd; }}
    h1 {{ margin:0 0 6px; font-size:24px; letter-spacing:0; }}
    h2 {{ margin:0 0 14px; font-size:17px; }}
    main {{ padding:22px 28px 40px; display:grid; gap:18px; }}
    section {{ border:1px solid var(--line); border-radius:8px; padding:16px; background:#fff; }}
    .cards {{ display:grid; grid-template-columns:repeat(auto-fit,minmax(170px,1fr)); gap:12px; }}
    .card {{ border:1px solid var(--line); border-radius:8px; padding:12px; background:var(--panel); }}
    .label {{ color:var(--muted); font-size:12px; text-transform:uppercase; }}
    .value {{ font-size:22px; font-weight:700; margin-top:4px; }}
    .grid {{ display:grid; grid-template-columns:1fr 1fr; gap:16px; align-items:start; }}
    .file-grid {{ display:grid; grid-template-columns:1fr 1fr; gap:16px; }}
    .file-panel pre {{ height:260px; max-height:none; margin:0; }}
    .file-panel h2 {{ margin-top:0; }}
    .chart-heading {{ display:flex; align-items:center; justify-content:space-between; gap:12px; margin-bottom:14px; }}
    .chart-heading h2 {{ margin:0; }}
    .chart-control {{ display:flex; align-items:center; gap:8px; color:var(--muted); font-size:13px; }}
    .chart-control select {{ border:1px solid var(--line); border-radius:6px; background:#fff; color:var(--ink); padding:5px 8px; font:inherit; }}
    .movie-grid {{ display:grid; grid-template-columns:1fr 1fr; gap:14px; }}
    .view-title {{ margin:0 0 8px; color:var(--muted); font-size:13px; font-weight:650; }}
    table {{ border-collapse:collapse; width:100%; font-size:13px; }}
    th, td {{ border-bottom:1px solid var(--line); padding:7px 8px; text-align:left; vertical-align:top; }}
    th {{ background:var(--panel); font-weight:650; }}
    svg {{ width:100%; height:300px; display:block; border:1px solid var(--line); border-radius:6px; background:#fff; }}
    svg.convergence-chart {{ height:380px; }}
    svg text {{ font-family:inherit; }}
    canvas {{ width:100%; aspect-ratio:1.65 / 1; height:auto; display:block; border:1px solid var(--line); border-radius:6px; background:#0f1720; }}
    pre {{ overflow:auto; white-space:pre-wrap; background:#111827; color:#e5e7eb; border-radius:6px; padding:12px; max-height:360px; }}
    .controls {{ display:flex; gap:10px; align-items:center; margin-top:10px; flex-wrap:wrap; }}
    button {{ border:1px solid var(--line); background:#fff; border-radius:6px; padding:7px 11px; cursor:pointer; }}
    input[type=range] {{ flex:1; min-width:200px; }}
    .muted {{ color:var(--muted); }}
    .ok {{ color:var(--accent); }}
    .warn {{ color:var(--warn); }}
    @media (max-width: 860px) {{ .grid, .movie-grid, .file-grid {{ grid-template-columns:1fr; }} header, main {{ padding-left:16px; padding-right:16px; }} }}
  </style>
</head>
<body>
  <header>
    <h1>{title}</h1>
    <div class="muted">Working directory: <code>{html.escape(data["workdir"])}</code></div>
    <div class="muted">Time: {html.escape(data["generated_at"])}</div>
  </header>
  <main>
    <section>
      <div class="cards" id="cards"></div>
    </section>
    <section>
      <h2>NEB Movie</h2>
      <div class="movie-grid">
        <div>
          <div class="view-title">Front view (x-y)</div>
          <canvas id="frontCanvas" width="1200" height="720"></canvas>
        </div>
        <div>
          <div class="view-title">Side view (x-z)</div>
          <canvas id="sideCanvas" width="1200" height="720"></canvas>
        </div>
      </div>
      <div class="controls">
        <button id="playBtn">Play</button>
        <input id="frameSlider" type="range" min="0" value="0">
        <span id="frameLabel" class="muted"></span>
      </div>
    </section>
    <section>
      <h2>NEB Energy Profile</h2>
      <svg id="profileChart"></svg>
    </section>
    <section class="grid">
      <div>
        <div class="chart-heading">
          <h2>Energy Convergence</h2>
          <label class="chart-control">Show
            <select id="energyStepWindow">
              <option value="all">all steps</option>
              <option value="50">last 50</option>
              <option value="30" selected>last 30</option>
              <option value="10">last 10</option>
            </select>
          </label>
        </div>
        <svg id="energyChart" class="convergence-chart"></svg>
      </div>
      <div>
        <div class="chart-heading">
          <h2>Force Convergence</h2>
          <label class="chart-control">Show
            <select id="forceStepWindow">
              <option value="all">all steps</option>
              <option value="50">last 50</option>
              <option value="30" selected>last 30</option>
              <option value="10">last 10</option>
            </select>
          </label>
        </div>
        <svg id="forceChart" class="convergence-chart"></svg>
      </div>
    </section>
    <section>
      <h2>Image Summary</h2>
      <div id="imageTable"></div>
    </section>
    <section class="file-grid">
      <div class="file-panel"><h2>INCAR</h2><pre>{raw_incar}</pre></div>
      <div class="file-panel"><h2>KPOINTS</h2><pre>{raw_kpoints}</pre></div>
      <div class="file-panel"><h2>POTCAR</h2><pre>{raw_potcar}</pre></div>
      <div class="file-panel"><h2>Slurm Script</h2><pre>{raw_slurm}</pre></div>
    </section>
  </main>
  <script id="neb-data" type="application/json">{payload}</script>
  <script>
    const DATA = JSON.parse(document.getElementById('neb-data').textContent);
    const fmt = (x, d=4) => x === null || x === undefined || Number.isNaN(x) ? 'n/a' : Number(x).toFixed(d);
    const cards = [
      ['Images', DATA.images.length],
      ['Ionic steps', Math.max(...DATA.images.map(x => x.ionic_steps || 0), 0)],
      ['Barrier', fmt(Math.max(...DATA.profile.map(x => x.relative_energy ?? 0)), 4) + ' eV'],
      ['Latest max force', fmt(DATA.status.latest_max_force, 4) + ' eV/A'],
      ['Force target', fmt(DATA.status.ediffg_abs, 4) + ' eV/A'],
      ['Status', DATA.status.force_converged ? 'Converged' : 'Running / not converged']
    ];
    document.getElementById('cards').innerHTML = cards.map(([k,v]) => `<div class="card"><div class="label">${{k}}</div><div class="value ${{String(v).startsWith('Running') ? 'warn' : ''}}">${{v}}</div></div>`).join('');

    function table(headers, rows) {{
      return `<table><thead><tr>${{headers.map(h=>`<th>${{h}}</th>`).join('')}}</tr></thead><tbody>${{rows.map(r=>`<tr>${{r.map(c=>`<td>${{c}}</td>`).join('')}}</tr>`).join('')}}</tbody></table>`;
    }}

    document.getElementById('imageTable').innerHTML = table(
      ['Image', 'Ionic step', 'Final energy (eV)', 'Relative energy vs 00 (eV)', 'Final max force (eV/A)'],
      DATA.images.map((img, i) => [
        img.image,
        img.ionic_steps,
        fmt(img.final_energy, 6),
        fmt(DATA.profile[i]?.relative_energy, 6),
        fmt(img.force_series.at(-1)?.force ?? img.force_series.at(-1)?.max_force, 6)
      ])
    );

    function drawChart(svgId, series, opts={{}}) {{
      const svg = document.getElementById(svgId);
      const W = 900, H = opts.height || 340, m = {{l:96, r:28, t:28, b:92}};
      svg.setAttribute('viewBox', `0 0 ${{W}} ${{H}}`);
      svg.innerHTML = '';
      const all = series.flatMap(s => s.points).filter(p => Number.isFinite(p.x) && Number.isFinite(p.y));
      if (!all.length) {{ svg.innerHTML = '<text x="20" y="36" fill="#65717f">No data found</text>'; return; }}
      const xmin = Math.min(...all.map(p=>p.x)), xmax = Math.max(...all.map(p=>p.x));
      const thresholdValues = (opts.thresholds || []).map(t => t.value).filter(Number.isFinite);
      const ymin0 = Math.min(...all.map(p=>p.y), ...thresholdValues), ymax0 = Math.max(...all.map(p=>p.y), ...thresholdValues);
      function niceTicks(min, max, targetCount=5) {{
        const span = Math.max(max - min, 1e-12);
        const rawStep = span / Math.max(targetCount - 1, 1);
        const power = Math.pow(10, Math.floor(Math.log10(rawStep)));
        const error = rawStep / power;
        const nice = error >= 7.5 ? 10 : error >= 3.5 ? 5 : error >= 1.5 ? 2 : 1;
        const step = nice * power;
        const niceMin = Math.floor(min / step) * step;
        const niceMax = Math.ceil(max / step) * step;
        const ticks = [];
        for (let v = niceMin; v <= niceMax + step * 0.5; v += step) ticks.push(Number(v.toPrecision(12)));
        return {{ticks, min:niceMin, max:niceMax}};
      }}
      const tickInfo = niceTicks(ymin0, ymax0, opts.yTickCount || 5);
      const ymin = tickInfo.min, ymax = tickInfo.max;
      const sx = x => m.l + (x - xmin) / (xmax - xmin || 1) * (W - m.l - m.r);
      const sy = y => H - m.b - (y - ymin) / (ymax - ymin || 1) * (H - m.t - m.b);
      const line = (x1,y1,x2,y2,stroke='#d8dee8', extra='') => `<line x1="${{x1}}" y1="${{y1}}" x2="${{x2}}" y2="${{y2}}" stroke="${{stroke}}" ${{extra}} />`;
      svg.innerHTML += line(m.l, H-m.b, W-m.r, H-m.b) + line(m.l, m.t, m.l, H-m.b);
      tickInfo.ticks.forEach(y => {{
        const yy = sy(y);
        const decimals = Math.abs(y) >= 100 ? 1 : Math.abs(y) >= 10 ? 1 : 2;
        svg.innerHTML += line(m.l, yy, W-m.r, yy, '#eef2f6') + `<text x="${{m.l-18}}" y="${{yy+5}}" text-anchor="end" font-size="13" fill="#65717f">${{fmt(y, decimals)}}</text>`;
      }});
      const defaultTickCount = opts.xTickCount || Math.min(6, Math.max(2, Math.floor(xmax - xmin) + 1));
      const xTicks = opts.xTicks || Array.from({{length: defaultTickCount}}, (_, i) => {{
        const denom = Math.max(defaultTickCount - 1, 1);
        const value = xmin + i * (xmax - xmin || 1) / denom;
        return {{value, label: String(Math.round(value))}};
      }});
      xTicks.forEach(t => {{
        if (!Number.isFinite(t.value)) return;
        const xx = sx(t.value);
        svg.innerHTML += line(xx, m.t, xx, H-m.b, '#f3f6fa');
        svg.innerHTML += line(xx, H-m.b, xx, H-m.b+7, '#6b7a8a');
        const labelText = String(t.label);
        svg.innerHTML += `<text x="${{xx}}" y="${{H-m.b+25}}" text-anchor="middle" font-size="13" fill="#17202a">${{labelText}}</text>`;
      }});
      const thresholdLabels = [];
      (opts.thresholds || []).forEach(t => {{
        if (!Number.isFinite(t.value)) return;
        const yy = sy(t.value);
        svg.innerHTML += line(m.l, yy, W-m.r, yy, t.color || '#b54708', 'stroke-width="2" stroke-dasharray="6 5"');
        thresholdLabels.push({{
          label:t.label,
          color:t.color || '#b54708',
          lineX:W-m.r-14,
          lineY:yy,
          textX:W-m.r-70,
          textY:H-m.b+58
        }});
      }});
      const colors = ['#1b7f83','#b54708','#5165b5','#7c3f00','#2d6a4f','#7b2cbf','#0f766e','#a21caf'];
      const plotted = [];
      series.forEach((s, idx) => {{
        const pts = s.points.filter(p => Number.isFinite(p.x) && Number.isFinite(p.y));
        if (!pts.length) return;
        plotted.push({{name:s.name, color:colors[idx%colors.length]}});
        const d = pts.map((p,i) => `${{i ? 'L' : 'M'}} ${{sx(p.x)}} ${{sy(p.y)}}`).join(' ');
        svg.innerHTML += `<path d="${{d}}" fill="none" stroke="${{colors[idx%colors.length]}}" stroke-width="2"/>`;
        pts.forEach(p => svg.innerHTML += `<circle cx="${{sx(p.x)}}" cy="${{sy(p.y)}}" r="3" fill="${{colors[idx%colors.length]}}"><title>${{s.name}}: ${{fmt(p.y,5)}}</title></circle>`);
        if (opts.pointLabels) {{
          pts.forEach((p, pointIndex) => {{
            const xx = sx(p.x);
            const yy = sy(p.y);
            const label = opts.pointLabelFormatter ? opts.pointLabelFormatter(p, pointIndex, s) : fmt(p.y, 3);
            const above = pointIndex % 2 === 0;
            const yText = yy + (above ? -12 : 18);
            const labelWidth = Math.max(34, String(label).length * 7 + 8);
            svg.innerHTML += `<rect x="${{xx-labelWidth/2}}" y="${{yText-12}}" width="${{labelWidth}}" height="16" rx="3" fill="white" opacity="0.9"/>`;
            svg.innerHTML += `<text x="${{xx}}" y="${{yText}}" text-anchor="middle" font-size="13" font-weight="650" fill="${{colors[idx%colors.length]}}">${{label}}</text>`;
          }});
        }}
      }});
      thresholdLabels.forEach(t => {{
        const w = Math.max(52, String(t.label).length * 8 + 12);
        const boxX = t.textX - w;
        const boxY = t.textY - 14;
        const arrowStartX = boxX + w * 0.5;
        const arrowStartY = boxY + 18;
        const angle = Math.atan2(t.lineY - arrowStartY, t.lineX - arrowStartX);
        const arrowLeftX = t.lineX - 8 * Math.cos(angle - 0.45);
        const arrowLeftY = t.lineY - 8 * Math.sin(angle - 0.45);
        const arrowRightX = t.lineX - 8 * Math.cos(angle + 0.45);
        const arrowRightY = t.lineY - 8 * Math.sin(angle + 0.45);
        svg.innerHTML += `<line x1="${{arrowStartX}}" y1="${{arrowStartY}}" x2="${{t.lineX}}" y2="${{t.lineY}}" stroke="${{t.color}}" stroke-width="1.8"/>`;
        svg.innerHTML += `<polygon points="${{t.lineX}},${{t.lineY}} ${{arrowLeftX}},${{arrowLeftY}} ${{arrowRightX}},${{arrowRightY}}" fill="${{t.color}}"/>`;
        svg.innerHTML += `<rect x="${{boxX}}" y="${{boxY}}" width="${{w}}" height="18" rx="3" fill="white" opacity="0.94" stroke="#eef2f6"/>`;
        svg.innerHTML += `<text x="${{t.textX-6}}" y="${{t.textY}}" text-anchor="end" font-size="13" font-weight="700" fill="${{t.color}}">${{t.label}}</text>`;
      }});
      svg.innerHTML += `<text x="${{W/2}}" y="${{H-10}}" text-anchor="middle" font-size="13" fill="#65717f">${{opts.xLabel || 'Step / Image'}}</text>`;
      svg.innerHTML += `<text x="24" y="${{H/2}}" transform="rotate(-90 24 ${{H/2}})" text-anchor="middle" font-size="13" fill="#65717f">${{opts.yLabel || ''}}</text>`;
      const legendX = W - m.r - 138;
      const legendY = m.t + 14;
      const legendH = Math.max(24, plotted.length * 19 + 10);
      if (plotted.length) {{
        svg.innerHTML += `<rect x="${{legendX-10}}" y="${{legendY-18}}" width="148" height="${{legendH}}" rx="5" fill="white" opacity="0.82" stroke="#eef2f6"/>`;
      }}
      plotted.forEach((s, idx) => {{
        const y = legendY + idx * 19;
        svg.innerHTML += `<line x1="${{legendX}}" y1="${{y-4}}" x2="${{legendX+22}}" y2="${{y-4}}" stroke="${{s.color}}" stroke-width="3"/>`;
        svg.innerHTML += `<text x="${{legendX+28}}" y="${{y}}" font-size="13" fill="#17202a">${{s.name}}</text>`;
      }});
    }}

    drawChart('profileChart', [{{name:'Relative energy', points: DATA.profile.map((p,i)=>({{x:i, y:p.relative_energy}}))}}], {{xLabel:'Image', yLabel:'Relative energy (eV)', xTicks: DATA.profile.map((p,i)=>({{value:i, label:p.image}})), pointLabels:true, pointLabelFormatter:p => fmt(p.y, 3) + ' eV'}});
    const energyStepWindow = document.getElementById('energyStepWindow');
    function drawEnergyChart() {{
      const selected = energyStepWindow?.value || '30';
      const windowSize = selected === 'all' ? null : Number(selected);
      const labelSuffix = selected === 'all' ? 'all steps' : `last ${{selected}}`;
      drawChart('energyChart', DATA.images.map(img => {{
        const full = img.energy_series.map(r => ({{x:r.step, y:r.e0}}));
        return {{name:'Image ' + img.image, points: windowSize ? full.slice(-windowSize) : full}};
      }}), {{xLabel:`Ionic step (${{labelSuffix}})`, yLabel:'Energy E0 (eV)', xTickCount:5, height:390}});
    }}
    energyStepWindow?.addEventListener('change', drawEnergyChart);
    drawEnergyChart();
    const movingImages = DATA.images.length > 2 ? DATA.images.slice(1, -1) : DATA.images;
    const forceStepWindow = document.getElementById('forceStepWindow');
    function drawForceChart() {{
      const selected = forceStepWindow?.value || '30';
      const windowSize = selected === 'all' ? null : Number(selected);
      const labelSuffix = selected === 'all' ? 'all steps' : `last ${{selected}}`;
      drawChart('forceChart', movingImages.map(img => {{
        const full = img.force_series.map(r => ({{x:r.step, y:r.force ?? r.max_force}}));
        return {{name:'Image ' + img.image, points: windowSize ? full.slice(-windowSize) : full}};
      }}), {{xLabel:`Ionic step (${{labelSuffix}})`, yLabel:'Max force (eV/A)', xTickCount:5, height:390, thresholds:[{{value:DATA.status.ediffg_abs, label:'EDIFFG', color:'#b54708'}}]}});
    }}
    forceStepWindow?.addEventListener('change', drawForceChart);
    drawForceChart();

    const frontCanvas = document.getElementById('frontCanvas');
    const sideCanvas = document.getElementById('sideCanvas');
    const frontCtx = frontCanvas.getContext('2d');
    const sideCtx = sideCanvas.getContext('2d');
    const frames = DATA.movie_frames || [];
    const slider = document.getElementById('frameSlider');
    const label = document.getElementById('frameLabel');
    const playBtn = document.getElementById('playBtn');
    let frameIndex = 0, timer = null;
    slider.max = Math.max(frames.length - 1, 0);
    const elementStyle = {{
      Mo: {{color:'#8fb3ff', highlight:'#d7e6ff', radius:11, covalent:1.54}},
      S:  {{color:'#ffd166', highlight:'#fff1b8', radius:7.5, covalent:1.05}},
      O:  {{color:'#ef476f', highlight:'#ffc2d0', radius:8.5, covalent:0.66}},
      H:  {{color:'#f8fafc', highlight:'#ffffff', radius:5.5, covalent:0.31}},
      C:  {{color:'#9ca3af', highlight:'#e5e7eb', radius:6.8, covalent:0.76}},
      N:  {{color:'#60a5fa', highlight:'#dbeafe', radius:7.0, covalent:0.71}}
    }};
    const fallbackStyle = {{color:'#e5e7eb', highlight:'#ffffff', radius:6.5, covalent:0.75}};

    const allAtoms = frames.flatMap(f => f.atoms || []);
    const bounds = {{
      x:[Math.min(...allAtoms.map(a=>a.x)), Math.max(...allAtoms.map(a=>a.x))],
      y:[Math.min(...allAtoms.map(a=>a.y)), Math.max(...allAtoms.map(a=>a.y))],
      z:[Math.min(...allAtoms.map(a=>a.z)), Math.max(...allAtoms.map(a=>a.z))]
    }};

    function atomStyle(element) {{
      return elementStyle[element] || fallbackStyle;
    }}

    function buildBonds(atoms) {{
      const bonds = [];
      for (let i = 0; i < atoms.length; i++) {{
        for (let j = i + 1; j < atoms.length; j++) {{
          const a = atoms[i], b = atoms[j];
          const dx = a.x - b.x, dy = a.y - b.y, dz = a.z - b.z;
          const dist = Math.sqrt(dx*dx + dy*dy + dz*dz);
          const cutoff = Math.min(3.15, (atomStyle(a.e).covalent + atomStyle(b.e).covalent) * 1.28);
          if (dist > 0.45 && dist <= cutoff) bonds.push([a, b, dist]);
        }}
      }}
      return bonds;
    }}

    function drawAtom(ctx, x, y, atom) {{
      const style = atomStyle(atom.e);
      const grad = ctx.createRadialGradient(x - style.radius*0.35, y - style.radius*0.4, style.radius*0.15, x, y, style.radius);
      grad.addColorStop(0, style.highlight);
      grad.addColorStop(0.55, style.color);
      grad.addColorStop(1, '#1f2937');
      ctx.beginPath();
      ctx.arc(x, y, style.radius, 0, Math.PI*2);
      ctx.fillStyle = grad;
      ctx.fill();
      ctx.strokeStyle = 'rgba(255,255,255,0.5)';
      ctx.lineWidth = 1.3;
      ctx.stroke();
    }}

    function drawElementLegend(canvas, ctx, atoms) {{
      const elements = [...new Set(atoms.map(a => a.e))].sort();
      const rowH = 24;
      const w = 112;
      const h = elements.length * rowH + 18;
      const x0 = 18;
      const y0 = canvas.height - h - 18;
      ctx.fillStyle = 'rgba(15, 23, 32, 0.78)';
      ctx.strokeStyle = 'rgba(255,255,255,0.18)';
      ctx.lineWidth = 1;
      ctx.fillRect(x0, y0, w, h);
      ctx.strokeRect(x0, y0, w, h);
      ctx.font = '14px system-ui';
      ctx.fillStyle = '#cbd5e1';
      ctx.fillText('Elements', x0 + 12, y0 + 18);
      elements.forEach((element, i) => {{
        const y = y0 + 36 + i * rowH;
        drawAtom(ctx, x0 + 22, y - 5, {{e:element}});
        ctx.fillStyle = '#e5e7eb';
        ctx.fillText(element, x0 + 42, y);
      }});
    }}

    function drawAxisGuide(canvas, ctx, xKey, yKey) {{
      const x0 = canvas.width - 140;
      const y0 = 42;
      const w = 118;
      const h = 116;
      const origin = {{x:x0 + 32, y:y0 + 88}};
      const len = 44;
      const axes = [
        {{label:xKey, dx:len, dy:0, color:'#93c5fd'}},
        {{label:yKey, dx:0, dy:-len, color:'#fda4af'}}
      ];
      ctx.save();
      ctx.fillStyle = 'rgba(15, 23, 32, 0.78)';
      ctx.strokeStyle = 'rgba(255,255,255,0.18)';
      ctx.lineWidth = 1;
      ctx.fillRect(x0, y0, w, h);
      ctx.strokeRect(x0, y0, w, h);
      ctx.fillStyle = '#cbd5e1';
      ctx.font = '14px system-ui';
      ctx.fillText('Axes', x0 + 12, y0 + 20);
      ctx.lineWidth = 4;
      ctx.font = 'bold 18px system-ui';
      axes.forEach(axis => {{
        const x2 = origin.x + axis.dx;
        const y2 = origin.y + axis.dy;
        ctx.strokeStyle = axis.color;
        ctx.fillStyle = axis.color;
        ctx.beginPath();
        ctx.moveTo(origin.x, origin.y);
        ctx.lineTo(x2, y2);
        ctx.stroke();
        const angle = Math.atan2(axis.dy, axis.dx);
        ctx.beginPath();
        ctx.moveTo(x2, y2);
        ctx.lineTo(x2 - 10*Math.cos(angle - 0.45), y2 - 10*Math.sin(angle - 0.45));
        ctx.lineTo(x2 - 10*Math.cos(angle + 0.45), y2 - 10*Math.sin(angle + 0.45));
        ctx.closePath();
        ctx.fill();
        ctx.fillText(axis.label, x2 + (axis.dx ? 8 : -5), y2 + (axis.dy ? -8 : 5));
      }});
      ctx.fillStyle = '#cbd5e1';
      ctx.beginPath();
      ctx.arc(origin.x, origin.y, 4, 0, Math.PI*2);
      ctx.fill();
      ctx.restore();
    }}

    function drawProjection(canvas, ctx, atoms, xKey, yKey, title) {{
      ctx.clearRect(0,0,canvas.width,canvas.height);
      ctx.fillStyle = '#0f1720';
      ctx.fillRect(0,0,canvas.width,canvas.height);
      if (!atoms.length) {{
        ctx.fillStyle = '#cbd5e1'; ctx.fillText('No CONTCAR/POSCAR frames found', 30, 40); return;
      }}
      const [xmin, xmax] = bounds[xKey];
      const [ymin, ymax] = bounds[yKey];
      const scale = Math.min((canvas.width-80)/(xmax-xmin || 1), (canvas.height-80)/(ymax-ymin || 1));
      const ox = (canvas.width - (xmax-xmin)*scale) / 2;
      const oy = (canvas.height - (ymax-ymin)*scale) / 2;
      const point = a => ({{
        x: ox + (a[xKey] - xmin) * scale,
        y: canvas.height - oy - (a[yKey] - ymin) * scale
      }});
      ctx.lineCap = 'round';
      buildBonds(atoms).forEach(([a, b]) => {{
        const pa = point(a), pb = point(b);
        ctx.beginPath();
        ctx.moveTo(pa.x, pa.y);
        ctx.lineTo(pb.x, pb.y);
        ctx.strokeStyle = 'rgba(203,213,225,0.44)';
        ctx.lineWidth = 3.0;
        ctx.stroke();
      }});
      const depthKey = yKey === 'z' ? 'y' : 'z';
      atoms.slice().sort((a,b)=>(a[depthKey] || 0)-(b[depthKey] || 0)).forEach(a => {{
        const p = point(a);
        drawAtom(ctx, p.x, p.y, a);
      }});
      drawElementLegend(canvas, ctx, atoms);
      drawAxisGuide(canvas, ctx, xKey, yKey);
    }}

    function drawFrame(i) {{
      if (!frames.length) {{
        drawProjection(frontCanvas, frontCtx, [], 'x', 'y', '');
        drawProjection(sideCanvas, sideCtx, [], 'x', 'z', '');
        return;
      }}
      const atoms = frames[i].atoms;
      drawProjection(frontCanvas, frontCtx, atoms, 'x', 'y', `Front view: ${{frames[i].comment || ''}}`);
      drawProjection(sideCanvas, sideCtx, atoms, 'x', 'z', `Side view: ${{frames[i].comment || ''}}`);
      label.textContent = `Frame ${{i+1}} / ${{frames.length}}`;
    }}
    slider.addEventListener('input', e => {{ frameIndex = Number(e.target.value); drawFrame(frameIndex); }});
    playBtn.addEventListener('click', () => {{
      if (timer) {{ clearInterval(timer); timer = null; playBtn.textContent = 'Play'; return; }}
      playBtn.textContent = 'Pause';
      timer = setInterval(() => {{
        frameIndex = (frameIndex + 1) % Math.max(frames.length, 1);
        slider.value = frameIndex;
        drawFrame(frameIndex);
      }}, 160);
    }});
    drawFrame(0);
  </script>
</body>
</html>
"""


def collect(workdir: Path) -> dict[str, Any]:
    workdir = workdir.resolve()
    image_dirs = numeric_image_dirs(workdir)
    settings = parse_incar(read_text(workdir / "INCAR"))
    images = build_image_data(image_dirs)
    data = {
        "title": "NEB Monitor",
        "generated_at": datetime.now().strftime("%Y-%m-%d %H:%M:%S"),
        "workdir": str(workdir),
        "settings": settings,
        "images": images,
        "profile": make_relative_profile(images),
        "status": status_from_data(settings, images),
        "movie_frames": build_movie_frames(image_dirs),
        "raw_files": {
            "INCAR": read_text(workdir / "INCAR", 120_000),
            "KPOINTS": read_text(workdir / "KPOINTS", 20_000),
            "POTCAR": summarize_potcar(workdir / "POTCAR"),
            "slurm.sh": read_text(workdir / "slurm.sh", 40_000),
        },
    }
    return data


def main() -> None:
    parser = argparse.ArgumentParser(description="Generate a self-contained VASP NEB HTML monitor.")
    parser.add_argument("neb_working_directory", type=Path, help="Directory containing 00, 01, ... image folders.")
    parser.add_argument("-o", "--output", type=Path, help="Output HTML path. Default: <NEB dir>/neb_monitor.html")
    args = parser.parse_args()

    workdir = args.neb_working_directory
    if not workdir.exists() or not workdir.is_dir():
        raise SystemExit(f"NEB working directory not found: {workdir}")
    output = args.output or (workdir / "neb_monitor.html")
    data = collect(workdir)
    output.write_text(html_page(data), encoding="utf-8")
    print(f"Wrote {output}")
    print(f"Images: {len(data['images'])}, movie frames: {len(data['movie_frames'])}")


if __name__ == "__main__":
    main()
