"""Build the v1 catalog and its black-and-white fidelity figure."""
import argparse
import csv
import html
import json
from pathlib import Path

from omicsbench.registry import registry


ROOT = Path(__file__).resolve().parents[1]
METRICS = (
    ("spearman_logfc", "Spearman effect correlation"),
    ("top_k_jaccard", "Top-50 effect overlap"),
    ("distance_correlation", "Sample-distance correlation"),
    ("sign_concordance", "Effect-sign agreement"),
)


def read_json(path: Path):
    return json.loads(path.read_text(encoding="utf-8"))


def collect() -> list[dict]:
    rows = []
    for model, folder in registry(ROOT).values():
        if model.kind != "real" or model.status != "validated":
            continue
        validation = read_json(folder / model.validation.profile)
        evidence = read_json(folder / validation["deseq2_evidence"])
        reference = read_json(folder / "reference.json")
        attribution = read_json(folder / "attribution.json")
        rows.append(
            {
                "id": model.id,
                "title": model.title,
                "object_key": attribution["object_key"],
                "source_accession": model.source.accession,
                "organism": model.organism,
                "taxon_id": model.taxon_id,
                "archetype": model.archetype,
                "samples": len(model.samples),
                "pocket_features": evidence["design"]["pocket_features"],
                "contrast": evidence["design"]["contrast"],
                "reference": f"{reference['provider']} {reference['release']} {reference['assembly']}",
                "rights": model.rights.status,
                "license": "CC-BY-4.0",
                "metrics": evidence["metrics"],
                "thresholds": evidence["thresholds"],
            }
        )
    return rows


def write_catalog(rows: list[dict]) -> None:
    catalog = ROOT / "catalog"
    catalog.mkdir(exist_ok=True)
    payload = {
        "schema_version": "1.0",
        "release": "1.0.0",
        "biological_objects": len(rows),
        "source_studies": len({row["source_accession"] for row in rows}),
        "objects": rows,
    }
    (catalog / "collection.json").write_text(
        json.dumps(payload, indent=2) + "\n", encoding="utf-8", newline="\n"
    )
    fields = (
        "id",
        "title",
        "object_key",
        "source_accession",
        "organism",
        "taxon_id",
        "archetype",
        "samples",
        "pocket_features",
        "contrast",
        "reference",
        "rights",
        "license",
        "spearman_logfc",
        "top_k_jaccard",
        "distance_correlation",
        "sign_concordance",
    )
    with (catalog / "collection.tsv").open("w", encoding="utf-8", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=fields, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        for row in rows:
            flat = {key: row.get(key) for key in fields}
            flat.update(row["metrics"])
            writer.writerow(flat)


def write_figure(rows: list[dict]) -> None:
    width, height = 1200, 930
    panel_width, panel_height = 570, 370
    x_min, x_max = 0.55, 1.0
    parts = [
        f'<svg xmlns="http://www.w3.org/2000/svg" width="{width}" height="{height}" viewBox="0 0 {width} {height}">',
        '<rect width="100%" height="100%" fill="white"/>',
        '<style>text{fill:#111;font-family:Arial,sans-serif}.title{font-size:26px;font-weight:700}.subtitle{font-size:15px;fill:#444}.panel{font-size:17px;font-weight:700}.label{font-size:12px}.tick{font-size:11px;fill:#555}.value{font-size:10px;font-family:Consolas,monospace}.note{font-size:12px;fill:#444}</style>',
        '<text x="52" y="42" class="title">DESeq2 fidelity of the v1 biological collection</text>',
        '<text x="52" y="67" class="subtitle">Twelve pocket objects compared with their full source matrices; dots show the recorded metric value.</text>',
    ]
    for panel_index, (metric, label) in enumerate(METRICS):
        column, row_index = panel_index % 2, panel_index // 2
        left = 55 + column * 590
        top = 105 + row_index * 385
        plot_left = left + 118
        plot_right = left + panel_width - 20
        plot_top = top + 45
        plot_bottom = top + panel_height - 28
        parts.append(f'<text x="{left}" y="{top + 19}" class="panel">{html.escape(label)}</text>')
        for tick in (0.6, 0.7, 0.8, 0.9, 1.0):
            x = plot_left + (tick - x_min) / (x_max - x_min) * (plot_right - plot_left)
            parts.append(f'<line x1="{x:.1f}" y1="{plot_top}" x2="{x:.1f}" y2="{plot_bottom}" stroke="#ddd"/>')
            parts.append(f'<text x="{x:.1f}" y="{plot_bottom + 18}" text-anchor="middle" class="tick">{tick:.1f}</text>')
        threshold = rows[0]["thresholds"][metric]
        threshold_x = plot_left + (threshold - x_min) / (x_max - x_min) * (plot_right - plot_left)
        parts.append(
            f'<line x1="{threshold_x:.1f}" y1="{plot_top}" x2="{threshold_x:.1f}" y2="{plot_bottom}" stroke="#111" stroke-width="1.5" stroke-dasharray="5 4"/>'
        )
        parts.append(
            f'<text x="{threshold_x - 4:.1f}" y="{plot_top - 7}" text-anchor="end" class="tick">minimum {threshold:.2f}</text>'
        )
        spacing = (plot_bottom - plot_top) / len(rows)
        for index, item in enumerate(rows):
            y = plot_top + spacing * (index + 0.5)
            value = item["metrics"][metric]
            x = plot_left + (value - x_min) / (x_max - x_min) * (plot_right - plot_left)
            parts.append(f'<text x="{plot_left - 9}" y="{y + 4:.1f}" text-anchor="end" class="label">{item["id"]}</text>')
            parts.append(f'<line x1="{threshold_x:.1f}" y1="{y:.1f}" x2="{x:.1f}" y2="{y:.1f}" stroke="#999"/>')
            parts.append(f'<circle cx="{x:.1f}" cy="{y:.1f}" r="4.3" fill="#222"/>')
            anchor = "end" if x > plot_right - 35 else "start"
            offset = -7 if anchor == "end" else 7
            parts.append(f'<text x="{x + offset:.1f}" y="{y - 6:.1f}" text-anchor="{anchor}" class="value">{value:.3f}</text>')
    parts.extend(
        [
            '<text x="55" y="895" class="note">All metrics pass their predeclared minimum. Values come from each object’s DESeq2 1.50.2 evidence record.</text>',
            '<text x="55" y="915" class="note">Exact values and source accessions are in catalog/collection.tsv.</text>',
            "</svg>",
        ]
    )
    output = ROOT / "figures/v1-fidelity.svg"
    output.parent.mkdir(exist_ok=True)
    output.write_text("\n".join(parts) + "\n", encoding="utf-8", newline="\n")


def write_preview(rows: list[dict]) -> Path:
    from PIL import Image, ImageDraw, ImageFont

    width, height = 1200, 930
    image = Image.new("RGB", (width, height), "white")
    draw = ImageDraw.Draw(image)

    def font(size: int, bold: bool = False):
        names = ["arialbd.ttf" if bold else "arial.ttf", "DejaVuSans-Bold.ttf" if bold else "DejaVuSans.ttf"]
        for name in names:
            try:
                return ImageFont.truetype(name, size)
            except OSError:
                pass
        return ImageFont.load_default()

    draw.text((52, 20), "DESeq2 fidelity of the v1 biological collection", fill="#111", font=font(26, True))
    draw.text(
        (52, 55),
        "Twelve pocket objects compared with their full source matrices; dots show the recorded metric value.",
        fill="#444",
        font=font(15),
    )
    x_min, x_max = 0.55, 1.0
    for panel_index, (metric, label) in enumerate(METRICS):
        column, row_index = panel_index % 2, panel_index // 2
        left = 55 + column * 590
        top = 105 + row_index * 385
        plot_left, plot_right = left + 118, left + 550
        plot_top, plot_bottom = top + 45, top + 342
        draw.text((left, top), label, fill="#111", font=font(17, True))
        for tick in (0.6, 0.7, 0.8, 0.9, 1.0):
            x = plot_left + (tick - x_min) / (x_max - x_min) * (plot_right - plot_left)
            draw.line((x, plot_top, x, plot_bottom), fill="#ddd", width=1)
            draw.text((x - 10, plot_bottom + 5), f"{tick:.1f}", fill="#555", font=font(11))
        threshold = rows[0]["thresholds"][metric]
        threshold_x = plot_left + (threshold - x_min) / (x_max - x_min) * (plot_right - plot_left)
        for y in range(plot_top, plot_bottom, 9):
            draw.line((threshold_x, y, threshold_x, min(y + 5, plot_bottom)), fill="#111", width=2)
        draw.text((threshold_x - 88, plot_top - 17), f"minimum {threshold:.2f}", fill="#555", font=font(11))
        spacing = (plot_bottom - plot_top) / len(rows)
        for index, item in enumerate(rows):
            y = plot_top + spacing * (index + 0.5)
            value = item["metrics"][metric]
            x = plot_left + (value - x_min) / (x_max - x_min) * (plot_right - plot_left)
            draw.text((plot_left - 80, y - 7), item["id"], fill="#111", font=font(12))
            draw.line((threshold_x, y, x, y), fill="#999", width=1)
            draw.ellipse((x - 4, y - 4, x + 4, y + 4), fill="#222")
            value_x = x - 38 if x > plot_right - 35 else x + 7
            draw.text((value_x, y - 14), f"{value:.3f}", fill="#111", font=font(10))
    draw.text(
        (55, 884),
        "All metrics pass their predeclared minimum. Values come from each object's DESeq2 1.50.2 evidence record.",
        fill="#444",
        font=font(12),
    )
    draw.text((55, 904), "Exact values and source accessions are in catalog/collection.tsv.", fill="#444", font=font(12))
    output = ROOT / "build/v1-fidelity-preview.png"
    output.parent.mkdir(exist_ok=True)
    image.save(output)
    return output


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--preview-png", action="store_true", help="Render a local PNG for visual inspection.")
    args = parser.parse_args()
    collection = collect()
    write_catalog(collection)
    write_figure(collection)
    if args.preview_png:
        print(f"Preview: {write_preview(collection)}")
    print(f"Wrote {len(collection)} catalog records and the v1 fidelity figure.")
