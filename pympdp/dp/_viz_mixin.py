from __future__ import annotations
from typing import TYPE_CHECKING

if TYPE_CHECKING:
    # Only for type checkers to know what `self` is
    from .dp import DP as _DP

import json
import sys
from pathlib import Path

import numpy as np

try:
    from pympdp.logger import logger
    from pympdp.dubins import dubins_shortest_path, circline, plotdubins
    from pympdp.dp.cell import Cell
except ModuleNotFoundError as exc:  # pragma: no cover - developer convenience
    if exc.name not in {
        "pympdp",
        "pympdp.logger",
        "pympdp.dubins",
        "pympdp.dp",
        "pympdp.dp.cell",
    }:
        raise
    repo_root = Path(__file__).resolve().parents[2]
    if str(repo_root) not in sys.path:
        sys.path.insert(0, str(repo_root))
    from pympdp.logger import logger
    from pympdp.dubins import dubins_shortest_path, circline, plotdubins
    from pympdp.dp.cell import Cell

class _VizMixin:
    def visualize_dp_matrix(self: "_DP", output_path=None, open_in_browser=True, show_optimal_path=False, samples_per_segment=80):
        """Render an interactive HTML dashboard that reflects the current DP matrix."""
        from pathlib import Path
        import html
        import webbrowser

        if not self.dp_matrix:
            raise ValueError("DP matrix is empty – run the solver before visualizing")

        output_path = Path(output_path) if output_path else Path.cwd() / "dp_matrix_visualization.html"

        cell_positions = {}
        for row_idx, row in enumerate(self.dp_matrix):
            for col_idx, cell in enumerate(row):
                cell_positions[id(cell)] = (row_idx, col_idx)

        def to_float(value):
            if value is None:
                return None
            try:
                val = float(value)
            except (TypeError, ValueError):
                return None
            if np.isnan(val) or np.isinf(val):
                return None
            return val

        def resolve_best_path():
            positions = []
            if self.best_path:
                for entry in self.best_path:
                    if isinstance(entry, Cell):
                        pos = cell_positions.get(id(entry))
                    elif isinstance(entry, (tuple, list)) and len(entry) >= 2:
                        pos = (int(entry[0]), int(entry[1]))
                    else:
                        pos = None
                    if pos is not None:
                        positions.append(pos)
            if positions:
                return positions

            if not self.dp_matrix or not self.dp_matrix[-1]:
                return []

            final_row = self.dp_matrix[-1]
            finite_cells = [cell for cell in final_row if np.isfinite(cell.l())]
            candidate = min(finite_cells or final_row, key=lambda cell: cell.l())

            visited = set()
            chain = []
            current = candidate
            while current is not None and id(current) not in visited:
                visited.add(id(current))
                pos = cell_positions.get(id(current))
                if pos is not None:
                    chain.append(pos)
                current = current.prev()

            chain.reverse()
            return chain

        best_path_positions = resolve_best_path()
        best_path_ids = {f"cell-{row}-{col}" for row, col in best_path_positions}
        default_target_id = f"cell-{best_path_positions[-1][0]}-{best_path_positions[-1][1]}" if best_path_positions else ""
        default_target_attr = html.escape(default_target_id)

        def build_optimal_path_payload():
            if not show_optimal_path:
                return None
            try:
                angles = self.best_angles(self.points)
            except Exception as exc:  # pragma: no cover - defensive
                logger.warning("Unable to recover optimal angles for visualization: %s", exc)
                return None
            if not angles or len(angles) != len(self.points):
                return None

            sampled_points: list[tuple[float, float]] = []
            node_points: list[dict[str, object]] = []
            segment_curves = []
            total_length = 0.0
            for idx, (px, py) in enumerate(self.points):
                node_points.append({"x": float(px), "y": float(py), "label": f"#{idx}"})

            for idx in range(len(self.points) - 1):
                x0, y0 = self.points[idx]
                x1, y1 = self.points[idx + 1]
                th0 = angles[idx]
                th1 = angles[idx + 1]
                try:
                    curve, _, seg_lengths = dubins_shortest_path(x0, y0, th0, x1, y1, th1, self.k_max)
                except Exception as exc:  # pragma: no cover - defensive
                    logger.warning("Failed to compute Dubins segment for visualization: %s", exc)
                    curve = None

                segment_samples: list[tuple[float, float]] = []
                if curve is not None:
                    segment_curves.append(curve)
                    total_length += float(sum(seg_lengths))
                    for arc in (curve.a1, curve.a2, curve.a3):
                        if arc.L <= 0:
                            continue
                        steps = max(4, int(max(arc.L * samples_per_segment, samples_per_segment)))
                        for step in range(steps):
                            s = arc.L * step / max(steps - 1, 1)
                            x, y, _ = circline(s, arc.x0, arc.y0, arc.th0, arc.k)
                            segment_samples.append((x, y))
                    segment_samples.append((curve.a3.xf, curve.a3.yf))
                else:
                    segment_samples = [(float(x0), float(y0)), (float(x1), float(y1))]

                if sampled_points:
                    segment_samples = segment_samples[1:]
                sampled_points.extend((float(px), float(py)) for px, py in segment_samples)

            if not sampled_points:
                return None

            xs = [pt[0] for pt in sampled_points]
            ys = [pt[1] for pt in sampled_points]
            pad_source = max((max(xs) - min(xs)), (max(ys) - min(ys)))
            pad = (pad_source if pad_source > 0 else 1.0) * 0.05
            bounds = {
                "min_x": min(xs) - pad,
                "max_x": max(xs) + pad,
                "min_y": min(ys) - pad,
                "max_y": max(ys) + pad,
            }
            image_data_url = None
            if segment_curves:
                try:
                    from io import BytesIO
                    import base64
                    import matplotlib.pyplot as plt

                    fig, ax = plt.subplots(figsize=(6, 6))
                    ax.set_aspect('equal', 'box')
                    plt.sca(ax)
                    for curve in segment_curves:
                        plotdubins(curve, show=False)
                    ax.set_xlabel('x')
                    ax.set_ylabel('y')
                    ax.set_title('Optimal Dubins Path')
                    for px, py in self.points:
                        ax.plot(px, py, 'ko', markersize=4)
                    buf = BytesIO()
                    fig.savefig(buf, format='png', bbox_inches='tight')
                    buf.seek(0)
                    image_data_url = 'data:image/png;base64,' + base64.b64encode(buf.getvalue()).decode('ascii')
                except Exception as exc:  # pragma: no cover - defensive
                    logger.warning("Unable to render Dubins plot: %s", exc)
                finally:
                    try:
                        plt.close(fig)
                    except Exception:
                        pass

            return {
                "path": [{"x": x, "y": y} for x, y in sampled_points],
                "nodes": node_points,
                "bounds": bounds,
                "totalLength": total_length,
                "imageDataUrl": image_data_url,
            }

        path_payload = build_optimal_path_payload()
        path_payload_attr = html.escape(json.dumps(path_payload)) if path_payload else ""
        path_panel_class = "path-panel" + ("" if path_payload else " hidden")

        rows_html = []
        for row_idx, (row, point) in enumerate(zip(self.dp_matrix, self.points)):
            px = to_float(point[0])
            py = to_float(point[1])
            coord_label = f"({px:.2f}, {py:.2f})" if px is not None and py is not None else "(?, ?)"
            row_cells = [f'<th class="row-header">#{row_idx}<br><span class="coord">{html.escape(coord_label)}</span></th>']

            for col_idx, cell in enumerate(row):
                cell_id = f"cell-{row_idx}-{col_idx}"
                classes = ["cell"]
                if cell_id in best_path_ids:
                    classes.append("best")

                angle_val = to_float(cell.th())
                angle_deg_val = float(np.degrees(angle_val)) if angle_val is not None else None
                length_val = to_float(cell.l())

                angle_label = "θ —" if angle_val is None else f"θ {angle_val:.3f} rad"
                angle_deg_label = "" if angle_deg_val is None else f"{angle_deg_val:.1f}°"
                length_label = "L —" if length_val is None else f"L {length_val:.3f}"

                prev_id = ""
                prev_cell = cell.prev()
                if prev_cell is not None:
                    pos = cell_positions.get(id(prev_cell))
                    if pos is not None:
                        prev_id = f"cell-{pos[0]}-{pos[1]}"

                attrs = {
                    "id": cell_id,
                    "class": " ".join(classes),
                    "data-row": str(row_idx),
                    "data-col": str(col_idx),
                    "data-point-x": "" if px is None else f"{px:.6f}",
                    "data-point-y": "" if py is None else f"{py:.6f}",
                    "data-angle": "" if angle_val is None else f"{angle_val:.6f}",
                    "data-angle-deg": "" if angle_deg_val is None else f"{angle_deg_val:.6f}",
                    "data-length": "" if length_val is None else f"{length_val:.6f}",
                    "data-prev-id": prev_id,
                }
                attr_html = " ".join(f'{key}="{html.escape(str(value))}"' for key, value in attrs.items())

                content_bits = [f'<span class="angle">{html.escape(angle_label)}</span>']
                if angle_deg_label:
                    content_bits.append(f'<span class="angle-deg">{html.escape(angle_deg_label)}</span>')
                content_bits.append(f'<span class="length">{html.escape(length_label)}</span>')

                row_cells.append(f'<td><button {attr_html}>{"".join(content_bits)}</button></td>')

            rows_html.append('<tr>' + ''.join(row_cells) + '</tr>')

        html_content = f"""<!DOCTYPE html>
    <html lang=\"en\">
    <head>
        <meta charset=\"utf-8\">
        <title>DP Matrix Visualization</title>
        <style>
            :root {{
                font-family: 'Segoe UI', Roboto, Helvetica, Arial, sans-serif;
                color: #0b1e34;
                background: #f8fbff;
                --cell-scale: 1;
                --cell-min-width: 160px;
                --cell-padding-y: 0.9rem;
                --cell-padding-x: 1rem;
            }}
            body {{
                margin: 2rem;
            }}
            h1 {{
                margin-bottom: 0.25rem;
            }}
            p.meta {{
                margin-top: 0;
                color: #4a6278;
            }}
            table.dp-table {{
                border-collapse: separate;
                border-spacing: 0;
                width: 100%;
                background: #ffffff;
                box-shadow: 0 6px 24px rgba(12, 30, 51, 0.12);
                border-radius: 12px;
                overflow: hidden;
            }}
            th.row-header {{
                padding: 1rem 1.25rem;
                background: linear-gradient(135deg, #0b7285, #228be6);
                color: #ffffff;
                width: 190px;
                text-align: left;
                font-weight: 600;
                vertical-align: top;
            }}
            th.row-header .coord {{
                display: block;
                font-size: 0.85rem;
                opacity: 0.9;
                margin-top: 0.35rem;
            }}
            td {{
                padding: calc(0.375rem * var(--cell-scale));
            }}
            button.cell {{
                all: unset;
                display: flex;
                flex-direction: column;
                gap: 0.25rem;
                align-items: flex-start;
                justify-content: center;
                min-width: calc(var(--cell-min-width) * var(--cell-scale));
                padding: calc(var(--cell-padding-y) * var(--cell-scale)) calc(var(--cell-padding-x) * var(--cell-scale));
                border-radius: 10px;
                background: #edf2ff;
                border: 2px solid transparent;
                cursor: pointer;
                transition: transform 0.15s ease, box-shadow 0.2s ease, border 0.2s ease;
                position: relative;
            }}
            button.cell.best {{
                background: linear-gradient(135deg, rgba(224, 243, 255, 0.95), #d0ebff);
                border-color: rgba(34, 139, 230, 0.6);
                box-shadow: 0 12px 26px rgba(34, 139, 230, 0.22);
            }}
            button.cell.best .angle,
            button.cell.best .angle-deg,
            button.cell.best .length {{
                color: #07364a;
            }}
            button.cell:hover {{
                transform: translateY(-2px);
                box-shadow: 0 10px 20px rgba(34, 139, 230, 0.20);
            }}
            button.cell.hover-source {{
                border-color: rgba(34, 139, 230, 0.55);
                box-shadow: 0 8px 18px rgba(34, 139, 230, 0.25);
            }}
            button.cell.active {{
                background: linear-gradient(135deg, #0b7285, #1b9aaa);
                color: #ffffff;
                border-color: #0b7285;
                box-shadow: 0 12px 24px rgba(11, 114, 133, 0.35);
            }}
            button.cell.trail {{
                background: rgba(34, 139, 230, 0.12);
                border-color: rgba(34, 139, 230, 0.35);
                color: #0b1e34;
            }}
            button.cell.hover-prev {{
                border-color: rgba(11, 114, 133, 0.6);
                box-shadow: 0 10px 20px rgba(11, 114, 133, 0.25);
            }}
            .back-arrow {{
                position: absolute;
                left: 0;
                top: 0;
                height: 4px;
                background: linear-gradient(90deg, rgba(11, 114, 133, 0.0), rgba(11, 114, 133, 0.95));
                border-radius: 999px;
                transform-origin: 0 50%;
                pointer-events: none;
                z-index: 1000;
                display: none;
            }}
            .back-arrow-head {{
                position: absolute;
                top: 50%;
                width: 0;
                height: 0;
                border-top: 6px solid transparent;
                border-bottom: 6px solid transparent;
                border-left: 12px solid rgba(11, 114, 133, 0.95);
                transform: translateY(-50%);
            }}
            button.cell .angle {{
                font-size: calc(1.05rem * var(--cell-scale));
                font-weight: 600;
            }}
            button.cell .angle-deg {{
                font-size: calc(0.85rem * var(--cell-scale));
                opacity: 0.85;
            }}
            button.cell .length {{
                font-size: calc(0.9rem * var(--cell-scale));
                font-weight: 500;
            }}
            .controls {{
                margin: 1.25rem 0 1.5rem;
                display: flex;
                align-items: center;
                gap: 1rem;
                color: #1c3144;
            }}
            .controls label {{
                font-weight: 600;
                font-size: 0.95rem;
            }}
            .controls input[type="range"] {{
                flex: 1 1 220px;
                accent-color: #0b7285;
            }}
            .controls .control-value {{
                font-variant-numeric: tabular-nums;
                font-weight: 600;
                min-width: 3rem;
            }}
            .details {{
                margin-top: 1.75rem;
                padding: 1.5rem;
                background: #ffffff;
                border-radius: 12px;
                box-shadow: 0 6px 18px rgba(12, 30, 51, 0.08);
            }}
            .details h2 {{
                margin-top: 0;
                margin-bottom: 0.75rem;
            }}
            .details ol {{
                margin: 0;
                padding-left: 1.35rem;
                color: #1c3144;
            }}
            .details li {{
                margin-bottom: 0.5rem;
            }}
            .details li:last-child {{
                margin-bottom: 0;
            }}
            .details em {{
                color: #587089;
            }}
            .path-panel {{
                margin-top: 1.75rem;
                padding: 1.5rem;
                background: #ffffff;
                border-radius: 12px;
                box-shadow: 0 6px 18px rgba(12, 30, 51, 0.08);
            }}
            .path-panel.hidden {{
                display: none;
            }}
            .path-panel img {{
                width: 100%;
                max-height: 360px;
                object-fit: contain;
                display: block;
            }}
            .path-panel h2 {{
                margin-top: 0;
                margin-bottom: 0.75rem;
            }}
            .path-panel .path-meta {{
                margin-bottom: 0.75rem;
                color: #4a6278;
            }}
        </style>
    </head>
    <body data-default-target=\"{default_target_attr}\">
        <h1>Dynamic Programming Matrix</h1>
        <p class=\"meta\">Click a node to highlight the stored optimal chain. Hover to preview the predecessor.</p>
        <div class=\"controls\">
            <label for=\"cell-scale\">Cell size</label>
            <input type=\"range\" id=\"cell-scale\" min=\"0.6\" max=\"1.6\" step=\"0.1\" value=\"1\">
            <span class=\"control-value\" id=\"cell-scale-value\">1.0x</span>
        </div>
        <table class=\"dp-table\">
            <tbody>
                {''.join(rows_html)}
            </tbody>
        </table>
        <div class=\"details\" id=\"details\"><em>Select a cell to explore its optimal path.</em></div>
        <div class=\"{path_panel_class}\" id=\"path-panel\" data-path-payload=\"{path_payload_attr}\">
            <h2>Optimal Path</h2>
            <p class=\"path-meta\" id=\"path-meta\">Visualizes the current best sequence of Dubins segments.</p>
            <img id=\"path-image\" alt=\"Optimal Dubins path\" />
        </div>
        <script>
            (function() {{
                const cells = Array.from(document.querySelectorAll('button.cell'));
                const details = document.getElementById('details');
                const pathPanel = document.getElementById('path-panel');
                const defaultTargetId = document.body.dataset.defaultTarget;
                const scaleControl = document.getElementById('cell-scale');
                const scaleValueLabel = document.getElementById('cell-scale-value');
                let hoverPrevCell = null;
                let hoverSourceCell = null;
                let arrowSourceCell = null;
                let arrowTargetCell = null;
                let backArrowEl = null;

                function clearHover() {{
                    if (hoverPrevCell) {{
                        hoverPrevCell.classList.remove('hover-prev');
                        hoverPrevCell = null;
                    }}
                    if (hoverSourceCell) {{
                        hoverSourceCell.classList.remove('hover-source');
                        hoverSourceCell = null;
                    }}
                    removeBackArrow();
                }}

                function collectChain(start) {{
                    const chain = [];
                    const seen = new Set();
                    let current = start;
                    while (current && !seen.has(current.id)) {{
                        chain.push(current);
                        seen.add(current.id);
                        const prevId = current.dataset.prevId;
                        if (!prevId) break;
                        current = document.getElementById(prevId);
                    }}
                    return chain;
                }}

                function formatNumber(value, digits) {{
                    if (!value) return '—';
                    const numeric = Number(value);
                    if (!Number.isFinite(numeric)) return '—';
                    return numeric.toFixed(digits);
                }}

                function formatAngle(value) {{
                    const base = formatNumber(value, 3);
                    return base === '—' ? base : base + ' rad';
                }}

                function formatDegrees(value) {{
                    const base = formatNumber(value, 1);
                    return base === '—' ? '' : base + '°';
                }}

                function ensureBackArrow() {{
                    if (!backArrowEl) {{
                        backArrowEl = document.createElement('div');
                        backArrowEl.className = 'back-arrow';
                        const head = document.createElement('span');
                        head.className = 'back-arrow-head';
                        backArrowEl.appendChild(head);
                        document.body.appendChild(backArrowEl);
                    }}
                    return backArrowEl;
                }}

                function updateBackArrowPosition() {{
                    if (!arrowSourceCell || !arrowTargetCell) {{
                        return;
                    }}
                    const arrowEl = ensureBackArrow();
                    const sourceRect = arrowSourceCell.getBoundingClientRect();
                    const targetRect = arrowTargetCell.getBoundingClientRect();
                    const x1 = sourceRect.left + sourceRect.width / 2 + window.scrollX;
                    const y1 = sourceRect.top + sourceRect.height / 2 + window.scrollY;
                    const x2 = targetRect.left + targetRect.width / 2 + window.scrollX;
                    const y2 = targetRect.top + targetRect.height / 2 + window.scrollY;
                    const dx = x2 - x1;
                    const dy = y2 - y1;
                    const distance = Math.hypot(dx, dy);

                    if (!Number.isFinite(distance) || distance < 2) {{
                        removeBackArrow();
                        return;
                    }}

                    const headLength = 14;
                    const shaftLength = Math.max(distance - headLength, 0);
                    arrowEl.style.width = shaftLength + 'px';
                    arrowEl.style.transform = 'translate(' + x1 + 'px, ' + y1 + 'px) rotate(' + Math.atan2(dy, dx) + 'rad)';
                    arrowEl.style.display = 'block';

                    const head = arrowEl.querySelector('.back-arrow-head');
                    if (head) {{
                        head.style.transform = 'translate(' + shaftLength + 'px, -50%)';
                    }}
                }}

                function showBackArrow(source, previous) {{
                    arrowSourceCell = source;
                    arrowTargetCell = previous;
                    ensureBackArrow();
                    updateBackArrowPosition();
                    requestAnimationFrame(updateBackArrowPosition);
                }}

                function removeBackArrow() {{
                    arrowSourceCell = null;
                    arrowTargetCell = null;
                    if (backArrowEl) {{
                        backArrowEl.style.display = 'none';
                    }}
                }}

                function updateDetails(chain) {{
                    if (!chain.length) {{
                        details.innerHTML = '<em>Select a cell to explore its optimal path.</em>';
                        return;
                    }}

                    const ordered = chain.slice().reverse();
                    const rows = ordered.map((cell, index) => {{
                        const label = index === ordered.length - 1 ? 'Selected' : 'Step ' + (index + 1);
                        const row = cell.dataset.row ?? '—';
                        const col = cell.dataset.col ?? '—';
                        const px = formatNumber(cell.dataset.pointX, 2);
                        const py = formatNumber(cell.dataset.pointY, 2);
                        const angle = formatAngle(cell.dataset.angle);
                        const angleDeg = formatDegrees(cell.dataset.angleDeg);
                        const length = formatNumber(cell.dataset.length, 3);
                        const degPart = angleDeg ? ' / ' + angleDeg : '';
                        return '<li><strong>' + label + '</strong> → row ' + row + ', col ' + col + ', point (' + px + ', ' + py + '), θ ' + angle + degPart + ', L ' + length + '</li>';
                    }}).join('');

                    details.innerHTML = '<h2>Path Details</h2><ol>' + rows + '</ol>';
                }}

                function activate(target) {{
                    cells.forEach(btn => btn.classList.remove('active', 'trail'));
                    clearHover();
                    if (!target) {{
                        updateDetails([]);
                        return;
                    }}
                    target.classList.add('active');
                    const chain = collectChain(target);
                    chain.slice(1).forEach(cell => cell.classList.add('trail'));
                    updateDetails(chain);
                }}

                function applyScale(value) {{
                    const numeric = Number(value);
                    const clamped = Number.isFinite(numeric) ? Math.min(Math.max(numeric, 0.6), 2) : 1;
                    document.documentElement.style.setProperty('--cell-scale', clamped.toString());
                    if (scaleValueLabel) {{
                        scaleValueLabel.textContent = clamped.toFixed(1) + 'x';
                    }}
                    if (arrowSourceCell && arrowTargetCell) {{
                        updateBackArrowPosition();
                        requestAnimationFrame(updateBackArrowPosition);
                    }}
                    renderOptimalPathPanel();
                }}

                function renderOptimalPathPanel() {{
                    if (!pathPanel) {{
                        return;
                    }}
                    const payloadRaw = pathPanel.dataset.pathPayload;
                    if (!payloadRaw) {{
                        pathPanel.classList.add('hidden');
                        return;
                    }}
                    let payload;
                    try {{
                        payload = JSON.parse(payloadRaw);
                    }} catch (error) {{
                        console.warn('Unable to parse optimal path payload', error);
                        pathPanel.classList.add('hidden');
                        return;
                    }}
                    const imageEl = pathPanel.querySelector('#path-image');
                    const metaEl = pathPanel.querySelector('#path-meta');
                    const imageUrl = payload.imageDataUrl || '';
                    if (!payload || !imageUrl) {{
                        pathPanel.classList.add('hidden');
                        return;
                    }}
                    pathPanel.classList.remove('hidden');
                    if (imageEl) {{
                        imageEl.src = imageUrl;
                    }}
                    if (metaEl) {{
                        const lengthText = typeof payload.totalLength === 'number'
                            ? `Length ≈ ${{payload.totalLength.toFixed(4)}}`
                            : 'Visualizes the current best sequence of Dubins segments.';
                        metaEl.textContent = lengthText;
                    }}
                }}

                if (scaleControl) {{
                    applyScale(scaleControl.value || '1');
                    scaleControl.addEventListener('input', event => {{
                        applyScale(event.currentTarget.value);
                    }});
                }} else {{
                    renderOptimalPathPanel();
                }}

                ['scroll', 'resize'].forEach(eventName => {{
                    window.addEventListener(eventName, () => {{
                        if (arrowSourceCell && arrowTargetCell) {{
                            updateBackArrowPosition();
                        }}
                    }}, {{ passive: true }});
                }});

                renderOptimalPathPanel();

            cells.forEach(cell => {{
                    cell.addEventListener('mouseenter', event => {{
                        const target = event.currentTarget;
                        clearHover();
                        target.classList.add('hover-source');
                        hoverSourceCell = target;
                        const prevId = target.dataset.prevId;
                        if (prevId) {{
                            const prevCell = document.getElementById(prevId);
                            if (prevCell) {{
                                prevCell.classList.add('hover-prev');
                                hoverPrevCell = prevCell;
                                showBackArrow(target, prevCell);
                            }} else {{
                                removeBackArrow();
                            }}
                        }} else {{
                            removeBackArrow();
                        }}
                    }});

                    cell.addEventListener('mouseleave', () => {{
                        clearHover();
                    }});

                    cell.addEventListener('click', event => {{
                        event.stopPropagation();
                        activate(event.currentTarget);
                    }});
                }});

                document.addEventListener('click', event => {{
                    if (!event.target.closest('button.cell')) {{
                        activate(null);
                    }}
                }});

                document.addEventListener('keydown', event => {{
                    if (event.key === 'Escape') {{
                        activate(null);
                    }}
                }});

                if (defaultTargetId) {{
                    const defaultCell = document.getElementById(defaultTargetId);
                    if (defaultCell) {{
                        activate(defaultCell);
                    }}
                }}
            }})();
        </script>
    </body>
    </html>
    """

        output_path.write_text(html_content, encoding="utf-8")

        if open_in_browser:
            webbrowser.open(output_path.as_uri())

        logger.info("DP matrix visualization saved to %s", output_path)
