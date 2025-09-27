import time
import numpy as np

from logger import logger
from cell import Cell
from utility import circles
from dubins import dubins_shortest_path

class DP:
    def __init__(self, points, fixed_angles, k_max, discretizations, refinements, def_thetas=None):
        if def_thetas:
            for i in range(len(def_thetas)):
                if fixed_angles[i] and def_thetas[i] is None:
                    raise AssertionError("if def_thetas is provided, all angles that are fixed must have a value in def_thetas")

        assert len(points) == len(fixed_angles), "Length of points and fixed_angles must be the same"
        assert def_thetas is None or len(def_thetas) == len(points), "Length of def_thetas must be the same as points if provided"
        assert def_thetas is None or all((not fixed) or (theta is not None) for fixed, theta in zip(fixed_angles, def_thetas)), "if def_thetas is provided, all angles that are fixed must have a value in def_thetas"
        assert def_thetas is not None or (def_thetas is None and all(not fixed_angles[i] for i in range(len(fixed_angles)))), "if def_thetas is None, all fixed_angles must be False"

        self.points = points
        self.fixed_angles = fixed_angles
        self.k_max = k_max
        self.discretizations = discretizations
        self.refinements = refinements

        self.def_thetas = def_thetas if def_thetas is not None else [0.0 for _ in points]

        self.dp_matrix = [[Cell() for _ in range(discretizations)] for _ in range(len(points))]
        self.best_path = None

        self.set_sampling_angles(hrange=2*np.pi)


    def set_first_angle(self, angle):
        """
        Set the angle of the first point in the DP matrix.
        :param angle: Angle in radians
        """
        if self.def_thetas is None:
            self.def_thetas = [0.0 for _ in self.points]
        
        self.def_thetas[0] = angle
        self.fixed_angles[0] = True
        self.dp_matrix[0] = [Cell(angle, 0, None)]

    
    def set_final_angle(self, angle):
        """
        Set the angle of the final point in the DP matrix.
        :param angle: Angle in radians
        """
        if self.def_thetas is None:
            self.def_thetas = [0.0 for _ in self.points]
        
        self.def_thetas[-1] = angle
        self.fixed_angles[-1] = True
        self.dp_matrix[-1] = [Cell(angle, 0, None)]


    def _reset_matrix(self):
        """
        Reset the dynamic programming matrix to an empty state.
        """
        for i in range(len(self.dp_matrix)):
            for j in range(len(self.dp_matrix[i])):
                self.dp_matrix[i][j] = Cell()


    def set_sampling_angles(self, hrange=2*np.pi):
        """
        Set the sampling angles for dynamic programming.
        :param hrange: Range of angles to sample. Default is 2pi for the first round, 3/2pi for subsequent rounds.
        """
        assert self.def_thetas is not None
        sampling_angles = [np.linspace(th-hrange/2, th+hrange/2, self.discretizations, endpoint=(hrange != 2*np.pi)).tolist() for th in self.def_thetas]

        for i in range(len(self.points)):
            if self.fixed_angles[i]:
                self.dp_matrix[i] = [Cell(self.def_thetas[i], 0, None)]
            else:
                self.dp_matrix[i] = [Cell(th, 0, None) for th in sampling_angles[i]]
                if i > 0:
                    th_prev, th_curr = self.guess_initial_angles(i)

                    # Adding guessed angles to the previous point if they are not already present
                    if not self.fixed_angles[i - 1]:
                        tmp_prev_angles = [cell.th() for cell in self.dp_matrix[i - 1]]
                        for th in th_prev:
                            if not any(np.isclose(th.th(), angle) for angle in tmp_prev_angles):
                                self.dp_matrix[i - 1].append(th)

                    # Adding guessed angles to the current point if they are not already present
                    tmp_curr_angles = [cell.th() for cell in self.dp_matrix[i]]
                    for th in th_curr:
                        if not any(np.isclose(th.th(), angle) for angle in tmp_curr_angles):
                            self.dp_matrix[i].append(th)


    def guess_initial_angles(self, i):
        """
        Guess initial angles for dynamic programming based on previous and current angles.
        :param i: Index of the current point
        """
        if i == 0:
            return [], []
        th = np.arctan2(self.points[i][1] - self.points[i - 1][1], self.points[i][0] - self.points[i - 1][0])
        th_prev = [Cell(th)]
        th_curr = [Cell(th)]

        XC, YX = circles(self.points[i - 1][0], self.points[i - 1][1], self.points[i][0], self.points[i][1], 1.0 / self.k_max)

        for xc, yc in zip(XC, YX):
            thp = np.arctan2(yc - self.points[i - 1][1], xc - self.points[i - 1][0])
            thc = np.arctan2(self.points[i][1] - yc, self.points[i][0] - xc)
            if np.isclose(thp, th_prev[-1].th()):
                th_prev.append(Cell(thp))
            if np.isclose(thc, th_curr[-1].th()):
                th_curr.append(Cell(thc))

        logger.info(f"curr {[cell.th() for cell in th_curr]}")
        logger.info(f"prev {[cell.th() for cell in th_prev]}")

        return th_prev, th_curr


    def solve_dp(self):
        # First round is on the house
        self.solve_dp_inner()

        # Extract optimal angles
        opt_angles = self.best_angles(self.points)
        logger.debug(f"Initial optimal angles: {opt_angles}")

        # Refinement rounds
        for r in range(self.refinements):
            logger.info(f"Refinement round {r+1}/{self.refinements}")
            self._reset_matrix()
            self.def_thetas = opt_angles
            self.set_sampling_angles(hrange = np.pi)
            self.solve_dp_inner()
            opt_angles = self.best_angles(self.points)
            logger.debug(f"Optimal angles for refinement {r+1}: {opt_angles}")

        logger.info(f"Optimal angles: {opt_angles}")

        return opt_angles


    def solve_dp_inner(self):
        for idx in range(0, len(self.points)-1):
            for j, cell_j in enumerate(self.dp_matrix[idx + 1]):
                cell_j._length = np.inf
                cell_j._next = None

                for i, cell_i in enumerate(self.dp_matrix[idx]):
                    # Compute Dubins path from (idx-1, cell_i) to (idx, cell_j)
                    logger.debug(f"Computing path between point {idx} and point {idx - 1} with angles {i} and {j}")
                    _, _, lengths = dubins_shortest_path(
                        self.points[idx][0],     self.points[idx][1],     cell_i.th(),
                        self.points[idx + 1][0], self.points[idx + 1][1], cell_j.th(),
                        self.k_max
                    )

                    curr_length = sum(lengths) + (cell_i.l() if idx > 0 and not np.isinf(cell_i.l()) else 0)

                    if curr_length < cell_j.l():
                        cell_j._length = curr_length
                        cell_j._next = self.dp_matrix[idx][i]
                    else:
                        logger.debug(f"Skipping path from point {idx + 1} angle {j} to point {idx} angle {i} with length {curr_length} because a better path with length {cell_i.l()} already exists")

                # logger.debug(f"Best from point {idx + 1} angle {j} is to point {idx} angle {i} with length {cell_j.l()}")
        logger.debug("Finished DP")


    def best_angles(self, points):
        # Backtrack to find the best angles
        best_angles = []
        for i in range(len(points)-1, -1, -1):
            if not self.dp_matrix[i]:
                continue
            best_cell = min(self.dp_matrix[i], key=lambda cell: cell.l())
            best_angles.append(best_cell.th())
            # Backtrack through the best previous cells
            while best_cell.prev():
                best_cell = best_cell.prev()
        return best_angles[::-1]
    

    def print_dp_matrix(self):
        import prettytable
        table = prettytable.PrettyTable()
        max_n_ths = max(len(row) for row in self.dp_matrix)
        table.field_names = ["Point Index"] + [f"Angle {i}" for i in range(max_n_ths)]
        for i, row in enumerate(self.dp_matrix):
            point_str = f"({self.points[i][0]:.2f}, {self.points[i][1]:.2f})"
            angles_str = []
            for cell in row:
                angle_str = f"{cell.th():.2f}" if cell.th() is not None else "None"
                length_str = f"{cell.l():.2e}" if cell.l() is not None else "None"
                id_prev = None
                if i > 0 and cell.prev() is not None:
                    for j, prev_cell in enumerate(self.dp_matrix[i - 1]):
                        if cell.prev() is not None and np.isclose(cell.prev().th(), prev_cell.th()):
                            id_prev = f"{j:4d}"
                            break

                angles_str.append(f"{angle_str}, {length_str}, {id_prev}")
            
            # Fill in empty cells if the row has fewer angles than max_n_ths
            while len(angles_str) < max_n_ths:
                angles_str.append("None, None, None")
            table.add_row([f"#{i} {point_str}"] + angles_str)

            # prev_angle = cell.prev().th() if cell.prev() else None
            # table.add_row([i, f"{cell.th():.4f}" if cell.th() is not None else "None", f"{cell.l():.4f}" if cell.l() is not None else "None", f"{prev_angle:.4f}" if prev_angle is not None else "None"])
        print(table)


    def visualize_dp_matrix(self, output_path=None, open_in_browser=True):
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
    <script>
        (function() {{
            const cells = Array.from(document.querySelectorAll('button.cell'));
            const details = document.getElementById('details');
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
            }}

            if (scaleControl) {{
                applyScale(scaleControl.value || '1');
                scaleControl.addEventListener('input', event => {{
                    applyScale(event.currentTarget.value);
                }});
            }}

            ['scroll', 'resize'].forEach(eventName => {{
                window.addEventListener(eventName, () => {{
                    if (arrowSourceCell && arrowTargetCell) {{
                        updateBackArrowPosition();
                    }}
                }}, {{ passive: true }});
            }});

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


import argparse
if __name__ == "__main__":
    args = argparse.ArgumentParser(description="Test DP class")
    args.add_argument('--debug', action='store_true', help='Enable debug logging')
    args = args.parse_args()

    if args.debug:
        logger.set_debug()

    points = [(0, 0), (1, 1), (2, 0)]
    def_thetas = [-np.pi, 0.0, np.pi]

    # points = [(0, 0), (1, 1), (2, 0)]

    dp_instance = DP(points, fixed_angles=[True, False, True], k_max=2, discretizations=90, refinements=1, def_thetas=def_thetas)
    now = time.time()
    dp_instance.solve_dp()
    logger.info(f"Solved in {time.time()-now:.4f} seconds")
    # dp_instance.print_dp_matrix()
    dp_instance.visualize_dp_matrix()


    # from dubins import plotdubins
    # dub1, _, _ = dubins_shortest_path(0, 0, -np.pi, 1, 1, 0.0, 2)
    # dub2, _, _ = dubins_shortest_path(1, 1, 0.0, 2, 0, np.pi, 2)
    # plotdubins(dub1)
    # plotdubins(dub2, color1='m', color2='c', color3='y', show=True)
