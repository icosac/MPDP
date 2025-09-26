"""
Create a web-based visualization of the DP matrix. When this function is called, it will generate an HTML file. 
The HTML file should contain the table representation of the DP matrix. Each cell, should display the angle and the length until that point. 
When a cell is clicked, it should show arrows to the previous cells that lead to it, representing the transitions in the DP algorithm.
"""

def visualize_dp_matrix(dp_matrix, points):
    import webbrowser
    import os

    html_content = """
    <html>
    <head>
        <title>DP Matrix Visualization</title>
        <style>
            table { border-collapse: collapse; width: 100%; }
            th, td { border: 1px solid black; padding: 8px; text-align: center; }
            th { background-color: #f2f2f2; }
            .cell { cursor: pointer; }
            .arrow { color: red; font-weight: bold; }
        </style>
        <script>
            function showArrows(cellId) {
                var arrows = document.getElementsByClassName('arrow');
                for (var i = 0; i < arrows.length; i++) {
                    arrows[i].style.display = 'none';
                }
                var cell = document.getElementById(cellId);
                var prevCells = cell.getAttribute('data-prev').split(',');
                for (var i = 0; i < prevCells.length; i++) {
                    if (prevCells[i]) {
                        document.getElementById(prevCells[i]).style.display = 'inline';
                    }
                }
            }
            function hideArrows() {
                var arrows = document.getElementsByClassName('arrow');
                for (var i = 0; i < arrows.length; i++) {
                    arrows[i].style.display = 'none';
                }
            }
            document.addEventListener('click', function(event) {
                if (!event.target.classList.contains('cell')) {
                    hideArrows();
                }
            });
        </script>
    </head>
    <body>
        <h1>DP Matrix Visualization</h1>
        <table>
            <tr>
                <th>Point Index</th>
                <th>Angle (rad)</th>
                <th>Length</th>
                <th>Previous Cells</th>
            </tr>
    """

    for idx, row in enumerate(dp_matrix):
        for jdx, cell in enumerate(row):
            cell_id = f"cell_{idx}_{jdx}"
            prev_ids = []
            prev_cell = cell.prev()
            while prev_cell:
                prev_idx = dp_matrix.index(prev_cell)
                prev_jdx = dp_matrix[prev_idx].index(prev_cell)
                prev_ids.append(f"cell_{prev_idx}_{prev_jdx}")
                prev_cell = prev_cell.prev()
            prev_ids_str = ','.join(prev_ids)
            html_content += f"""
            <tr>
                <td>{idx}</td>
                <td class="cell" id="{cell_id}" data-prev="{prev_ids_str}" onclick="showArrows('{cell_id}')">{cell.th():.2f}</td>
                <td>{cell.l():.2f}</td>
                <td>{', '.join(prev_ids)}</td>
            </tr>
            """
            for prev_id in prev_ids:
                html_content += f'<div class="arrow" id="{prev_id}" style="display:none;">&#8592;</div>'
    html_content += """
        </table>
    </body>
    </html>
    """
    file_path = os.path.abspath("dp_matrix_visualization.html")
    with open(file_path, "w") as f:
        f.write(html_content)

    webbrowser.open(f"file://{file_path}")
    print(f"DP matrix visualization saved to {file_path} and opened in web browser.")
