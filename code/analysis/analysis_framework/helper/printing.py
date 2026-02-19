from rich.console import Console
from rich.table import Table
from rich.panel import Panel
from rich import box
import numpy as np

console = Console()


def pretty_print_gate(matrix, title, fidelity_indices):
    """
    matrix: 4x4 numpy array
    fidelity_indices: list of (i, j) tuples used for fidelity calculation
    """

    matrix = np.array(matrix)

    # ---- Create Table ----
    table = Table(
        title=title,
        box=box.SIMPLE_HEAVY,
        show_lines=True,
        header_style="bold magenta",
    )

    basis = ["00", "01", "10", "11"]

    table.add_column("⟨out| in⟩", justify="center", style="bold")
    for b in basis:
        table.add_column(b, justify="right")

    # ---- Fill Table ----
    for i, row_label in enumerate(basis):
        row = [row_label]

        for j in range(4):
            value = matrix[i, j]

            # Highlight fidelity-relevant entries
            if (i, j) in fidelity_indices:
                row.append(f"[bold green]{value: .4f}[/bold green]")
            else:
                row.append(f"{value: .4f}")

        table.add_row(*row)

    # ---- Fidelity ----
    fidelity = sum(matrix[i, j] for i, j in fidelity_indices) / 4

    console.print(table)
    console.print(
        Panel.fit(
            f"[bold cyan]Fidelity = {fidelity:.6f}[/bold cyan]",
            border_style="cyan",
        )
    )
