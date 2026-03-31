from rich.console import Console
from rich.table import Table
from rich.panel import Panel
from rich import box
import numpy as np

console = Console()


def pm(val, err, fmt="{:.2f}"):
    # returns "val ± err" with consistent formatting
    return f"{fmt.format(val)} ± {fmt.format(err)}"


def pretty_print_gate(
    matrix,
    title,
    fidelity_indices,
    atom_basis,
    photon_basis,
):
    """
    matrix: 4x4 array
    fidelity_indices: list of (i, j) tuples
    atom_basis: ["|0⟩", "|1⟩"] (or custom)
    photon_basis: ["H", "V"] (or ["R", "L"], etc.)
    """

    matrix = np.array(matrix)

    # Construct tensor-product basis labels
    joint_basis = [f"{a} ⊗ {p}" for a in atom_basis for p in photon_basis]

    # ---- Create Table ----
    table = Table(
        title=title,
        box=box.SIMPLE_HEAVY,
        show_lines=True,
        header_style="bold magenta",
    )

    table.add_column("⟨out| in⟩", justify="center", style="bold")

    for label in joint_basis:
        table.add_column(label, justify="right")

    # ---- Fill Table ----
    for i, row_label in enumerate(joint_basis):
        row = [row_label]

        for j in range(4):
            value = matrix[i, j]

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
