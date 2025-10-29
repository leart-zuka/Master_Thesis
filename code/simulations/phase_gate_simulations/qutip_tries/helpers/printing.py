from rich.console import Console
from rich.panel import Panel
from rich.text import Text
from rich.table import Table
from rich import box


def mini_bar(value, width=20, fill_char="█", empty_char="░"):
    """
    Return a small bar string for a normalized value in [0, inf).
    Values >1 will be capped to 1.0 for the bar (but the numeric label still shows >100%).
    """
    # protect against division by zero or NaN
    try:
        frac = float(value)
    except Exception:
        frac = 0.0
    capped = max(0.0, min(frac, 1.0))
    filled = int(round(capped * width))
    bar = f"{fill_char * filled}{empty_char * (width - filled)}"
    pct = f"{value:.7f}"
    # color the fill based on thresholds
    if value >= 1.0:
        return f"[bold red]{bar}[/bold red] {pct}"
    elif value >= 0.8:
        return f"[bold yellow]{bar}[/bold yellow] {pct}"
    else:
        return f"[green]{bar}[/green] {pct}"


# Helper for nicer formatting
def fmt_complex(z: complex) -> str:
    """Format complex numbers cleanly with color-coded sign."""
    re, im = z.real, z.imag
    sign = "+" if im >= 0 else "-"
    return f"[yellow]{re:.4f}[/yellow] [bright_black]{sign}[/bright_black] [cyan]{abs(im):.4f}i[/cyan]"


def fmt_mag(x: float) -> str:
    """Color magnitude based on value strength."""
    if x > 0.9:
        color = "green"
    elif x > 0.7:
        color = "yellow"
    else:
        color = "red"
    return f"[{color}]{x:.4f}[/{color}]"


def fmt_power(x: float) -> str:
    """Scientific notation with faint style."""
    return f"[bright_blue]{x:.3}[/bright_blue]"


def fmt_phase(x: float) -> str:
    """Color by phase sign."""
    color = "cyan" if x >= 0 else "magenta"
    return f"[{color}]{x:.3f}[/{color}]"


def print_data(title, basis, results):
    console = Console()

    table = Table(
        title=f"[bold cyan]Reflection Coefficient — Summary for {title}[/bold cyan]",
        box=box.ROUNDED,
        header_style="bold magenta",
        border_style="bright_black",
        row_styles=["none", "dim"],
    )

    table.add_column("Quantity", style="bold cyan", no_wrap=True)

    cols = list(basis.keys())
    for col_label in cols:
        pretty = col_label.replace("pi", "π").replace("V", "V")
        table.add_column(f"[bold]{pretty}[/bold]", justify="right")

    # Add rows
    table.add_row(
        "Reflection coefficient (r)", *[fmt_complex(results[c]["r"]) for c in cols]
    )
    # table.add_row("Magnitude |r|", *[fmt_mag(results[c]["|r|"]) for c in cols])
    table.add_row(
        "Reflectivity R=|r|² ", *[fmt_power(results[c]["|r|^2"]) for c in cols]
    )
    table.add_row("Phase (rad)", *[fmt_phase(results[c]["phase_rad"]) for c in cols])
    console.print(table)


def display_fidelity(fidelity: float, title: str = "CNOT Fidelity"):
    """
    Nicely display a process fidelity value with Rich and emojis.

    Args:
        fidelity (float): The process fidelity to display.
        title (str, optional): Title of the panel. Defaults to "Quantum Gate Metrics".
    """
    console = Console()
    fidelity_text = Text.assemble(
        ("✨ Process Fidelity: ", "bold magenta"),
        (f"{fidelity:.4f} ", "bold green"),
        "🎉✅",
    )
    console.print(Panel(fidelity_text, title=title, border_style="cyan"))
