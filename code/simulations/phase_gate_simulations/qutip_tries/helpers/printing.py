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
def fmt_complex(z: complex, err=None) -> str:
    """
    Format complex numbers cleanly with color-coded sign.
    Optionally include an error (err) which may be a real scalar or a complex.
    Examples:
      - fmt_complex(1+0.2j) -> "1.0000 + 0.2000i"
      - fmt_complex(1+0.2j, 0.01) -> "1.0000 + 0.2000i ± 0.0100"
      - fmt_complex(1+0.2j, 0.005+0.002j) -> "1.0000 + 0.2000i ± (0.0050 + 0.0020i)"
    """
    if z is None:
        return "-"
    re, im = z.real, z.imag
    sign = "+" if im >= 0 else "-"
    main = f"[yellow]{re:.4f}[/yellow] [bright_black]{sign}[/bright_black] [cyan]{abs(im):.4f}i[/cyan]"

    if err is None:
        return main

    # Try to detect complex-like error (has real and imag) first
    try:
        if isinstance(err, complex) or (hasattr(err, "real") and hasattr(err, "imag")):
            e = complex(err)
            e_re, e_im = e.real, e.imag
            e_sign = "+" if e_im >= 0 else "-"
            err_str = f"[yellow]{e_re:.4f}[/yellow] [bright_black]{e_sign}[/bright_black] [cyan]{abs(e_im):.4f}i[/cyan]"
            return f"{main} ± ({err_str})"
        else:
            # treat as scalar
            e_val = float(err)
            return f"{main} ± [bright_black]{e_val:.4f}[/bright_black]"
    except Exception:
        # fallback: just str() the error
        return f"{main} ± {err}"


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


def fmt_phase(x: float, err=None) -> str:
    """Color by phase sign."""
    color = "cyan" if x >= 0 else "magenta"
    if err is None:
        return f"[{color}]{x:.3f}[/{color}]"
    return f"[{color}]{x:.3f}[/{color}] ± [{color}]{err:.3f}[/{color}]"


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
        "Reflection coefficient (r)",
        *[fmt_complex(results[c]["r"], results[c]["dr"]) for c in cols],
    )
    table.add_row(
        "Reflectivity R=|r|² ", *[fmt_power(results[c]["|r|^2"]) for c in cols]
    )
    table.add_row(
        "Phase (rad)",
        *[
            fmt_phase(results[c]["phase_rad"], results[c]["phase_rad_err"])
            for c in cols
        ],
    )
    console.print(table)


def display_fidelity(
    fidelity: float,
    error: float | None = None,
    title: str = "CNOT Fidelity",
):
    """
    Nicely display a process fidelity value with Rich and emojis in the format:
    value ± error

    Args:
        fidelity (float): The process fidelity value.
        error (float, optional): Statistical error for the fidelity.
        title (str, optional): Title for the displayed panel.
    """
    console = Console()

    if error is not None:
        text_str = f"{fidelity:.4f} ± {error:.4f}"
        emoji = "📏✨"
    else:
        text_str = f"{fidelity:.4f}"
        emoji = "🎉✅"

    fidelity_text = Text.assemble(
        ("Process Fidelity: ", "bold magenta"), (text_str + " ", "bold green"), emoji
    )

    console.print(Panel(fidelity_text, title=title, border_style="cyan"))
