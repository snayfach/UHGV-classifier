from typing import Any

import rich_click as click
from rich.containers import Renderables
from rich.rule import Rule
from rich.text import Text
from rich_click.rich_command import RichCommand, RichGroup
from rich_click.rich_context import RichContext
from rich_click.rich_help_formatter import RichHelpFormatter
from rich_click.rich_panel import RichOptionPanel

import uhgv
from uhgv.subcommands import classify as classify_module
from uhgv.subcommands import download_database as download_module
from uhgv.utility import get_n_available_cpus

UHGV_HELP_CONFIG: dict[str, Any] = {
    "theme": "modern",
    "text_markup": "rich",
    "max_width": None,
    "padding_usage": (0, 0, 1, 0),
    "padding_helptext": (0, 0, 1, 0),
    "style_command_help": "italic",
    "options_table_column_types": ["opt_short", "opt_long", "help"],
    "options_table_help_sections": [
        "required",
        "help",
        "envvar",
        "default",
        "deprecated",
        "metavar",
    ],
    "command_groups": {
        "uhgv": [
            {
                "name": "Commands",
                "commands": ["download-database", "classify"],
            },
        ]
    },
}

uhgv_rich_config = click.rich_config(UHGV_HELP_CONFIG)
CONTEXT_SETTINGS = dict(help_option_names=["-h", "--help"])


class FlushPanel:
    """Panel-like renderable without Rich's invisible border gutter."""

    def __init__(
        self,
        renderable: Any,
        *,
        title: Any = None,
        border_style: Any = "",
        **_: Any,
    ):
        self.renderable = renderable
        self.title = title
        self.border_style = border_style

    def __rich_console__(self, console: Any, options: Any):
        renderables: list[Any] = []
        if self.title is not None:
            renderables.append(self.title)
        renderables.extend(
            [
                Rule(style=self.border_style or ""),
                self.renderable,
                Text(""),
            ]
        )
        yield from console.render(Renderables(renderables), options)


class FlushOptionPanel(RichOptionPanel):
    """Render Arguments and Options flush with the left edge of the help output."""

    panel_class = FlushPanel

    def get_table(self, command: Any, ctx: Any, formatter: RichHelpFormatter) -> Any:
        if self.name == formatter.config.arguments_panel_title:
            help_sections = formatter.config.options_table_help_sections
            formatter.config.options_table_help_sections = [
                section for section in help_sections if section != "required"
            ]
            try:
                table = super().get_table(command, ctx, formatter)
            finally:
                formatter.config.options_table_help_sections = help_sections
        else:
            table = super().get_table(command, ctx, formatter)

        table.show_edge = False
        return table


class UHGVHelpFormatter(RichHelpFormatter):
    option_panel_class = FlushOptionPanel


class UHGVRichContext(RichContext):
    formatter_class = UHGVHelpFormatter


class UHGVRichCommand(RichCommand):
    context_class = UHGVRichContext


class UHGVRichGroup(RichGroup):
    context_class = UHGVRichContext
    command_class = UHGVRichCommand


@click.group(cls=UHGVRichGroup, context_settings=CONTEXT_SETTINGS)
@click.version_option(version=uhgv.__version__, prog_name="UHGV")
@uhgv_rich_config
def cli():
    """
    uhgv-classifier: classification of viral genomes into UHGV taxa-like clusters.
    """


@cli.command(context_settings=CONTEXT_SETTINGS)
@click.argument(
    "output",
    type=click.Path(file_okay=False),
    help="Directory to write the uhgv-classify database.",
)
@click.option(
    "--keep",
    is_flag=True,
    default=False,
    show_default=True,
    help="Do not delete the compressed database file.",
)
@click.option(
    "--quiet",
    is_flag=True,
    default=False,
    show_default=True,
    help="Suppress logging messages.",
)
def download_database(output, keep, quiet):
    """
    Download the reference database required for using the [cyan]classify[/cyan]
    subcommand.
    """
    download_module.main(output=output, quiet=quiet, keep=keep)


@cli.command(context_settings=CONTEXT_SETTINGS)
@click.argument(
    "input",
    type=click.Path(exists=True),
    help="Input FASTA file containing the genomes to classify.",
)
@click.argument(
    "output",
    type=click.Path(),
    help="Output directory to write the classification results.",
)
@click.argument(
    "database",
    type=click.Path(exists=True),
    help="Path to the uhgv-classify database.",
)
@click.option(
    "-s",
    "--sensitivity",
    type=click.Choice(["fast", "sensitive", "very-sensitive"], case_sensitive=False),
    default="sensitive",
    show_default=True,
    help="DIAMOND search sensitivity.",
)
@click.option(
    "-t",
    "--threads",
    type=int,
    default=get_n_available_cpus(),
    show_default=True,
    help="Number of threads to use.",
)
@click.option(
    "-p",
    "--splits",
    type=int,
    default=None,
    help="Number of BLASTN jobs to spawn in parallel. "
    "Defaults to the number of threads.",
)
@click.option(
    "--cleanup",
    is_flag=True,
    default=False,
    show_default=True,
    help="Remove the tmp temporary directory after the pipeline finishes.",
)
@click.option(
    "--continue",
    "continue_",
    is_flag=True,
    default=False,
    show_default=True,
    help="Continue where the program left off.",
)
@click.option(
    "--quiet",
    is_flag=True,
    default=False,
    show_default=True,
    help="Suppress logging messages.",
)
def classify(
    input, output, database, sensitivity, threads, splits, continue_, quiet, cleanup
):
    """
    Classify new genomes into UHGV taxa-like clusters.
    """
    classify_module.main(
        input=input,
        outdir=output,
        dbdir=database,
        sens=sensitivity,
        threads=threads,
        splits=splits,
        continue_=continue_,
        quiet=quiet,
        cleanup=cleanup,
    )


if __name__ == "__main__":
    cli()
