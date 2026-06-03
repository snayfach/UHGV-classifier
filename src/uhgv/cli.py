import rich_click as click

import uhgv
from uhgv.modules import classify as classify_module
from uhgv.modules import download as download_module
from uhgv.utility import get_n_available_cpus

CONTEXT_SETTINGS = dict(help_option_names=["-h", "--help"])
click.rich_click.THEME = "modern"
click.rich_click.USE_RICH_MARKUP = True
click.rich_click.GROUP_ARGUMENTS_OPTIONS = True
click.rich_click.SHOW_METAVARS_COLUMN = False
click.rich_click.APPEND_METAVARS_HELP = True
click.rich_click.MAX_WIDTH = None

click.rich_click.COMMAND_GROUPS = {
    "uhgv": [
        {
            "name": "Commands",
            "commands": ["download", "classify"],
        },
    ]
}

click.rich_click.OPTION_GROUPS = {
    "uhgv download": [
        {
            "name": "Main options",
            "options": ["destination", "--keep", "--quiet"],
        },
    ],
    "uhgv classify": [
        {
            "name": "Input and output",
            "options": ["input", "output", "database"],
        },
        {
            "name": "Additional options",
            "options": ["-s", "-t", "-p", "--continue", "--quiet"],
        },
    ],
}


@click.group(context_settings=CONTEXT_SETTINGS)
@click.version_option(version=uhgv.__version__, prog_name="UHGV")
def cli():
    """
    uhgv-classifier: classification of viral genomes into UHGV taxa-like clusters.
    """


@cli.command(context_settings=CONTEXT_SETTINGS)
@click.argument(
    "destination",
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
def download(destination, keep, quiet):
    """
    Download the UHGV genome database required for the [yellow]classify[/yellow] module.
    """
    download_module.main(destination=destination, quiet=quiet, keep=keep)


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
def classify(input, output, database, sensitivity, threads, splits, continue_, quiet):
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
    )


if __name__ == "__main__":
    cli()
