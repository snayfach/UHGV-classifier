import rich_click as click

import uhgv
from uhgv.modules import classify as classify_module
from uhgv.modules import download as download_module
from uhgv.utility import get_n_available_cpus

CONTEXT_SETTINGS = dict(help_option_names=["-h", "--help"])
click.rich_click.THEME = "modern"
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
            "name": "Input/Output options",
            "options": ["-i", "-o", "-d"],
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
    UHGV-toolkit: taxonomic classifier for the Unified Human Gut Virome Catalog

    Read the documentation at:
    https://github.com/snayfach/UHGV
    """


@cli.command(context_settings=CONTEXT_SETTINGS)
@click.argument("destination", type=click.Path(file_okay=False))
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
    Download the UHGV genome database required for the [cyan]classify[/cyan] module.

    The database will be saved in the [u]DESTINATION[/u] directory.
    """
    download_module.main(destination=destination, quiet=quiet, keep=keep)


@cli.command(context_settings=CONTEXT_SETTINGS)
@click.option(
    "-i",
    "--input",
    required=True,
    type=click.Path(exists=True),
    help="Path to input genomes in FASTA format.",
)
@click.option(
    "-o",
    "--outdir",
    required=True,
    type=click.Path(),
    help="Path to the output directory.",
)
@click.option(
    "-d",
    "--dbdir",
    required=True,
    type=click.Path(exists=True),
    help="Path to the uhgv-classify database directory.",
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
def classify(input, outdir, dbdir, sensitivity, threads, splits, continue_, quiet):
    """
    Classify new genomes into phylogenetic groups from the UHGV.
    """
    classify_module.main(
        input=input,
        outdir=outdir,
        dbdir=dbdir,
        sens=sensitivity,
        threads=threads,
        splits=splits,
        continue_=continue_,
        quiet=quiet,
    )


if __name__ == "__main__":
    cli()
