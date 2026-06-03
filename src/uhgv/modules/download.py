import os
import shutil
import subprocess as sp
import sys
import time
from functools import partial
from urllib.request import urlopen

from rich.progress import (
    BarColumn,
    DownloadColumn,
    Progress,
    TextColumn,
    TimeRemainingColumn,
    TransferSpeedColumn,
)

from uhgv import utility


class DatabaseDownloader:
    def __init__(self, destination, console):
        self.url = "https://portal.nersc.gov/UHGV/toolkit/"
        self.destination = destination
        self.console = console
        self.version = (
            urlopen(self.url + "CURRENT_RELEASE.txt").read().decode("utf-8").strip()
        )
        self.filename = self.version + ".tar.gz"
        self.database_url = self.url + self.filename
        self.output_file = os.path.join(self.destination, self.filename)

    def _copy_url(self, task_id, progress):
        progress.console.log(
            f"Requesting [blue link={self.database_url}]{self.database_url}[/blue link]."
        )
        response = urlopen(self.database_url)
        progress.update(task_id, total=int(response.info()["Content-length"]))
        with open(self.output_file, "wb") as dest_file:
            progress.start_task(task_id)
            for data in iter(partial(response.read, 32768), b""):
                dest_file.write(data)
                progress.update(task_id, advance=len(data))

    def download(self):
        progress = Progress(
            TextColumn("{task.fields[filename]}", justify="right", style="green"),
            BarColumn(bar_width=None),
            "[progress.percentage]{task.percentage:>3.1f}%",
            "|",
            DownloadColumn(),
            "|",
            TransferSpeedColumn(),
            "|",
            TimeRemainingColumn(elapsed_when_finished=True),
            console=self.console.console,
            transient=True,
        )
        with progress:
            task_id = progress.add_task("download", filename=self.filename, start=False)
            self._copy_url(task_id, progress)

    def extract(self):
        shutil.unpack_archive(self.output_file, self.destination, "gztar")

    def blastn_makedb(self):
        self.dbdir = os.path.join(
            self.destination, self.filename.replace(".tar.gz", "")
        )
        logpath = f"{self.dbdir}/genomes.log"
        cmd = [
            "makeblastdb",
            "-in",
            f"{self.dbdir}/genomes.fna",
            "-out",
            f"{self.dbdir}/genomes",
            "-dbtype",
            "nucl",
        ]
        with open(logpath, "w") as log:
            result = sp.run(cmd, stdout=sp.DEVNULL, stderr=log)
        if result.returncode != 0:
            msg = "\nError: BLASTN database failed to build\n"
            msg += f"See log for details: {logpath}"
            sys.exit(msg)

    def diamond_makedb(self):
        self.dbdir = os.path.join(
            self.destination, self.filename.replace(".tar.gz", "")
        )
        logpath = f"{self.dbdir}/proteins.log"
        cmd = [
            "diamond",
            "makedb",
            "--in",
            f"{self.dbdir}/proteins.faa",
            "--db",
            f"{self.dbdir}/proteins",
        ]
        with open(logpath, "w") as log:
            result = sp.run(cmd, stdout=sp.DEVNULL, stderr=log)
        if result.returncode != 0:
            msg = "\nError: DIAMOND database failed to build\n"
            msg += f"See log for details: {logpath}"
            sys.exit(msg)


def main(destination, quiet=False, keep=False):
    program_start = time.time()
    console = utility.ConsoleLogger(quiet)
    if not os.path.exists(destination):
        os.makedirs(destination)

    console.log("UHGV download")

    utility.check_executables(["makeblastdb", "diamond"])

    console.log("Checking latest version of database...")
    db = DatabaseDownloader(destination, console)

    db.download()

    with console.status("Extracting..."):
        db.extract()

    console.log("Building BLASTN database...")
    db.blastn_makedb()

    console.log("Building DIAMOND database...")
    db.diamond_makedb()

    if not keep:
        os.remove(db.output_file)

    console.log("Run time: %s seconds" % round(time.time() - program_start, 2))
    console.log("Peak mem: %s GB" % round(utility.max_mem_usage(), 2))
