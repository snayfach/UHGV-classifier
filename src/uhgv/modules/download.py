import os
import shutil
import subprocess as sp
import sys
import time
import urllib.request

from uhgv import utility


class DatabaseDownloader:
    def __init__(self, destination):
        self.url = "https://portal.nersc.gov/UHGV/toolkit/"
        self.destination = destination
        self.version = (
            urllib.request.urlopen(self.url + "CURRENT_RELEASE.txt")
            .read()
            .decode("utf-8")
            .strip()
        )
        self.filename = self.version + ".tar.gz"
        self.output_file = os.path.join(self.destination, self.filename)

    def download(self):
        database_url = self.url + self.filename
        with urllib.request.urlopen(database_url) as response:
            with open(self.output_file, "wb") as fout:
                shutil.copyfileobj(response, fout)

    def extract(self):
        shutil.unpack_archive(self.output_file, self.destination, "gztar")
        os.remove(self.output_file)

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


def main(destination, quiet=False):
    program_start = time.time()
    console = utility.ConsoleLogger(quiet)
    if not os.path.exists(destination):
        os.makedirs(destination)

    console.log("UHGV download")

    utility.check_executables(["makeblastdb", "diamond"])

    console.log("Checking latest version of database...")
    db = DatabaseDownloader(destination)

    with console.status("Downloading..."):
        db.download()

    with console.status("Extracting..."):
        db.extract()

    console.log("Building BLASTN database...")
    db.blastn_makedb()

    console.log("Building DIAMOND database...")
    db.diamond_makedb()

    console.log("Run time: %s seconds" % round(time.time() - program_start, 2))
    console.log("Peak mem: %s GB" % round(utility.max_mem_usage(), 2))
