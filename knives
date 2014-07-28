#!/usr/bin/env python

from __future__ import print_function
from csv import DictReader
import datetime
from distutils.spawn import find_executable
from docopt import docopt
import gzip
import logging
import os
from os import path
import subprocess as sp
import tempfile

CLI_DOC = """
USAGE:
    knives.py [-K SICKLE -C SCYTHE -E SEQQS -P PAIRS -q QUAL -t TYPE -l LEN -s SKIP \
            -L LOGDIR -S -n ] -a CONTAM -i IN1 -I IN2 -p PREFIX -r PRIOR

OPTIONS:
    -K SICKLE   Location of the sickle executable.
    -C SCYTHE   Location of the scythe executable.
    -P PAIRS    Location of the pairs executable.
    -E SEQQS    Location of the seqqs executable.
    -q QUAL     Minimum quality score for windowed quality trimming
                [default: 20]
    -t TYPE     Quality score encoding. One of phred, sanger, solexa or
                illumina. See the help message of scythe for an explanation.
                [default: sanger]
    -l LEN      Minimum length of sequence remaining after trimming.
                [default: 64]
    -s SKIP     Skip steps. Comma seperated list of program names. Must be one
                or more of: scythe, sickle, seqqs.
    -S          Don't log stderr of programs.
    -n          Remove any sequence with Ns.
    -L LOGDIR   Directory to store logs. [Default: .]
    -a CONTAM   Fasta file containing adaptor sequences, for scythe.
    -i IN1      Input R1 file.
    -I IN2      Input R2 file.
    -p PREFIX   Prefix of all ouput files, including seqqs reports and log
                files.
    -r PRIOR    Prior probabilty of adaptor contamination. See scythe docs.
"""

# declare global vars
SICKLE = None
SCYTHE = None
SEQQS = None
PAIRS = None
CMD_LOG = None
STDERR_LOG = None
NOW = datetime.datetime.now().strftime("%Y-%m-%d_%H:%M:%S")

def setup_logs(opts):
    """Setup logging and store loggers in global vars."""
    global CMD_LOG
    global STDERR_LOG
    if not path.exists(opts["-L"]):
        raise RuntimeError("Log directory not found.")
    # Logger for the exact comands run
    CMD_LOG = logging.getLogger("knives_cmds")
    cmd_lfn = path.join(opts["-L"], opts["-p"] + "_cmds_" + NOW + ".log")
    cmd_fh = logging.FileHandler(cmd_lfn)
    cmd_fh.setLevel(logging.INFO)
    cmd_fmt = logging.Formatter('%(asctime)s:\t%(message)s',
                                "%Y-%m-%d %H:%M:%S")
    cmd_fh.setFormatter(cmd_fmt)
    CMD_LOG.addHandler(cmd_fh)
    CMD_LOG.setLevel(logging.INFO)
    # Logger for the stderr of each process
    STDERR_LOG = logging.getLogger("knives_stderr")
    if not opts["-S"]:
        stderr_lfn = path.join(opts["-L"], opts["-p"] + "_stderrs_" + NOW + ".log")
        stderr_fh = logging.FileHandler(stderr_lfn)
        stderr_fh.setLevel(logging.DEBUG)
        stderr_fmt = logging.Formatter(
            '%(asctime)20s - %(levelname)8s: %(message)s',
            "%Y-%m-%d %H:%M:%S")
        stderr_fh.setFormatter(stderr_fmt)
        STDERR_LOG.addHandler(stderr_fh)
    else:
        STDERR_LOG.addHandler(logging.NullHandler())
    STDERR_LOG.setLevel(logging.INFO)

def find_progs(opts):
    """Finds EXEs of the 4 programs needed; stores paths in global vars."""
    global SICKLE
    global SCYTHE
    global SEQQS
    global PAIRS
    # Find sickle
    if opts["-K"]:
        SICKLE = opts["-K"]
    else:
        SICKLE = find_executable("sickle")
    if not SICKLE or \
            not (os.path.isfile(SICKLE) and os.access(SICKLE, os.X_OK)):
        raise RuntimeError("Can't find sickle")
    # Find Scythe
    if opts["-C"]:
        SCYTHE = opts["-C"]
    else:
        SCYTHE = find_executable("scythe")
    if not SCYTHE or \
            not (os.path.isfile(SCYTHE) and os.access(SCYTHE, os.X_OK)):
        raise RuntimeError("Can't find scythe")
    # Find seqqs
    if opts["-E"]:
        SEQQS = opts["-E"]
    else:
        SEQQS = find_executable("seqqs")
    if not SEQQS or \
            not (os.path.isfile(SEQQS) and os.access(SEQQS, os.X_OK)):
        raise RuntimeError("Can't find seqqs")
    # Find pairs
    if opts["-P"]:
        PAIRS = opts["-P"]
    else:
        PAIRS = find_executable("pairs")
    if not PAIRS or \
            not (os.path.isfile(PAIRS) and os.access(PAIRS, os.X_OK)):
        raise RuntimeError("Can't find sickle")

def main(opts):
    global SICKLE
    global SCYTHE
    global SEQQS
    global PAIRS
    global CMD_LOG
    global STDERR_LOG
    if not path.exists(opts['-i']) or  not path.exists(opts['-I']):
        print("Input files must exist!")
        exit(1)
    if not path.exists(opts['-a']):
        print("Adapter file must exist!")
        exit(1)
    if os.sep in opts['-p']:
        if not path.isdir(path.dirname(opts['-p'])):
            print("Prefix is invalid path.")
            print("dirname(p) must exist if path contains {}.".format(os.sep))
            exit(1)
    setup_logs(opts)
    find_progs(opts)
    # interleave reads
    interleave_cli = [
        PAIRS,  # EXE path
        "join",  # join command to interleave 2 fqs
        "-t",  # write /1 /2 tags to header
        opts['-i'],  # R1 file
        opts['-I'],  # R2 file
    ]
    interleave_stderr = tempfile.TemporaryFile("rw+b")
    CMD_LOG.info(" ".join(interleave_cli))
    interleave_proc = sp.Popen(interleave_cli, stdout=sp.PIPE,
                               stderr=interleave_stderr)
    # initial QC
    seqqs_il_cli = [
        SEQQS,  # exe path
        "-i",  # interleaved input
        "-e",  # enable streaming
        "-q",  # Qual encoding is:
        opts['-t'],
        "-p",  # prefix is:
        opts["-p"] + "_initial_",
        "-",  # Use stdin
    ]
    CMD_LOG.info(" ".join(seqqs_il_cli))
    seqqs_il_stderr = tempfile.TemporaryFile("rw+b")
    seqqs_il_proc = sp.Popen(seqqs_il_cli, stdin=interleave_proc.stdout,
                             stdout=sp.PIPE, stderr=seqqs_il_stderr)
    interleave_proc.stdout.close()
    # scythe
    scythe_cli = [
        SCYTHE,  # exe path
        "-p",  # prior is:
        opts["-r"],
        "-a",  # Adaptor file is:
        opts["-a"],
        "-M", opts['-l'],
        "/dev/stdin",  # Use stdin
    ]
    CMD_LOG.info(" ".join(scythe_cli))
    scythe_stderr = tempfile.TemporaryFile("rw+b")
    scythe_proc = sp.Popen(scythe_cli, stdin=seqqs_il_proc.stdout,
            stdout=sp.PIPE, stderr=scythe_stderr)
    seqqs_il_proc.stdout.close()
    # QC after adaptor removal
    seqqs_na_cli = [
        SEQQS,  # exe path
        "-i",  # interleaved input
        "-e",  # enable streaming
        "-q",  # Qual encoding is:
        opts['-t'],
        "-p",  # prefix is:
        opts["-p"] + "_noadapt",
        "-",  # Use stdin
    ]
    CMD_LOG.info(" ".join(seqqs_na_cli))
    seqqs_na_stderr = tempfile.TemporaryFile("rw+b")
    seqqs_na_proc = sp.Popen(seqqs_na_cli, stdin=scythe_proc.stdout,
                             stdout=sp.PIPE, stderr=seqqs_na_stderr)
    scythe_proc.stdout.close()
    # Sickle
    sickle_cli = [
        SICKLE,  # exe path
        "pe",  # paired end mode
        "-t",  # Qual encoding is:
        opts['-t'],
        "-c", "/dev/stdin",  # Use stdin for input
        "-M", "/dev/stdout",  # And stdout for output
        "-q",  # min window quality is:
        opts["-q"],
        "-l",  # min length is:
        opts["-l"],
    ]
    if opts['-n']:
        sickle_cli.append("-n")  # Remove any read w/ Ns
    CMD_LOG.info(" ".join(sickle_cli))
    sickle_stderr = tempfile.TemporaryFile("rw+b")
    sickle_proc = sp.Popen(sickle_cli, stdin=seqqs_na_proc.stdout,
                           stdout=sp.PIPE, stderr=sickle_stderr)
    seqqs_il_proc.stdout.close()
    # QC after qual trimming
    seqqs_qt_cli = [
        SEQQS,  # exe path
        "-i",  # interleaved input
        "-e",  # enable streaming
        "-q",  # Qual encoding is:
        opts['-t'],
        "-p",  # prefix is:
        opts["-p"] + "_qualtrim",
        "-",  # Use stdin
    ]
    CMD_LOG.info(" ".join(seqqs_qt_cli))
    seqqs_qt_stderr = tempfile.TemporaryFile("rw+b")
    seqqs_qt_proc = sp.Popen(seqqs_qt_cli, stdin=sickle_proc.stdout,
                             stdout=sp.PIPE, stderr=seqqs_qt_stderr)
    sickle_proc.stdout.close()
    # deinterleave
    deinterleave_cli = [
        PAIRS,  # exe path
        "split",  # interleaved input
        "-1", opts['-p'] + "_R1.knives.fastq",
        "-2", opts['-p'] + "_R2.knives.fastq",
        "-u", opts['-p'] + "_singles.knives.fastq",
        "-",  # Use stdin
    ]
    CMD_LOG.info(" ".join(deinterleave_cli))
    deinterleave_stderr = tempfile.TemporaryFile("rw+b")
    deinterleave_stdout = tempfile.TemporaryFile("rw+b")
    deinterleave_proc = sp.Popen(deinterleave_cli, stdin=seqqs_qt_proc.stdout,
                                 stdout=deinterleave_stdout,
                                 stderr=deinterleave_stderr)
    sickle_proc.stdout.close()
    deinterleave_proc.wait()
    # reset al stderr files
    interleave_stderr.seek(0)
    seqqs_il_stderr.seek(0)
    scythe_stderr.seek(0)
    seqqs_na_stderr.seek(0)
    sickle_stderr.seek(0)
    seqqs_qt_stderr.seek(0)
    deinterleave_stderr.seek(0)
    STDERR_LOG.info(interleave_stderr.read())
    STDERR_LOG.info(seqqs_il_stderr.read())
    STDERR_LOG.info(scythe_stderr.read())
    STDERR_LOG.info(seqqs_na_stderr.read())
    STDERR_LOG.info(sickle_stderr.read())
    STDERR_LOG.info(seqqs_qt_stderr.read())
    STDERR_LOG.info(deinterleave_stderr.read())

if __name__ == "__main__":
    opts = docopt(CLI_DOC)
    main(opts)
