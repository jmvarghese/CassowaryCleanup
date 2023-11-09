# Resumable Deep Storage Converter for Sequence Data Files

The module `compress_files.py` uses parallel subprocesses to compress FASTQ and FASTA format files with `pigz`. It also uses `samtools view` to convert human-readable SAM format files to binary BAM format files. Be sure to have `pigz` and `samtools` installed and on your command line's `$PATH` before running the module. Users should consult the documentation in `docs/` for more information about the module's API. 

To use the module, simply specify the directory where the module should search and watch it go:

```
python3 compress_files.py --search_dir /path/to/directory
```

The module will then search through the provided directory and any of its subdirectories for files that could be compressed or converted to save disk space. The module records the files it processes as it goes, which makes it possible for it to be resumed in case it is interrupted in the middle of a very long queue of files.

Users should also be aware that the module has a dry run option available via the flag `--dry-run`/`-d`. When this option is set to true, the script will do everything except the compression/conversion. This option can be used to assess how many files are available to be compressed or converted via the log file.
