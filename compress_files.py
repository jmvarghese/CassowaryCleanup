#!/usr/bin/env python3

"""
"""

import os
import subprocess
import logging
import multiprocessing
import argparse
from icecream import ic
from typing import Tuple, Set, List

def parse_command_line_args() -> str:
    """Parse command line arguments."""
    parser = argparse.ArgumentParser()
    parser.add_argument("--search_dir", "-d",
                        type=str,
                        required=True,
                        help="Directory to search for uncompressed/human readable sequence data.")
    args = parser.parse_args()
    return args.search_dir


# Function to find files
def find_files(directory: str, extensions: Tuple[str, ...]) -> List[str]:
    """ """

    file_list: List[str] = []
    
    for root, _, files in os.walk(directory):
        for file in files:
            if file.endswith(extensions) and not file.startswith("._"):
                ic(file)
                file_list.append(os.path.join(root, file))
    
    return file_list


def process_file(
    file: str,
    preprocessed_files: Set[str],
    extensions1: Tuple[str],
    extensions2: Tuple[str],
    extensions3: Tuple[str],
    log: logging.Logger,
) -> None:
    """ """

    if file not in preprocessed_files:
        if file.endswith(extensions1) or file.endswith(extensions2):
            log.info(f"Compressing: {file}")
            subprocess.run(["pigz", file], check=True)
        elif file.endswith(extensions3):
            log.info(f"Converting to BAM: {file}")
            subprocess.run(
                ["samtools", "view", "-b", "-o", f"{file}.bam", file], check=True
            )
        preprocessed_files.add(file)


def main() -> None:
    """ """

    # Set the directory to start the search
    base_directory = parse_command_line_args()

    # Define file extensions
    fasta_extensions = tuple([".fasta", ".fa"])
    fastq_extensions = tuple([".fastq", ".fq"])
    sam_extensions = tuple([".sam"])

    # Check for existing logs
    processed_files = set()
    if os.path.exists("processing_log.txt"):
        with open("processing_log.txt", "r", encoding="utf-8") as log_file:
            for line in log_file:
                processed_files.add(line.strip())

    # Logging setup
    logging.basicConfig(filename="processing_log.txt", level=logging.INFO)
    logger = logging.getLogger("file_processor")

    # Get lists of files
    with open("file_list.txt", "a", encoding="utf-8") as list_handle:
        fasta_files = find_files(base_directory, fasta_extensions)
        fastq_files = list(find_files(base_directory, fastq_extensions))
        sam_files = list(find_files(base_directory, sam_extensions))
        for item in fasta_files:
            list_handle.write(f"{item}\n")
        for item in fastq_files:
            list_handle.write(f"{item}\n")
        for item in sam_files:
            list_handle.write(f"{item}\n")

    # Create a pool of workers for parallel processing
    pool = multiprocessing.Pool(processes=multiprocessing.cpu_count())

    # Process files in parallel
    pool.starmap(
        process_file,
        [
            (
                file,
                processed_files,
                fasta_extensions,
                fastq_extensions,
                sam_extensions,
                logger,
            )
            for file in fasta_files + fastq_files + sam_files
        ],
    )

    # Close the pool
    pool.close()
    pool.join()


if __name__ == "__main__":
    main()
