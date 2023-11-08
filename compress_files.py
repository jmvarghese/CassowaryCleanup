#!/usr/bin/env python3

"""
"""

import os
import subprocess
import multiprocessing
import argparse
from typing import Tuple, Set, List
from filelock import FileLock, Timeout
import time


def parse_command_line_args() -> Tuple[str, str]:
    """Parse command line arguments."""
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--search_dir",
        "-d",
        type=str,
        required=True,
        help="Directory to search for uncompressed/human readable sequence data.",
    )
    parser.add_argument(
        "--logging_file",
        "-l",
        type=str,
        required=False,
        default="compression.log",
        help="Directory to search for uncompressed/human readable sequence data.",
    )
    args = parser.parse_args()
    return args.search_dir, args.logging_file


# Function to find files
def find_files(directory: str, extensions: Tuple[str, ...]) -> List[str]:
    """ """

    file_list: List[str] = []

    for root, _, files in os.walk(directory):
        for file in files:
            if file.endswith(extensions) and not file.startswith("._"):
                # ic(file)
                file_list.append(os.path.join(root, file))

    return file_list


def process_file(
    file: str,
    preprocessed_files: Set[str],
    extensions1: Tuple[str],
    extensions2: Tuple[str],
    extensions3: Tuple[str],
    log_path: str,
) -> None:
    """ """
    t1 = time.time()
    orig_size = os.stat(file).st_size / (1024*1024)
    
    if file in preprocessed_files:
        return

    if file.endswith(extensions1) or file.endswith(extensions2):
        subprocess.run(["pigz", file], check=True)
    elif file.endswith(extensions3):
        subprocess.run(
            ["samtools", "view", "-b", "-o", f"{file}.bam", file], check=True
        )
    
    time_elapsed = time.time() - t1
    fin_size = os.stat(f"{file}.gz").st_size / (1024*1024)

    lock_path = f"{log_path}.lock"
    lock = FileLock(lock_path, timeout=60)
    lock.acquire()
    try:
        with open(log_path, "a", encoding="utf-8") as log:
            log.write(f"{file}\t{time_elapsed}\t{orig_size}\t{fin_size}\n")
    finally:
        lock.release()


def main() -> None:
    """ """

    # Set the directory to start the search
    base_directory, log_path = parse_command_line_args()

    # Define file extensions
    fasta_extensions = tuple([".fasta", ".fa"])
    fastq_extensions = tuple([".fastq", ".fq"])
    sam_extensions = tuple([".sam"])

    # Check for existing logs
    processed_files = set()
    if os.path.exists("compression.log"):
        with open("compression.log", "r", encoding="utf-8") as log_file:
            for line in log_file:
                processed_files.add(line.strip())
    else:
        with open("compression.log", "w", encoding="utf-8") as log:
            log.write("File\tTime (ms)\tOriginal Size\tCompressed Size\n")
        
    # Get lists of files
    fasta_files = find_files(base_directory, fasta_extensions)
    fastq_files = find_files(base_directory, fastq_extensions)
    sam_files = find_files(base_directory, sam_extensions)

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
                log_path,
            )
            for file in fasta_files + fastq_files + sam_files
        ],
    )

    # Close the pool
    pool.close()
    pool.join()


if __name__ == "__main__":
    main()
