import os
import subprocess
import logging
import multiprocessing

# Function to find files
def find_files(directory, extensions):
	for root, _, files in os.walk(directory):
		for file in files:
			if file.endswith(tuple(extensions)):
				yield os.path.join(root, file)

def process_file(file, processed_files):
	if file not in processed_files:
		if file.endswith(fasta_extensions) or file.endswith(fastq_extensions):
			logger.info(f"Compressing: {file}")
			subprocess.run(["pigz", file])
		elif file.endswith(sam_extensions):
			logger.info(f"Converting to BAM: {file}")
			subprocess.run(["samtools", "view", "-b", "-o", f"{file}.bam", file])
		processed_files.add(file)

if __name__ == '__main':
	# Set the directory to start the search
	base_directory = "/Users/jmvarghese/Documents/Compress Files"
	# Define file extensions
	fasta_extensions = (".fasta", ".fa")
	fastq_extensions = (".fastq", ".fq")
	sam_extensions = (".sam")

	# Check for existing logs
	processed_files = set()
	if os.path.exists("processing_log.txt"):
		with open("processing_log.txt", "r") as log_file:
			for line in log_file:
				processed_files.add(line.strip())

	# Logging setup
	logging.basicConfig(filename="processing_log.txt", level=logging.INFO)
	logger = logging.getLogger("file_processor")

	# Get lists of files
	f = open("file_list.txt", "w")
	fasta_files = list(find_files(base_directory, fasta_extensions))
	fastq_files = list(find_files(base_directory, fastq_extensions))
	sam_files = list(find_files(base_directory, sam_extensions))
	f.write(fasta_files. fastq_files, sam_files)
	f.close()
	
	# Create a pool of workers for parallel processing
	pool = multiprocessing.Pool(processes=multiprocessing.cpu_count())

	# Process files in parallel
	pool.starmap(process_file, [(file, processed_files) for file in fasta_files + fastq_files + sam_files])

	# Close the pool
	pool.close()
	pool.join()