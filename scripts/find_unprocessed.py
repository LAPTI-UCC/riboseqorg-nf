import csv
import subprocess
import os

def get_processed_runs():
    # Get the list of processed sqlite files using the shell command
    cmd = "ls data/Processed/sqlites/*.sqlite | xargs -n 1 basename | sed 's/\.sqlite$//'"
    process = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    output, _ = process.communicate()
    
    # Convert bytes to string and split into list
    processed = output.decode().strip().split('\n')
    return set(processed)

def get_unprocessed_runs(csv_path, processed_runs):
    unprocessed = []
    
    # Read the CSV file
    with open(csv_path, 'r') as f:
        csv_reader = csv.reader(f)
        next(csv_reader)  # Skip header if present
        for row in csv_reader:
            run_id = row[0]  # Assuming run ID is in first column
            if run_id not in processed_runs:
                unprocessed.append(run_id)
    
    return unprocessed

def main():
    # Get processed runs
    processed_runs = get_processed_runs()
    
    # Get unprocessed runs
    csv_path = "data/RDP-Human.csv"
    unprocessed_runs = get_unprocessed_runs(csv_path, processed_runs)

    with open("data/unprocessed_runs.csv", "w") as f:
        f.write("Run,study_accession\n")
        for run in unprocessed_runs:
            f.write(f"{run},{run[:6]}\n")


if __name__ == "__main__":
    main()