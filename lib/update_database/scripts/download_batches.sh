#!/bin/bash

# Loop through each file matching the pattern
for file in rcsb_pdb_ids_1c5c0f0a76ee120adda649c4ecfc740b_*.txt; do
    # Extract the batch number from the filename
    batch_number=$(echo "$file" | sed -E 's/.*_([0-9]+)-[0-9]+\.txt/\1/')

    # Create a directory for this batch
    mkdir -p "batch$batch_number"
    
    # Move the file to the new directory and rename it
    mv "$file" "batch$batch_number/batch$batch_number.txt"
    
    # Go into the directory
    cd "batch$batch_number"
    
    # Run the download script in the background
    bash ../batch_download.sh -f "batch$batch_number.txt" &
    
    # Return to the original directory
    cd ..
done

echo "All downloads started in the background."

