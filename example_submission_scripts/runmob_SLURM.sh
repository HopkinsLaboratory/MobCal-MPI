#!/bin/sh

# Create a list of input files
ls *.mfj > mobfit.in
input="mobfit.in"
suffix=".mfj"
count=0

while IFS= read -r filename; do
  # Extract the base filename without the suffix
  base_filename="${filename%$suffix}"

  # Create a temporary input file
  echo "$base_filename" > temp.in
  tail -n +2 "$filename" >> temp.in

  # Overwrite the original input file
  cp temp.in "$filename"
  rm temp.in

  # Create a mobcal run file
  echo "$base_filename.mfj" > "$base_filename.run"
  echo "$base_filename.mout" >> "$base_filename.run"

  # Create a runfile with .sh extension for job submission
  runfile="$base_filename.sh"
  echo "#!/bin/bash" > "$runfile"
  echo "#SBATCH --account=def-shopkins" >> "$runfile" #change this to your groups account name
  echo "#SBATCH --nodes=1" >> "$runfile"
  echo "#SBATCH --ntasks-per-node=8" >> "$runfile"
  echo "#SBATCH --mem=8000M" >> "$runfile"
  echo "#SBATCH --time=0-03:00" >> "$runfile"
  echo "srun MobCal_MPI_202.exe $base_filename.run" >> "$runfile"

  # Submit the job using sbatch
  sbatch "$runfile"

  # Increment the count
  count=$((count + 1))
done < "$input"
