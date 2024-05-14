#!/bin/sh
 ls *.mfj > mobfit.in 
 input="mobfit.in"
 suffix=".mfj"
 while IFS= read -r var
 do
  inp=${var%$suffix}
  echo "running $inp.mfj"
  mpirun -np 4 ./MobCal_MPI_203.exe "$inp.mfj" "$inp.mout" &
  sleep 0.5s
 done < "$input"
 rm $input
