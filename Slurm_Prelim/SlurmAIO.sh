#!/bin/bash

if grep -q "void ${1%.*}(TString infile, TString outfile)" $1; then
	echo "Creating slurm job from input code $1"
else
	echo "Error: Code requires two input arguments, TString infile and TString outfile, set up to specify the input file path and output file name, respectively."
	exit 0
fi

read -p "Enter filepath to data directory/directories: " search_dir

echo "Searching for .hipo files..."
i=0
declare -a ListOfFiles
if [[ "${search_dir: -1}" == "/" ]]; then
  search_dir=${search_dir%/}
fi

if [[ "${search_dir: -2}" == "/*" ]]; then
  for entry in $search_dir; do
    if [[ "${entry: -5}" == ".hipo" ]]; then
      ((i++))
	    ListOfFiles+=($entry)
    fi
    for entry2 in $entry/*; do
      if [[ "${entry2: -5}" == ".hipo" ]]; then
        ((i++))
				ListOfFiles+=($entry2)
      fi
    done
  done
  for entry in ${search_dir::-2}; do
    if [[ "${entry: -5}" == ".hipo" ]]; then
      ((i++))
      ListOfFiles+=($entry)
    fi
  done

else
  for entry in $search_dir/*; do
    if [[ "${entry: -5}" == ".hipo" ]]; then
      ((i++))
	  ListOfFiles+=($entry)
    fi
  done
fi

if [[ i == 0 ]]; then
  echo "Error: no .hipo files found at given filepath"
  exit 0
else
	run='blank'
	i2=0
	declare -a ListOfRuns
	for ((j=0; j<i; j++)); do
	    run2=${ListOfFiles[$j]}
	    run2=${run2%.evio.*}
	    if [[ $run2 != $run ]]; then
	      run=$run2
				((i2++))
	      ListOfRuns+=($run".*")
	    fi
	done
  echo "Found $i .hipo files in $i2 runs"
fi

read -p "Enter number of runs to process: " totalruns
if [[ totalruns == "all" ]]; then
	totalruns=i2
elif [[ totalruns -lt 1 ]] || [[ totalruns -gt i ]]; then
	echo "invalid number"
	exit 0
fi
wait=1

while [[ $wait == 1 ]]; do
	read -p "Enter .root output filename: " OutputFile
  OutputFile=${OutputFile%.}


	if [[ -f $OutputFile.root ]]; then
		read -p "Warning: $OutputFile.root already exists. Replace?: " response
		if [[ ${response,,}=="yes" ]] || [[ ${response,,}=="y" ]] || [[ ${response,,}=="ok" ]]; then
			wait=0
		elif [[ ${response,,}=="exit" ]] || [[ ${response,,}=="quit" ]] || [[ ${response,,}=="stop" ]]; then
			echo "Aborting"
			exit 0
		fi
  else
    wait=0
	fi
done


variable=$(grep -n "void ${1%.*}(TString infile, TString outfile)" $1)
line=${variable%%:*}
totallines=$(cat $1 | wc -l)

for ((j=0; j<totalruns; j++)); do


  filename="${1%.*}"_"$j"."${1##*.}"
	if [[ -f $filename ]]; then
  	rm $filename
	fi
  touch $filename

  more $1 | head -n $((line-1)) >> $filename
  echo 'void '${1%.*}'_'$j'(TString infile, TString outfile){' >> $filename
  more $1 | tail -n $((totallines-line)) >> $filename

  filename2=process_"${1%.*}"_"$j".sh

	if [[ -f $filename2 ]]; then
  	rm $filename2
	fi
  touch $filename2

  echo 'gROOT->ProcessLine(".L '${1%.*}'_'$j'.C");' >> $filename2
  echo 'gROOT->ProcessLine("'${1%.*}'_'$j'(\"'${ListOfRuns[$j]}'\",\"'$OutputFile'_'$j'.root\")");' >> $filename2

  filename3=slurm_"${1%.*}"_"$j".sh

	if [[ -f $filename3 ]]; then
  	rm $filename3
	fi
  touch $filename3

  echo '#!/bin/bash' >> $filename3
  echo '#SBATCH --nodes=1' >> $filename3
  echo '#SBATCH --ntasks=1' >> $filename3
  echo '#SBATCH --mem-per-cpu=3000' >> $filename3
  echo '#SBATCH --job-name='${1%.*}'_'$j'' >> $filename3
  echo '#SBATCH --partition=production' >> $filename3
  echo '#SBATCH --account=clas12' >> $filename3
  echo '/bin/more '$filename2' | clas12root -b' >> $filename3

  echo submitting job $filename3
  sbatch $filename3

  sleep 5

done
