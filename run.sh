#!/usr/bin/env bash

if [[ $# -eq 0 ]] ; then
    echo 'Specify directory!'
    exit 0
fi

timestamp() {
   date +%d-%m-%Y_%H-%M-%S
}

destination="${1}"
mkdir $destination

echo "Running..."
echo "Creating dir"
dirname=$(timestamp)

mkdir "${destination}/${dirname}"
mkdir "${destination}/${dirname}/res"

echo "Compiling..."
programname="$(timestamp).out"
g++ -o ${programname} main.cpp

echo "Moving to exec folder..."
cp params.txt "${destination}/${dirname}"

cp p.txt "${destination}/${dirname}"
cp u.txt "${destination}/${dirname}"
cp v.txt "${destination}/${dirname}"

mv ${programname} "${destination}/${dirname}"
cd ${destination}/${dirname} || { echo "No such dir"; exit 1; }

echo "Executing in background..."
./${programname} &

process_id=`/bin/ps -fu $USER| grep ${programname} | grep -v "grep" | awk '{print $2}'`

echo "PID of the process: ${process_id}"
echo "PID of the process: ${process_id}" >> processPID.txt
