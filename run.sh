#!/usr/bin/env bash

timestamp() {
   date +%d-%m-%Y_%H-%M-%S
}

echo "Running..."
echo "Creating dir"
dirname=$(timestamp)

mkdir "../txt_results/${dirname}"
mkdir "../txt_results/${dirname}/res"

echo "Compiling..."
programname="$(timestamp).out"
g++ -o ${programname} main.cpp

echo "Moving to exec folder..."
cp params.txt "../txt_results/${dirname}"
mv ${programname} "../txt_results/${dirname}"
cd ../txt_results/${dirname} || { echo "No such dir"; exit 1; }

echo "Executing in background..."
./${programname} &

process_id=`/bin/ps -fu $USER| grep ${programname} | grep -v "grep" | awk '{print $2}'`

echo "PID of the process: ${process_id}"
echo "PID of the process: ${process_id}" >> processPID.txt
