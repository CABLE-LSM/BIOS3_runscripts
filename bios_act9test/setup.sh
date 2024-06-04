#!/bin/bash
set -e

mkdir outputs restart input logs

while read path; do
    ln -s $path input
done < inputs.config

while read path; do
    ln -s $path/* input
done < restart.config

while read exe_path; do
    ln -s $exe_path .
done < exe.config
