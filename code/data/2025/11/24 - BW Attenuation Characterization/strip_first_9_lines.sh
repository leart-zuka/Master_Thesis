#!/bin/bash

chmod -x ./*.csv
chown -R lz ./*.csv 

if [ ! -d "./cleaned/" ]; then
    mkdir ./cleaned/
fi

for file in ./*.csv; do
    if [ -f "$file" ]; then
        tail -n +9 "${file}" > "$(pwd)/cleaned/$(basename "$file")"
    fi
done
