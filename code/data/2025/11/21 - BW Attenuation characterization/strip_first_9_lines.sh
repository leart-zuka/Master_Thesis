#!/bin/bash

if [ ! -d "./cleaned/" ]; then
    mkdir ./cleaned/
fi

for file in ./*.csv; do
    if [ -f "$file" ]; then
        tail -n +9 "${file}" > "$(pwd)/cleaned/$(basename "$file")"
    fi
done
