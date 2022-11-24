#!/usr/bin/env bash

for f in *.pdf; do
    DIR="$(dirname "${f}")" ; basename="${f%%.*}"
    pdftoppm -f 1 -l 1 -rx 350 -ry 350 "$f" $"$basename" -png -singlefile
    # inkscape "$f" --export-dpi=600 --export-area-drawing --export-filename="$basename.png"
done