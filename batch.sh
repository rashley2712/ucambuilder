#!/bin/bash
objectdbcreator.py $1 -d2
postprocessor.py $1 -d2 --xyls
imageaverager.py $1 -d2 -n300
create_html.py $1
wcssolver.py $1 -d2
mergecolours.py $1 -d3

