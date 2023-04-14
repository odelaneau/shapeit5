#!/bin/bash

for FILE in $(find . -name "*.pdf"); do
	convert -quality 100 -density 300 $FILE $(dirname $FILE)/$(basename $FILE .pdf)\.png
done
