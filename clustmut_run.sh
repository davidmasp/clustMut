#!/bin/sh

mode=$1

Rscript -e "source(file = system.file('exec/clustmut_$mode.R', package = 'clustMut'))" $@

