#!/bin/bash

cd $SUNBEAM_DIR

snakefmt workflow/rules/
snakefmt workflow/Snakefile

black workflow/scripts/ sunbeamlib/ tests/
