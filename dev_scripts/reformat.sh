#!/bin/bash

cd $SUNBEAM_DIR

snakefmt workflow/rules/
snakefmt workflow/Snakefile

black workflow/scripts/ src/sunbeamlib/ tests/e2e/*.py tests/unit/rules/*.py tests/unit/sunbeamlib/*.py tests/*.py
