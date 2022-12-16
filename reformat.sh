#!/bin/bash

snakefmt workflow/rules/
snakefmt workflow/Snakefile
snakefmt tests/

black workflow/scripts/ sunbeamlib/ tests/
