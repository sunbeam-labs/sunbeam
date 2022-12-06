#!/bin/bash

snakefmt rules/
snakefmt Snakefile
snakefmt tests/

black scripts/ sunbeamlib/ tests/
