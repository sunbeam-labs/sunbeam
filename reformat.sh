#!/bin/bash

snakefmt rules/
snakefmt Snakefile

black scripts/ sunbeamlib/ tests/