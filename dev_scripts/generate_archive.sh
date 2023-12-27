#!/bin/bash
# $1 is the version for this release (e.g. 4.0.0 or 4.1.3-rc.2)

tar -czv --exclude '__pycache__' -f sunbeam$1.tar.gz etc/ extensions/.placeholder src/sunbeamlib/ stable_env/ tests/ workflow/ environment.yml install.sh README.md pyproject.toml MANIFEST.in pytest.ini