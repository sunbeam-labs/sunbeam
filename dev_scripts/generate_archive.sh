#!/bin/bash
# $1 is the version for this release (e.g. 4.0.0 or 4.1.3-rc.2)

tar -czvf sunbeam$1.tar.gz etc/ extensions/.placeholder sunbeamlib/ tests/ workflow/ environment.yml install.sh Readme.md setup.py MANIFEST.in
