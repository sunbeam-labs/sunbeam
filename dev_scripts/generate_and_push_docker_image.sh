#!/bin/bash

# NOT YET AT THE POINT WHERE IT SHOULD BE RUN AS A SCRIPT
# THIS SHOULD JUST PROVIDE GUIDANCE FOR DOING IT MANUALLY

sunbeam init --data_fp tests/data/reads/ projects/docker/
sunbeam run --profile projects/docker/ --containerize > Dockerfile

docker build -f Dockerfile -t sunbeam4.2.0 .

docker images # Find image ID in output

docker tag image_id ctbushman/sunbeam:4.2.0

docker push ctbushman/sunbeam:4.2.0