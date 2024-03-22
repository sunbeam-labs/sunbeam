#!/bin/bash

# NOT YET AT THE POINT WHERE IT SHOULD BE RUN AS A SCRIPT
# THIS SHOULD JUST PROVIDE GUIDANCE FOR DOING IT MANUALLY

docker build -f Dockerfile -t sunbeamlabs/sunbeam:4.2.0 .

docker images # Find image ID in output

docker tag image_id sunbeamlabs/sunbeam:4.2.0 # If not already tagged properly

docker run --rm -it sunbeamlabs/sunbeam:4.2.0 /bin/bash # Test the image

docker push sunbeamlabs/sunbeam:4.2.0