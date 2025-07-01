#!/bin/bash

# Kill all docker containers
docker ps -q | xargs -r docker kill

# Remove all docker containers
docker ps -a -q | xargs -r docker rm
