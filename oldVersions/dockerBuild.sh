#!/bin/bash
# One line script to just hold the docker build command
# This needs to be run as sudo as far as I know
# Additionally, this has to be run in the same directory as the Dockerfile!
sudo docker build -t rna-seq-docker - < Dockerfile