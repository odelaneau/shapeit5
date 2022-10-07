#!/bin/bash

#Upload docker image onto DNA Nexus
LAB=shapeit5_$(git log -1 --format=%cd --date=short)\_$(git rev-parse --short HEAD)

dx cd /docker/
dx rm $LAB\.tar.gz
dx upload $LAB\.tar.gz --path "/docker/"
