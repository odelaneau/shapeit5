#!/bin/bash

#Upload docker image onto DNA Nexus
LAB=shapeit5_$(date +'%d%m%Y')\_$(git rev-parse --short HEAD)
dx upload $LAB\.tar.gz --path "/docker/"
