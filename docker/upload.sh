#!/bin/bash

#Upload docker image onto DNA Nexus
LAB=shapeit5_$(git log -1 --format=%cd --date=short)\_$(git rev-parse --short HEAD)

dx upload $LAB\.tar.gz --path "/docker/"
