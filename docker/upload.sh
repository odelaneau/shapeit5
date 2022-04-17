#!/bin/bash

#Upload docker image onto DNA Nexus
dx rm "/docker/shapeit5_0.0.1.tar.gz" -f
dx upload shapeit5_0.0.1.tar.gz --path "/docker/"
