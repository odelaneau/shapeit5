#!/bin/bash

#Upload docker image onto DNA Nexus
dx rm "/docker/shapeit5_220824_c85eb0c.tar.gz" -f
dx upload shapeit5_220824_c85eb0c.tar.gz --path "/docker/"
