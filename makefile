projects = phase_common phase_rare switch ligate simulate xcftools

.PHONY: all docker $(projects)

all: $(projects)

$(projects):
	$(MAKE) -C $@

clean:
	for dir in $(projects); do \
	$(MAKE) $@ -C $$dir; \
	done
	rm -f static_bins/*
	rm -f docker/resources/*
	rm -f docker/shapeit5*.tar.gz

static_exe:
	for dir in $(projects); do \
	$(MAKE) $@ -C $$dir; \
	done

# we pin the date here, so container build time is not a moving target
DOCKER_BDATE=$(shell date +%Y%m%d)
DOCKER_TAG ?= $(shell git describe --tags --broken --dirty --all --long | \
			sed "s,heads/,," | sed "s,tags/,," \
		)_$(shell uname -m)_$(shell uname -s | \
			tr '[:upper:]' '[:lower:]' \
		)
DOCKER_PROJECTS ?= $(projects)

# A multi-stage build is a bit tedious here, because we are copying
# from the directory we are building in. We also want to make things
# as reproducible as possible, so we tag and keep the artifacts
docker:
	@echo "Building docker images: $(DOCKER_PROJECTS)"
	@for p in $(DOCKER_PROJECTS); do \
		docker build --progress=plain \
			-f docker/Dockerfile.dev \
			-t $$p:$(DOCKER_TAG) \
			--build-arg BDATE=$(DOCKER_BDATE) \
			--build-arg BVERSION=$(DOCKER_TAG) \
			--build-arg PROJECT=$$p \
		. ; \
	done
