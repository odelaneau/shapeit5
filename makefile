projects = phase_common phase_rare switch ligate simulate xcftools

.PHONY: all $(projects)

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

rgc:
	for dir in $(projects); do \
	$(MAKE) $@ -C $$dir; \
	done

