projects = phase_common phase_rare switch ligate simulate xcftools

.PHONY: all $(projects)

all: $(projects)

$(projects):
	$(MAKE) -C $@

mac:
	for dir in $(projects); do \
	$(MAKE) -C $$dir mac_apple_silicon; \
	done

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

