projects = phase_common phase_rare switch ligate

.PHONY: all $(projects) 

all: $(projects)

$(projects):
	+$(MAKE) -C $@ $(MAKECMDGOALS)

clean:
	for dir in $(projects); do \
	$(MAKE) $@ -C $$dir; \
	done
	rm -f static_bins/*
	rm -f docker/resources/*
	rm -f docker/shapeit5*.tar.gz

