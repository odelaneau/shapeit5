projects = phase_common phase_rare switch ligate

.PHONY: all $(projects) 

all: $(projects)

$(projects):
	+$(MAKE) -C $@ $(MAKECMDGOALS)

clean:
	for dir in $(projects); do \
	$(MAKE) $@ -C $$dir; \
	done
	rm static_bins/*
	rm docker/resources/*
	rm docker/shapeit5*.tar.gz

