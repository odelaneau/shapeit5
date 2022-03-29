projects = phase1 phase2 switch ligate

.PHONY: all $(projects)

all: $(projects)

$(projects):
	$(MAKE) -C $@

clean:
	for dir in $(projects); do \
	$(MAKE) $@ -C $$dir; \
	done
	rm static_bins/*
	rm docker/ressources/*
	rm docker/shapeit5_0.0.1.tar.gz

