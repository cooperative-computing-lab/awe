
include Makefile.config


EXES = executables

ASSIGN = $(EXES)/assign
CASSIGN = cassign/assign


.PHONY: assign
assign : $(ASSIGN)

$(ASSIGN) : $(CASSIGN)
	cp $(CASSIGN) $(ASSIGN)
	chmod +x $(ASSIGN)

$(CASSIGN) : cassign
	make -C $^


