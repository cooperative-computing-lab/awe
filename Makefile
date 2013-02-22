
include Makefile.config


EXES = executables

ASSIGN = $(EXES)/assign
CASSIGN = cassign/assign

PYTHON_SRC = awe/*.py trax/*.py executables/*.py


.PHONY: assign
assign : $(ASSIGN)

$(ASSIGN) : $(CASSIGN)
	cp $(CASSIGN) $(ASSIGN)
	chmod +x $(ASSIGN)

$(CASSIGN) : cassign
	make -C $^


.PHONY: build
build : $(PYTHON_SRC)
	python setup.py build


.PHONY: install
install : build
	python setup.py install --prefix $(PREFIX)
