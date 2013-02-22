
include Makefile.config


EXES = executables

ASSIGN = $(EXES)/awe-assign
CASSIGN = cassign/assign

PYTHON_SRC = awe/*.py executables/*.py

PY_SUBMODULES = trax.git


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


.PHONY: install.submodules
install.submodules : install.submodules.trax

.PHONY: install.submodules.trax
install.submodules.trax : trax.git
	cd $^
	python setup.py install --prefix $(PREFIX)

.PHONY: install
install : build install.submodules
	python setup.py install --prefix $(PREFIX)



.PHONY: clean
clean :
	make -C cassign clean
	rm -rf *-workers debug resample transactional.*
