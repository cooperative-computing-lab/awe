
include Makefile.config


EXES = executables

ASSIGN = $(EXES)/awe-assign
CASSIGN = cassign/assign

PYTHON_SRC = awe/*.py executables/*.py

PY_SUBMODULES = trax.git

PY_INSTALL_ARGS = --prefix $(PREFIX) --install-data $(PREFIX)/share/awe

AWE_DATAFILES = awe-generic-data.tar.bz2 awe-instance-data.tar.bz2


.PHONY: assign
assign : $(ASSIGN)

$(ASSIGN) : $(CASSIGN)
	cp $(CASSIGN) $(ASSIGN)
	chmod +x $(ASSIGN)

$(CASSIGN) : cassign
	make -C $^


$(AWE_DATAFILES) :
	$(eval dname = $(@:.tar.bz2=))
	tar cf $@ $(dname)

.PHONY: datafiles
datafiles : $(AWE_DATAFILES)


.PHONY: build
build : $(PYTHON_SRC)
	python setup.py build


.PHONY: install.submodules
install.submodules : install.submodules.trax

.PHONY: install.submodules.trax
install.submodules.trax : trax.git
	cd $^ && python setup.py install $(PY_INSTALL_ARGS)

.PHONY: install
install : build install.submodules  $(ASSIGN) $(AWE_DATAFILES)
	python setup.py install $(PY_INSTALL_ARGS)




.PHONY: clean.trax
clean.trax : trax.git
	rm -rf $^/build

.PHONY: clean
clean : clean.trax
	make -C cassign clean
	rm -rf *-workers debug resample transactional.*


.PHONY: update.trax
update.trax : trax.git
	git submodule init   $^
	git submodule update $^


.PHONY: update
update : update.trax
	git pull
