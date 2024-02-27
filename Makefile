prep_dirs:
	mkdir -p bin bin/exception bin/utility bin/nonltr
src prep_dirs:
	cd src && ${MAKE}

clean:
	cd src && ${MAKE} clean

all: src

