PLAT = linux
CPP = g++-mp-8

PACKAGES_DIR = ${HOME}/Packages/
ARPACKPP_INC = ${PACKAGES_DIR}/arpackpp/include
BOOST_INC = ${PACKAGES_DIR}/boost

LAPACK_LIB = -L/opt/local/lib/lapack/ -llapack -lblas
ARPACK_LIB = -L/opt/local/lib -larpack
# SUPERLU_LIB = -L/opt/local/lib -lsuperlu

CPP_WARNINGS = -fpermissive -Wall -ansi
# CPP_WARNINGS =  -pedantic-errors
CPP_DEBUG = -g
CPP_OPTIM = -O2 -std=c++11

CPP_INC = -I${ARPACKPP_INC} -I${BOOST_INC}
CPP_FLAGS = ${CPP_DEBUG} ${CPP_WARNINGS} ${CPP_OPTIM}
CPP_LIBS = ${LAPACK_LIB} ${ARPACK_LIB}

TARGET = freefermion

all:${TARGET}
$(TARGET):main.cpp
	${CPP} ${CPP_INC} ${CPP_FLAGS} ${CPP_LIBS} $^ -o $@
clean:
	rm -f ${TARGET} *.o log *.dat
	rm -rf *.dSYM
