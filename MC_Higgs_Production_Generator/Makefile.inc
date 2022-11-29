# PYTHIA configuration file.
# Generated on 2022年10月19日 星期三 12时52分54秒 CDT with the user supplied options:
# --with-fastjet3
# --with-root
# --with-hepmc3

# Install directory prefixes.
PREFIX_BIN=/Users/wangyijie/pythia8307/bin
PREFIX_INCLUDE=/Users/wangyijie/pythia8307/include
PREFIX_LIB=/Users/wangyijie/pythia8307/lib
PREFIX_SHARE=/Users/wangyijie/pythia8307/share/Pythia8

# Compilation flags (see ./configure --help for further documentation).
CXX=g++
CXX_COMMON=-O2 -std=c++11 -pedantic -W -Wall -Wshadow -fPIC -pthread
CXX_SHARED=-dynamiclib
CXX_SONAME=-Wl,-dylib_install_name,@rpath/
LIB_SUFFIX=.dylib
OBJ_COMMON=

EVTGEN_USE=false
EVTGEN_CONFIG=
EVTGEN_BIN=
EVTGEN_INCLUDE=
EVTGEN_LIB=

FASTJET3_USE=true
FASTJET3_CONFIG=fastjet-config
FASTJET3_BIN=
FASTJET3_INCLUDE=-I/usr/local/include
FASTJET3_LIB=-L/usr/local/lib -Wl,-rpath,/usr/local/lib -lfastjet

HEPMC2_USE=false
HEPMC2_CONFIG=
HEPMC2_BIN=
HEPMC2_INCLUDE=
HEPMC2_LIB=

HEPMC3_USE=true
HEPMC3_CONFIG=HepMC3-config
HEPMC3_BIN=
HEPMC3_INCLUDE=-I/usr/local/Cellar/hepmc3/3.2.5_2/include
HEPMC3_LIB=-L/usr/local/Cellar/hepmc3/3.2.5_2/lib -Wl,-rpath,/usr/local/Cellar/hepmc3/3.2.5_2/lib -lHepMC3

LHAPDF5_USE=false
LHAPDF5_CONFIG=
LHAPDF5_BIN=
LHAPDF5_INCLUDE=
LHAPDF5_LIB=

LHAPDF6_USE=false
LHAPDF6_CONFIG=
LHAPDF6_BIN=
LHAPDF6_INCLUDE=
LHAPDF6_LIB=

POWHEG_USE=false
POWHEG_CONFIG=
POWHEG_BIN=
POWHEG_INCLUDE=
POWHEG_LIB=

RIVET_USE=false
RIVET_CONFIG=
RIVET_BIN=
RIVET_INCLUDE=
RIVET_LIB=

ROOT_USE=true
ROOT_CONFIG=root-config
ROOT_BIN=/usr/local/Cellar/root/6.26.06_1/bin/
ROOT_INCLUDE=-I/usr/local/Cellar/root/6.26.06_1/include/root
ROOT_LIB=-L/usr/local/Cellar/root/6.26.06_1/lib/root -Wl,-rpath,/usr/local/Cellar/root/6.26.06_1/lib/root -lCore

GZIP_USE=false
GZIP_CONFIG=
GZIP_BIN=
GZIP_INCLUDE=
GZIP_LIB=

PYTHON_USE=false
PYTHON_CONFIG=
PYTHON_BIN=
PYTHON_INCLUDE=
PYTHON_LIB=

MG5MES_USE=false
MG5MES_CONFIG=
MG5MES_BIN=
MG5MES_INCLUDE=
MG5MES_LIB=

OPENMP_USE=false
OPENMP_CONFIG=
OPENMP_BIN=
OPENMP_INCLUDE=
OPENMP_LIB=

MPICH_USE=false
MPICH_CONFIG=
MPICH_BIN=
MPICH_INCLUDE=
MPICH_LIB=

HDF5_USE=false
HDF5_CONFIG=
HDF5_BIN=
HDF5_INCLUDE=
HDF5_LIB=

HIGHFIVE_USE=false
HIGHFIVE_CONFIG=
HIGHFIVE_BIN=
HIGHFIVE_INCLUDE=
HIGHFIVE_LIB=
