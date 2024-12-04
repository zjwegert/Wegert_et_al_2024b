## MPI and HPC dependent
CONFIGURE_FLAGS='--with-cc=mpicc --with-cxx=mpicxx --with-fc=mpif90'
I_MPI_CC=icx
I_MPI_CXX=icpx
I_MPI_F90=ifx

## Version
PETSC_VERSION="3.19.5"

## Install
CURR_DIR=$(pwd)
PACKAGE=petsc
INSTALL_ROOT=$HOME/bin
PETSC_INSTALL=$INSTALL_ROOT/$PACKAGE/$PETSC_VERSION
TAR_FILE=$PACKAGE-$PETSC_VERSION.tar.gz
URL="https://ftp.mcs.anl.gov/pub/petsc/release-snapshots/"
ROOT_DIR=/tmp
mkdir $ROOT_DIR
SOURCES_DIR=$ROOT_DIR/$PACKAGE-$PETSC_VERSION
BUILD_DIR=$SOURCES_DIR/build
wget -q $URL/$TAR_FILE -O $ROOT_DIR/$TAR_FILE
mkdir -p $SOURCES_DIR
tar xzf $ROOT_DIR/$TAR_FILE -C $SOURCES_DIR --strip-components=1
cd $SOURCES_DIR

./configure --prefix=$PETSC_INSTALL $CONFIGURE_FLAGS \
    --download-cmake --download-bison --download-mumps \
    --download-scalapack --download-parmetis --download-metis \
    --download-ptscotch --with-debugging --with-x=0 --with-shared=1 \
    --with-mpi=1 --with-64-bit-indices
make
make install
rm -rf $ROOT_DIR/$TAR_FILE $SOURCES_DIR
cd $CURR_DIR
