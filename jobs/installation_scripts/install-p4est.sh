## MPI and HPC dependent
CONFIGURE_FLAGS='--with-cc=mpicc --with-cxx=mpicxx --with-fc=mpif90'
I_MPI_CC=icx
I_MPI_CXX=icpx
I_MPI_F90=ifx

## Version
P4EST_VERSION="2.8.5"
PETSC_VERSION="3.19.5"

## Install
CURR_DIR=$(pwd)
PACKAGE=p4est
INSTALL_ROOT=$HOME/bin
P4EST_INSTALL=$INSTALL_ROOT/$PACKAGE/$P4EST_VERSION-$MPI_VERSION
TAR_FILE=$PACKAGE-$P4EST_VERSION.tar.gz
URL="https://github.com/p4est/p4est.github.io/raw/master/release"
ROOT_DIR=/tmp
SOURCES_DIR=$ROOT_DIR/$PACKAGE-$P4EST_VERSION
BUILD_DIR=$SOURCES_DIR/build
wget -q $URL/$TAR_FILE -O $ROOT_DIR/$TAR_FILE
mkdir -p $SOURCES_DIR
tar xzf $ROOT_DIR/$TAR_FILE -C $SOURCES_DIR --strip-components=1
cd $SOURCES_DIR
./configure --prefix=$P4EST_INSTALL $CONFIGURE_FLAGS --without-blas --without-lapack --enable-mpi --disable-dependency-tracking
make --quiet
make --quiet install
rm -rf $ROOT_DIR/$TAR_FILE $SOURCES_DIR
cd $CURR_DIR
