#!/bin/bash 

setup git v2_15_1 

GITHUB_LOCATION=$1
GENIE_VERSION=master
CONFIGURE_INCL=$3
CONFIGURE_G4=$4
GENIE_CONFIG_DIR=$5
CONFIGURE_HEPMC3=$6
RUN_LOCALLY=$7

git clone $GITHUB_LOCATION Generator 
if [ ! -z "$2" ] ; then 
  cd Generator
  GENIE_VERSION=$2
  git checkout "$GENIE_VERSION" 
  cd ..
fi
echo Requested Genie version $GENIE_VERSION

if [ "$CONFIGURE_INCL" = "true" ] ; then
    source /cvmfs/dune.opensciencegrid.org/products/dune/setup_dune.sh
    setup inclxx  v5_2_9_5b -q e20:prof
    setup lhapdf  v6_5_3    -q e20:prof:p3913
    setup log4cpp v1_1_3d   -q e20:prof
elif [ "$CONFIGURE_G4" = "true" ] ; then
    source /cvmfs/dune.opensciencegrid.org/products/dune/setup_dune.sh
    setup root v6_12_06a -q e17:debug                                               
    setup lhapdf v5_9_1k -q e17:debug                                               
    setup log4cpp v1_1_3a -q e17:debug                                              
    setup pdfsets v5_9_1b                                                           
    setup gdb v8_1                                                                  
    setup git v2_15_1                                                                       
    setup cmake v3_14_3
    setup geant4 v4_10_3_p03e -q debug:e17
else
    setup root v6_12_06a -q e17:debug 
    setup lhapdf v5_9_1k -q e17:debug 
    setup log4cpp v1_1_3a -q e17:debug 

fi

if [ ! -z "$6" ] ; then
    echo Requested HepMC3 library 
    setup hepmc3 v3_2_7 -f Linux64bit+3.10-2.17 -q debug:e26:p3915
fi 

export GENIEBASE=$(pwd)
export GENIE=$GENIEBASE/Generator 
export GENIE_BIN=$GENIE/bin
export GENIE_LIB=$GENIE/lib
export PYTHIA6=$PYTHIA_FQ_DIR/lib 
export LHAPDF5_INC=$LHAPDF_INC 
export LHAPDF5_LIB=$LHAPDF_LIB 
export XSECSPLINEDIR=$GENIEBASE/data 
export GENIE_REWEIGHT=$GENIEBASE/Reweight 
export PATH=$GENIE/bin:$GENIE_REWEIGHT/bin:$PATH 
export LD_LIBRARY_PATH=$GENIE_LIB:$GENIE_REWEIGHT/lib:$LD_LIBRARY_PATH 
unset GENIEBASE 

cd $GENIE 

if [ ! -z "$5" ] ; then
  echo Requested Genie config $GENIE_CONFIG_DIR 
  if [ -d "$GENIE_CONFIG_DIR" ]; then 
      GALGCONF=$GENIE_CONFIG_DIR
      export GALGCONF
      echo Genie config directory $GALGCONF
  else
      if [ -d "$GENIE/$GENIE_CONFIG_DIR" ]; then
	  GALGCONF=$GENIE/$GENIE_CONFIG_DIR
	  export GALGCONF
	  echo Genie config directory $GALGCONF
      fi
  fi
else
  unset GALGCONF
fi

if [ "$CONFIGURE_INCL" = "true" ] ; then
    echo Configuring INCL... 
    ./configure \
	--enable-gsl \
	--with-optimiz-level=O2 \
	--with-log4cpp-inc=${LOG4CPP_INC} \
	--with-log4cpp-lib=${LOG4CPP_LIB} \
	--with-libxml2-inc=${LIBXML2_INC} \
	--with-libxml2-lib=${LIBXML2_LIB} \
	--with-lhapdf5-lib=${LHAPDF_LIB} \
	--with-lhapdf5-inc=${LHAPDF_INC} \
	--with-pythia6-lib=${PYTHIA_LIB} \
	--enable-incl \
	--with-incl-inc=${INCLXX_FQ_DIR}/include/inclxx \
	--with-incl-lib=${INCLXX_FQ_DIR}/lib  \
	--with-boost-inc=${BOOST_FQ_DIR}/include \
	--with-boost-lib=${BOOST_FQ_DIR}/lib \
        --enable-lhapdf6 \
        --disable-lhapdf5
elif [ "$CONFIGURE_G4" = "true" ] ; then
    echo Configuring G4...
    ./configure \
	--enable-gsl \
	--with-optimiz-level=O2 \
	--with-log4cpp-inc=${LOG4CPP_INC} \
	--with-log4cpp-lib=${LOG4CPP_LIB} \
	--with-libxml2-inc=${LIBXML2_INC} \
	--with-libxml2-lib=${LIBXML2_LIB} \
	--with-lhapdf5-lib=${LHAPDF_LIB} \
	--with-lhapdf5-inc=${LHAPDF_INC} \
	--with-pythia6-lib=${PYTHIA_LIB} \
	--enable-geant4
elif [ ! -z "$6" ] ; then
    ./configure \
	--enable-gsl \
	--enable-rwght \
	--with-optimiz-level=O2 \
	--with-log4cpp-inc=${LOG4CPP_INC} \
	--with-log4cpp-lib=${LOG4CPP_LIB} \
	--with-libxml2-inc=${LIBXML2_INC} \
	--with-libxml2-lib=${LIBXML2_LIB} \
	--with-lhapdf5-lib=${LHAPDF_LIB} \
	--with-lhapdf5-inc=${LHAPDF_INC} \
	--with-pythia6-lib=${PYTHIA_LIB} \
	--enable-hepmc3
else 
    ./configure \
	--enable-gsl \
	--enable-rwght \
	--with-optimiz-level=O2 \
	--with-log4cpp-inc=${LOG4CPP_INC} \
	--with-log4cpp-lib=${LOG4CPP_LIB} \
	--with-libxml2-inc=${LIBXML2_INC} \
	--with-libxml2-lib=${LIBXML2_LIB} \
	--with-lhapdf5-lib=${LHAPDF_LIB} \
	--with-lhapdf5-inc=${LHAPDF_INC} \
	--with-pythia6-lib=${PYTHIA_LIB}    
fi

make -j4
