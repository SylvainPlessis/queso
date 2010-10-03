## Process this file with automake to produce Makefile.in

# Build in these directories:
# DIST_SUBDIRS = basic stats

AM_CPPFLAGS = -I$(top_srcdir)/src/core/inc \
	     -I$(top_srcdir)/src/misc/inc  \
	     -I$(top_srcdir)/src/basic/inc \
	     -I$(top_srcdir)/src/stats/inc \
	     -I$(top_srcdir)/src/interface \
	       $(BOOST_CPPFLAGS) \
	       $(GSL_CFLAGS) \
	       $(GRVY_CFLAGS) \
	       $(GLPK_CFLAGS) \
	       $(HDF5_CFLAGS)

BUILT_SOURCES = .license.stamp


# Dynamic library build

lib_LTLIBRARIES      = libqueso.la
libqueso_includedir  = $(prefix)/include
libqueso_la_LDFLAGS  = $(all_libraries) -release $(GENERIC_RELEASE)
libqueso_la_LDFLAGS += $(GRVY_LIBS) $(GSL_LIBS) $(EPETRA_LIBS)
libqueso_la_LDFLAGS += $(BOOST_PROGRAM_OPTIONS_LDFLAGS) $(BOOST_PROGRAM_OPTIONS_LIBS)
LIBS                 = $(EPETRA_LIBS)

# Sources from core/src

libqueso_la_SOURCES = \
	$(top_srcdir)/src/core/src/uqDefines.C \
	$(top_srcdir)/src/core/src/uqEnvironment.C \
	$(top_srcdir)/src/core/src/uqEnvironmentOptions.C \
	$(top_srcdir)/src/core/src/uqVector.C \
	$(top_srcdir)/src/core/src/uqMatrix.C \
	$(top_srcdir)/src/core/src/uqMPI.C \
	$(top_srcdir)/src/core/src/uqTrilinosVector.C \
	$(top_srcdir)/src/core/src/uqTrilinosMatrix.C

# Sources from core/src with gsl conditional 

if UQBT_GSL
libqueso_la_SOURCES += \
	$(top_srcdir)/src/core/src/uqGslVector.C \
	$(top_srcdir)/src/core/src/uqGslMatrix.C
endif

# Sources from misc/src

libqueso_la_SOURCES += \
	$(top_srcdir)/src/misc/src/uq1D1DFunction.C \
	$(top_srcdir)/src/misc/src/uq1DQuadrature.C \
	$(top_srcdir)/src/misc/src/uqComplexFft.C \
	$(top_srcdir)/src/misc/src/uqMiscellaneous.C \
	$(top_srcdir)/src/misc/src/uqRealFft.C

# Sources from misc/src withn gsl conditional

if UQBT_GSL
libqueso_la_SOURCES += $(top_srcdir)/src/misc/src/uqGslAsciiTable.C
endif

# Sources from basic/src

libqueso_la_SOURCES += \
	$(top_srcdir)/src/basic/src/uqSequenceStatisticalOptions.C \
	$(top_srcdir)/src/basic/src/uqTrilinosVectorSpace.C

# Sources from basic/src with gsl conditional

if UQBT_GSL
libqueso_la_SOURCES += $(top_srcdir)/src/basic/src/uqGslVectorSpace.C
endif

# Sources from stats/src

libqueso_la_SOURCES += \
	$(top_srcdir)/src/stats/src/uqFiniteDistribution.C \
	$(top_srcdir)/src/stats/src/uqMetropolisHastingsSG.C \
	$(top_srcdir)/src/stats/src/uqMetropolisHastingsSGOptions.C \
	$(top_srcdir)/src/stats/src/uqMLSampling.C \
	$(top_srcdir)/src/stats/src/uqMLSamplingOptions.C \
	$(top_srcdir)/src/stats/src/uqMLSamplingLevelOptions.C \
	$(top_srcdir)/src/stats/src/uqMonteCarloSGOptions.C \
	$(top_srcdir)/src/stats/src/uqStatisticalInverseProblemOptions.C \
	$(top_srcdir)/src/stats/src/uqStatisticalForwardProblemOptions.C

# Sources from interface 

libqueso_la_SOURCES += \
	$(top_srcdir)/src/interface/basic_classes.cpp \
	$(top_srcdir)/src/interface/basic_interface.cpp

# Headers to install from core/inc

libqueso_include_HEADERS = \
	$(top_srcdir)/src/core/inc/uqDefines.h \
	$(top_srcdir)/src/core/inc/uqEnvironment.h \
	$(top_srcdir)/src/core/inc/uqEnvironmentOptions.h \
	$(top_srcdir)/src/core/inc/uqGslMatrix.h \
	$(top_srcdir)/src/core/inc/uqGslVector.h \
	$(top_srcdir)/src/core/inc/uqMatrix.h \
	$(top_srcdir)/src/core/inc/uqMPI.h \
	$(top_srcdir)/src/core/inc/uqTrilinosMatrix.h \
	$(top_srcdir)/src/core/inc/uqTrilinosVector.h \
	$(top_srcdir)/src/core/inc/uqVector.h

# Headers to install from misc/inc

libqueso_include_HEADERS += \
	$(top_srcdir)/src/misc/inc/uq1D1DFunction.h \
	$(top_srcdir)/src/misc/inc/uq1DQuadrature.h \
	$(top_srcdir)/src/misc/inc/uq2dArrayOfStuff.h \
	$(top_srcdir)/src/misc/inc/uqArrayOfOneDGrids.h \
	$(top_srcdir)/src/misc/inc/uqArrayOfOneDTables.h \
	$(top_srcdir)/src/misc/inc/uqAsciiTable.h \
	$(top_srcdir)/src/misc/inc/uqCovCond.h \
	$(top_srcdir)/src/misc/inc/uqFft.h \
	$(top_srcdir)/src/misc/inc/uqMiscellaneous.h \
	$(top_srcdir)/src/misc/inc/uqOneDGrid.h

# Headers to install from basic/inc

libqueso_include_HEADERS += \
	$(top_srcdir)/src/basic/inc/uqArrayOfSequences.h \
	$(top_srcdir)/src/basic/inc/uqInstantiateIntersection.h \
	$(top_srcdir)/src/basic/inc/uqScalarFunction.h \
	$(top_srcdir)/src/basic/inc/uqScalarFunctionSynchronizer.h \
	$(top_srcdir)/src/basic/inc/uqScalarSequence.h \
	$(top_srcdir)/src/basic/inc/uqSequenceOfVectors.h \
	$(top_srcdir)/src/basic/inc/uqSequenceStatisticalOptions.h \
	$(top_srcdir)/src/basic/inc/uqVectorFunction.h \
	$(top_srcdir)/src/basic/inc/uqVectorFunctionSynchronizer.h \
	$(top_srcdir)/src/basic/inc/uqVectorSequence.h \
	$(top_srcdir)/src/basic/inc/uqVectorSet.h \
	$(top_srcdir)/src/basic/inc/uqVectorSpace.h \
	$(top_srcdir)/src/basic/inc/uqVectorSubset.h

# Headers to install from stats/inc 

libqueso_include_HEADERS += \
	$(top_srcdir)/src/stats/inc/uqFiniteDistribution.h \
	$(top_srcdir)/src/stats/inc/uqJointPdf.h \
	$(top_srcdir)/src/stats/inc/uqMarkovChainPositionData.h \
	$(top_srcdir)/src/stats/inc/uqMetropolisHastingsSG1.h \
	$(top_srcdir)/src/stats/inc/uqMetropolisHastingsSG2.h \
	$(top_srcdir)/src/stats/inc/uqMetropolisHastingsSGOptions.h \
	$(top_srcdir)/src/stats/inc/uqMLSampling1.h \
	$(top_srcdir)/src/stats/inc/uqMLSampling2.h \
	$(top_srcdir)/src/stats/inc/uqMLSampling3.h \
	$(top_srcdir)/src/stats/inc/uqMLSamplingOptions.h \
	$(top_srcdir)/src/stats/inc/uqMLSamplingLevelOptions.h \
	$(top_srcdir)/src/stats/inc/uqModelValidation.h \
	$(top_srcdir)/src/stats/inc/uqMonteCarloSG.h \
	$(top_srcdir)/src/stats/inc/uqMonteCarloSGOptions.h \
	$(top_srcdir)/src/stats/inc/uqScalarCdf.h \
	$(top_srcdir)/src/stats/inc/uqStatisticalForwardProblem.h \
	$(top_srcdir)/src/stats/inc/uqStatisticalForwardProblemOptions.h \
	$(top_srcdir)/src/stats/inc/uqStatisticalInverseProblem.h \
	$(top_srcdir)/src/stats/inc/uqStatisticalInverseProblemOptions.h \
	$(top_srcdir)/src/stats/inc/uqTKGroup.h \
	$(top_srcdir)/src/stats/inc/uqValidationCycle.h \
	$(top_srcdir)/src/stats/inc/uqVectorCdf.h \
	$(top_srcdir)/src/stats/inc/uqVectorMdf.h \
	$(top_srcdir)/src/stats/inc/uqVectorRealizer.h \
	$(top_srcdir)/src/stats/inc/uqVectorRV.h

# Headers to install from interface

libqueso_include_HEADERS += \
	$(top_srcdir)/src/interface/queso.h \
	$(top_srcdir)/src/interface/basic_classes.h \
	$(top_srcdir)/src/interface/basic_interface.h \
        $(top_srcdir)/src/interface/basic_int.h


#---------------------------------
# Embedded license header support
#---------------------------------

.license.stamp: $(top_srcdir)/LICENSE
	$(top_srcdir)/src/common/lic_utils/update_license.pl $(top_srcdir)/LICENSE $(cc_sources) $(h_sources)
	echo 'updated source license headers' >$@