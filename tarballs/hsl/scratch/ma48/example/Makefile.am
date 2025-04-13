AM_FCFLAGS = -I $(top_builddir)/src

TESTS_ENVIRONMENT = $(srcdir)/runtest.sh
EXTRA_DIST = runtest.sh
TESTS = \
	ma48ss \
	ma48ds \
	ma48ss1 \
	ma48ds1

check_PROGRAMS = \
	ma48ss \
	ma48ds \
	ma48ss1 \
	ma48ds1

ma48ss_SOURCES = ma48ss.f
ma48ss_LDADD = $(top_builddir)/src/libma48.a $(LAPACK_LIBS) $(BLAS_LIBS) $(METIS_LIBS) $(FCLIBS)
EXTRA_DIST += ma48ss.output ma48ss.data
ma48ds_SOURCES = ma48ds.f
ma48ds_LDADD = $(top_builddir)/src/libma48.a $(LAPACK_LIBS) $(BLAS_LIBS) $(METIS_LIBS) $(FCLIBS)
EXTRA_DIST += ma48ds.output ma48ds.data
ma48ss1_SOURCES = ma48ss1.f
ma48ss1_LDADD = $(top_builddir)/src/libma48.a $(LAPACK_LIBS) $(BLAS_LIBS) $(METIS_LIBS) $(FCLIBS)
EXTRA_DIST += ma48ss1.output ma48ss1.data
ma48ds1_SOURCES = ma48ds1.f
ma48ds1_LDADD = $(top_builddir)/src/libma48.a $(LAPACK_LIBS) $(BLAS_LIBS) $(METIS_LIBS) $(FCLIBS)
EXTRA_DIST += ma48ds1.output ma48ds1.data
