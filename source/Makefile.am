AUTOMAKE_OPTIONS = foreign

INCLUDES = -I$(includedir) -I$(OFFLINE_MAIN)/include -I$(ROOTSYS)/include

AM_LDFLAGS = -L$(libdir) -L$(OFFLINE_MAIN)/lib -L$(ROOTSYS)/lib

lib_LTLIBRARIES = \
  libspin_analyzer.la

noinst_HEADERS = \
	saModuleJPsi.h \
	saModuleSimpleDimuon.h \
	saModuleDimuonJpsiHaiwang.h \
	saModuleJpsiAN.h \
	saModuleDimuonMing.h \
	saModuleDimuonDarshana.h \
	saModuleDimuonDYDarshana.h \
	saModuleSngmuonMing.h \
	saModuleSngmuonFeng.h \
	saModSimpleHist.h \
	saModLoadSpinInfo.h \
  SpinAnalyzerLinkDef.h

include_HEADERS = \
	saFlag.h \
	saHist.h \
	saHistManager.h \
	saModuleBase.h \
	saEventProperty.h \
	saEventProperty_v1.h \
	saSpinAnalyzer.h 

libspin_analyzer_la_SOURCES = \
	saFlag.C \
	saHist.C \
	saHistManager.C \
	saEventProperty.C \
	saEventProperty_v1.C \
	saModuleBase.C \
	saModLoadSpinInfo.C \
	saSpinAnalyzer.C \
	saModuleSimpleDimuon.C \
	saModuleDimuonJpsiHaiwang.C \
	saModuleJpsiAN.C \
	saModuleDimuonMing.C \
	saModuleDimuonDarshana.C \
	saModuleDimuonDYDarshana.C \
	saModuleSngmuonMing.C \
	saModuleSngmuonFeng.C \
	saModuleJPsi.C \
	saModSimpleHist.C \
  SpinAnalyzerDict.C

libspin_analyzer_la_LIBADD = \
  -lmutoo_core \
  -lmutoo_base \
  -lmutoo_interface \
  -lfvtxoo_interface \
  -lrpcoo_interface \
  -lmuioo_interface \
  -lMWG_interface \
  -lMWG_tools \
  -lSubsysReco \
  -lreactionplane \
  -lMuidTriggerEmulator \
  -lphool \
  -luspin \
  -lpicodst_object 

# executables needed to make sure libraries are properly linked
noinst_PROGRAMS = \
  testexternals

testexternals_SOURCES = testexternals.C
testexternals_LDADD = libspin_analyzer.la

testexternals.C:
	echo "//*** this is a generated file. Do not commit, do not edit" > $@
	echo "int main()" >> $@
	echo "{" >> $@
	echo "  return 0;" >> $@
	echo "}" >> $@

# Rule for generating CINT dictionaries from class headers
SpinAnalyzerDict.C: \
	${include_HEADERS} \
	${noinst_HEADERS} 
	rootcint -f $@ -c $(DEFAULT_INCLUDES) $(INCLUDES) $^

clean-local:
	rm -f *Dict*
