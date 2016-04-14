#ifdef __CINT__
#pragma link off all globals;
#pragma link off all classes;
#pragma link off all functions;

#pragma link C++ class saModuleBase-!;
#pragma link C++ class saSpinAnalyzer-!;
#pragma link C++ class saModLoadSpinInfo-!;
#pragma link C++ class saModLoadSpinInfo::SpinQA-!;

#pragma link C++ class saEventProperty+;
#pragma link C++ class saEventProperty_v1+;

#pragma link C++ class saHist+;
#pragma link C++ class saHistManager-!;

// append user modues here
#pragma link C++ class saModuleSimpleDimuon-!;
#pragma link C++ class saModuleDimuonJpsiHaiwang-!;
#pragma link C++ class saModuleJpsiAN-!;
#pragma link C++ class saModuleDimuonMing-!;
#pragma link C++ class saModuleDimuonDYDarshana-!;
#pragma link C++ class saModuleDimuonDarshana-!;
#pragma link C++ class saModuleSngmuonMing-!;
#pragma link C++ class saModuleSngmuonFeng-!;
#pragma link C++ class saModuleJPsi-!;
#pragma link C++ class saModSimpleHist-!;

#pragma link C++ class saFlag<TArrayI>+;
#pragma link C++ class saFlag<TArrayS>+;
#pragma link C++ class saFlag<TArrayC>+;

#pragma link C++ class saFlagI+;
#pragma link C++ class saFlagS+;
#pragma link C++ class saFlagC+;

#endif /* __CINT__ */
