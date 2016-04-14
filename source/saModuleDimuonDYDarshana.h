#ifndef SAMODULEDIMUONDYDARSHANA_H_
#define SAMODULEDIMUONDYDARSHANA_H_

#include <vector>
#include <map>
#include <string>

#include <SubsysReco.h>
#include <DiMuonContainer.h>
#include <TFile.h>
#include <TTree.h>
#include <TH1.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TGraphErrors.h>
#include <TCanvas.h>
#include <DiMuon.h>
#include <FvtxPrimVtxContainer.h>
#include <SingleMuonContainer.h>
#include <SingleMuon.h>
#include <DSTReader_v1.h>
#include <PHPythiaHeader.h>
#include <PHPythiaContainerV2.h>
#include <PHPythia.h>
#include "saModuleBase.h"


class saModuleDimuonDYDarshana : public saModuleBase{

public:
  saModuleDimuonDYDarshana(const std::string &name = "DimuonDYDarshana");
  virtual
  ~saModuleDimuonDYDarshana();

	bool _use_2_hit_tracklet;
	bool _use_cut_tracklet_chi2;
	float _cut_tracklet_chi2;
	float _cut_tracklet_dcar;

  bool passesCuts() const;
  double _pT_bins[3];
  //double _MS_bins[3];
  int GetNumberofTracklets(DSTReader *fvtx_trk_map, const DiMuon *dimuon);
  enum EVENT_TYPE {OS_SIG_N, OS_SIG_S, SS_SIG_N, SS_SIG_S,BAD_EVENT};//, OS_ALL_, SS_ALL_
  enum ARM_TYPE {NORTH, SOUTH};

#ifndef __CINT__

  //! global initialization
  virtual int
  init(PHCompositeNode *topNode, sa_hist_mangager_ptr hm);

  //! Run initialization
  virtual int
  init_run(PHCompositeNode *topNode, sa_hist_mangager_ptr hm);

  //! event method
  virtual int
  event(PHCompositeNode *topNode, sa_hist_mangager_ptr hm);

  //! global termination
  virtual int
  end(PHCompositeNode *topNode, sa_hist_mangager_ptr hm);

#endif

//saHist * _h_mass;

saHist * _h_pT_os_N; // opposite sign, north
saHist * _h_pT_os_S; // opposite sign, south

saHist * _h_pT_sm_N; // same sign, north
saHist * _h_pT_sm_S; // same sign, south

//saHist * _h_pT_al_os; // opposite sign all, north
//saHist * _h_pT_al_ss; // opposite sign all, south
/*
saHist * _h_MS_os_N; // opposite sign, north
saHist * _h_MS_os_S; // opposite sign, south

saHist * _h_MS_sm_N; // same sign, north
saHist * _h_MS_sm_S; // same sign, south
*/
/*saHist * _h_MS_al_os; // opposite sign all, north
saHist * _h_MS_al_ss; // opposite sign all, south
*/
private:

  DSTReader *fvtx_trk_map;
  DiMuonContainer *dimuons;
  SingleMuonContainer *singlemuoncontainer;
  DiMuon* dimuon;
  int nTracklets;


};

#endif /* saModuleDimuonDYDarshana_H_ */
