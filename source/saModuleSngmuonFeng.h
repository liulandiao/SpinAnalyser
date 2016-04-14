// $Id: saModuleSngmuonFeng.h,v 1.4 2014/04/07 18:49:44 weifeng Exp $
// 
//
// for Run12 p+p 200GeV Transverse Single Spin Analyses: heavy flavor, light hadron
//
/////////////////////////////////////////////////////////////////////

#ifndef SAMODULESNGMUONFENG_H_
#define SAMODULESNGMUONFENG_H_

#include <vector>
#include <map>
#include <string>
#include <TFile.h>
#include <TH1.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TArrayD.h>
#include <TGraphErrors.h>
#include <TObjArray.h>
#include "saModuleBase.h"

class saHist;

class saModuleSngmuonFeng : public saModuleBase
{
public:
  saModuleSngmuonFeng(const std::string &name = "SngMuonFeng");
  virtual
  ~saModuleSngmuonFeng();

  int Verbosity() {return _verbosity;}
  void Verbosity(int ver) {_verbosity = ver;}

  void setHistPtBins(int nbins, double *bins);
  void setHistPzBins(int nbins, double *bins);
  void setHistAnBins(int nbins, double low, double high);


  const TH2D* calcFillAN_lh(saHist* sah_an, saHist* sah_anerr2, saHist* sah_lumi);
  const TH2D* calcFillAN_sqrt(saHist* sah_count, saHist* sah_cosphi, saHist* sah_lumi);

  const TH1D* calcAN_lh(saHist* sah_an, saHist* sah_anerr2, saHist* sah_lumi);
  const TH1D* calcAN_sqrt(saHist* sah_count, saHist* sah_cosphi, saHist* sah_lumi);
 
  const TObjArray* FitFillAN(TH2D *h2_AN);
  
  const TH1D* combineBYAN(const TH1D *h1_blue, const TH1D *h1_yellow);

  const TGraphErrors* makeTGE_pt(const TH1D *h1_AN, const TArrayD *pvx, const double offset = 0.0);  
  const TGraphErrors* makeTGE_pz(const TH1D *h1_AN_forward, const TH1D *h1_AN_backward, 
			   const TArrayD *pvx, const double offset = 0.0);

  const TArrayD* getVarMean(TH2D *h2_varmean, TH2D *h2_count);

  const TArrayD* getAvgPol(saHist *sah_lumi);

  enum SQRT_BIN
  {
    UPLEFT = 1,
    UPRIGHT,
    DOWNLEFT,
    DOWNRIGHT,
    SUM
  };

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

 protected:

  int findBin(float val, const int nbins, const float *cut_bins) const;
  int get_MaxLikelihood(const TH1D *h1) const;
  int find_MinVal_Bin(const TH1D *h1, int left, int right) const;
 
  int _nevts;

  int _nHistPtBins;
  double *_HistPtBins;

  int _nHistPzBins; 
  double *_HistPzBins;

  int _nHistAnBins;
  double _HistAnLow, _HistAnHigh;

  double _PhiSpinUp, _PhiSpinDown;
 
  int _verbosity;

};

#endif /* SAMODULESNGMUONFENG_H_ */
