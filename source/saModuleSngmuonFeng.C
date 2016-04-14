// $Id: saModuleSngmuonFeng.C,v 1.9 2014/05/24 19:54:24 weifeng Exp $  
//
// modified from Jin's saModuleSimpleDimuon  07/12/2013
// for Run12+ single muon analyses: hadron, heavy flavor and W
// 
////////////////////////////////////////////////////////////////////

#include <Fun4AllReturnCodes.h>
#include <PHCompositeNode.h>
#include <PHNodeIterator.h>
#include <getClass.h>
#include <MWGConsts.h>
#include <Tools.h>
#include <Fun4AllServer.h>
#include <Fun4AllHistoManager.h>
#include <TH1F.h>
#include <TrigRunLvl1.h>
#include <TH3.h>
#include <TH2.h>
#include <TF1.h>
#include <SyncObject.h>
#include <RunHeader.h>
#include <TMath.h>
#include <SingleMuonContainer.h>
#include <SingleMuon.h>
#include <RpcSingleMuonContainer.h>
#include <RpcSingleMuon.h>
#include <MutrRefitSingleMuonContainer.h>
#include <MutrRefitSingleMuon.h>

#include <PHGlobal.h>
#include <PHMuoTracksOut.h>
#include <TrigLvl1v1.h>
#include <MWGVertex.h>

#include <limits>
#include <vector>
#include <map>
#include <TGraphErrors.h>
#include <TFile.h>

#include "saHist.h"
#include "saFlag.h"
#include "saModuleSngmuonFeng.h"

using namespace std;

static const int _n_pt_cut_bins = 15;
static const int _n_p_cut_bins = 15;

static const float _pt_cut_bins[_n_pt_cut_bins] = {1.25, 1.5, 1.75, 2.0, 2.25, 2.5, 2.75, 3.0, 3.5, 4.0, 5.0, 6.0, 7.0, 10.0, 15.0};
static const float _p_cut_bins[_n_p_cut_bins] = {2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0, 12.0, 14.0, 16.0, 18.0, 20.0, 40.0};

#ifdef USEPTBINCUT
static const float _ddg0_cuts[2][3][_n_pt_cut_bins] = {
  {
    {14.0, 14.0, 14.0, 14.0, 14.0, 14.0, 12.0, 12.0, 12.0, 12.0, 12.0, 12.0, 12.0, 12.0, 12.0},
    {12.0, 12.0, 10.0, 10.0, 10.0, 10.0, 10.0, 10.0, 10.0, 10.0, 10.0, 10.0, 10.0, 10.0, 10.0},
    {8.0, 6.0, 6.0, 5.0, 5.0, 5.0, 5.0, 5.0, 5.0, 5.0, 5.0, 5.0, 5.0, 5.0, 5.0}
  },
  {
    {14.0, 14.0, 14.0, 14.0, 14.0, 14.0, 14.0, 14.0, 14.0, 14.0, 14.0, 14.0, 14.0, 14.0, 14.0},
    {14.0, 12.0, 10.0, 10.0, 10.0, 10.0, 10.0, 10.0, 10.0, 10.0, 10.0, 10.0, 10.0, 10.0, 10.0},
    {8.0, 6.0, 6.0, 5.0, 5.0, 5.0, 5.0, 5.0, 5.0, 5.0, 5.0, 5.0, 5.0, 5.0, 5.0}
  }
};

static const float _dg0_cuts[2][3][_n_pt_cut_bins] = {
  {
    {32.0, 20.0, 20.0, 20.0, 15.0, 15.0, 15.0, 15.0, 15.0, 15.0, 15.0, 15.0, 15.0, 15.0, 15.0},
    {37.0, 32.0, 20.0, 20.0, 20.0, 20.0, 20.0, 20.0, 20.0, 20.0, 20.0, 20.0, 20.0, 20.0, 20.0},
    {20.0, 15.0, 12.0, 12.0, 12.0, 12.0, 10.0, 10.0, 10.0, 10.0, 10.0, 10.0, 10.0, 10.0, 10.0}
  },
  {
    {16.0, 16.0, 15.0, 15.0, 15.0, 15.0, 15.0, 15.0, 15.0, 15.0, 15.0, 15.0, 15.0, 15.0, 15.0},
    {17.0, 16.0, 16.0, 16.0, 16.0, 15.0, 15.0, 15.0, 15.0, 15.0, 15.0, 15.0, 15.0, 15.0, 15.0},
      {12.0, 12.0, 10.0, 10.0, 10.0, 10.0, 10.0, 10.0, 10.0, 10.0, 10.0, 10.0, 10.0, 10.0, 10.0}
  }
};

static const float _VtxRad_cuts[2][3][_n_pt_cut_bins] = {
  {
    {120.0, 120.0, 100.0, 100.0, 100.0, 100.0, 100.0, 100.0, 100.0, 100.0, 100.0, 100.0,100.0, 100.0, 100.0},
    {120.0, 120.0, 100.0, 100.0, 100.0, 100.0, 100.0, 100.0, 100.0, 100.0, 100.0, 100.0,100.0, 100.0, 100.0},
    {100.0, 100.0, 80.0, 60.0, 60.0, 60.0, 60.0, 60.0, 60.0, 60.0, 60.0, 60.0, 60.0, 60.0, 60.0}
  },
  {
    {180.0, 160.0, 140.0, 140.0, 140.0, 140.0, 140.0, 140.0, 140.0, 140.0, 140.0, 140.0, 140.0, 140.0, 140.0},
    {180.0, 160.0, 140.0, 140.0, 140.0, 140.0, 140.0, 140.0, 140.0, 140.0, 140.0, 140.0, 140.0, 140.0, 140.0},
    {100.0, 100.0, 80.0, 60.0, 60.0, 60.0, 60.0, 60.0, 60.0, 60.0, 60.0, 60.0, 60.0, 60.0, 60.0}
  }
};

static const float _VtxChi2_cuts[2][3][_n_pt_cut_bins] = {
  {
    {0, 15.88, 17.38, 16.12, 14.38, 12.38, 10.12, 8.12, 6.62, 4.88, 3.62, 2.38, 1.62, 1.62, 0.88},
    {0, 15.88, 17.38, 16.12, 14.38, 12.38, 10.12, 8.12, 6.62, 4.88, 3.62, 2.38, 1.62, 1.62, 0.88},
    {0, 15.88, 17.38, 16.12, 14.38, 12.38, 10.12, 8.12, 6.62, 4.88, 3.62, 2.38, 1.62, 1.62, 0.88}
  },
  {
    {0, 14.38, 16.38, 17.62, 15.38, 12.88, 10.38, 8.38, 6.62, 4.88, 3.62, 2.12, 1.62, 1.38, 0.62},
    {0, 14.38, 16.38, 17.62, 15.38, 12.88, 10.38, 8.38, 6.62, 4.88, 3.62, 2.12, 1.62, 1.38, 0.62},
    {0, 14.38, 16.38, 17.62, 15.38, 12.88, 10.38, 8.38, 6.62, 4.88, 3.62, 2.12, 1.62, 1.38, 0.62}
  }
};

#endif

#ifndef USEPTBINCUT
static const float _dg0_cuts[2][3][_n_p_cut_bins] = {
  {
    {0, 35.44, 26.06, 17.27, 13.01, 10.51, 9.21, 8.41, 7.88, 7.41, 6.92, 6.72, 6.49, 6.51, 6.45},
    {0, 35.44, 26.06, 17.27, 13.01, 10.51, 9.21, 8.41, 7.88, 7.41, 6.92, 6.72, 6.49, 6.51, 6.45},
    {0, 35.44, 26.06, 17.27, 13.01, 10.51, 9.21, 8.41, 7.88, 7.41, 6.92, 6.72, 6.49, 6.51, 6.45}
  },
  {
    {0, 11.59, 8.96, 7.45, 6.78, 6.46, 6.33, 6.24, 6.20, 6.15, 6.07, 6.05, 6.03, 6.09, 6.04},
    {0, 11.59, 8.96, 7.45, 6.78, 6.46, 6.33, 6.24, 6.20, 6.15, 6.07, 6.05, 6.03, 6.09, 6.04},
    {0, 11.59, 8.96, 7.45, 6.78, 6.46, 6.33, 6.24, 6.20, 6.15, 6.07, 6.05, 6.03, 6.09, 6.04},
  }
};

static const float _ddg0_cuts[2][3][_n_p_cut_bins] = {
  {
    {0, 13.69, 9.63, 6.51, 5.12, 4.45, 4.08, 3.85, 3.70, 3.56, 3.44, 3.38, 3.34, 3.34, 3.28},
    {0, 13.69, 9.63, 6.51, 5.12, 4.45, 4.08, 3.85, 3.70, 3.56, 3.44, 3.38, 3.34, 3.34, 3.28},
    {0, 13.69, 9.63, 6.51, 5.12, 4.45, 4.08, 3.85, 3.70, 3.56, 3.44, 3.38, 3.34, 3.34, 3.28}
  },
  {
    {0, 11.32, 8.28, 6.08, 4.99, 4.39, 4.05, 3.85, 3.69, 3.55, 3.45, 3.36, 3.32,3.29, 3.23},
    {0, 11.32, 8.28, 6.08, 4.99, 4.39, 4.05, 3.85, 3.69, 3.55, 3.45, 3.36, 3.32,3.29, 3.23},
    {0, 11.32, 8.28, 6.08, 4.99, 4.39, 4.05, 3.85, 3.69, 3.55, 3.45, 3.36, 3.32,3.29, 3.23}
  }
};

static const float _deltaZ_lowcuts[2][3][_n_p_cut_bins] = {
  {
    {0, -1.32, -1.89, -2.22, -2.51, -2.80, -3.04, -2.99, -3.14, -3.15, -3.03, -2.86, -2.68, -2.63, -2.33},
    {0, -1.32, -1.89, -2.22, -2.51, -2.80, -3.04, -2.99, -3.14, -3.15, -3.03, -2.86, -2.68, -2.63, -2.33},
    {0, -1.32, -1.89, -2.22, -2.51, -2.80, -3.04, -2.99, -3.14, -3.15, -3.03, -2.86, -2.68, -2.63, -2.33}
  },
  {
    {0, -1.79, -2.14, -2.77, -3.25, -3.49, -3.75, -3.60, -3.60, -3.45, -3.09, -2.92, -2.72, -2.56, -2.34},
    {0, -1.79, -2.14, -2.77, -3.25, -3.49, -3.75, -3.60, -3.60, -3.45, -3.09, -2.92, -2.72, -2.56, -2.34},
    {0, -1.79, -2.14, -2.77, -3.25, -3.49, -3.75, -3.60, -3.60, -3.45, -3.09, -2.92, -2.72, -2.56, -2.34}
  }
};

static const float _deltaZ_highcuts[2][3][_n_p_cut_bins] = {
  {
    {0, 1.66, 1.91, 2.22, 2.47, 2.72, 3.00, 2.83, 3.05, 3.09, 2.80, 2.76, 2.57, 2.47, 2.27},
    {0, 1.66, 1.91, 2.22, 2.47, 2.72, 3.00, 2.83, 3.05, 3.09, 2.80, 2.76, 2.57, 2.47, 2.27},
    {0, 1.66, 1.91, 2.22, 2.47, 2.72, 3.00, 2.83, 3.05, 3.09, 2.80, 2.76, 2.57, 2.47, 2.27}
  },
  {
    {0, 1.43, 2.11, 2.75, 3.33, 3.59, 3.78, 3.77, 3.68, 3.53, 3.39, 3.06, 2.87, 2.82, 2.44},
    {0, 1.43, 2.11, 2.75, 3.33, 3.59, 3.78, 3.77, 3.68, 3.53, 3.39, 3.06, 2.87, 2.82, 2.44},
    {0, 1.43, 2.11, 2.75, 3.33, 3.59, 3.78, 3.77, 3.68, 3.53, 3.39, 3.06, 2.87, 2.82, 2.44}
  }
};

static const float _VtxChi2_cuts[2][3][_n_p_cut_bins] = {
  {
    {0, 15.88, 17.38, 16.12, 14.38, 12.38, 10.12, 8.12, 6.62, 4.88, 3.62, 2.38, 1.62, 1.62, 0.88},
    {0, 15.88, 17.38, 16.12, 14.38, 12.38, 10.12, 8.12, 6.62, 4.88, 3.62, 2.38, 1.62, 1.62, 0.88},
    {0, 15.88, 17.38, 16.12, 14.38, 12.38, 10.12, 8.12, 6.62, 4.88, 3.62, 2.38, 1.62, 1.62, 0.88}
    //    {0, 9.88, 11.28, 11.93, 11.12, 10.03, 8.0, 7.0, 5.5, 3.5, 2.5, 2.0, 1.5, 0.75, 0.55},
  },
  {
    {0, 14.38, 16.38, 17.62, 15.38, 12.88, 10.38, 8.38, 6.62, 4.88, 3.62, 2.12, 1.62, 1.38, 0.62},
    {0, 14.38, 16.38, 17.62, 15.38, 12.88, 10.38, 8.38, 6.62, 4.88, 3.62, 2.12, 1.62, 1.38, 0.62},
    {0, 14.38, 16.38, 17.62, 15.38, 12.88, 10.38, 8.38, 6.62, 4.88, 3.62, 2.12, 1.62, 1.38, 0.62}
  }
};

static const float _VtxRad_cuts[2][3][_n_pt_cut_bins] = {
  {
    {0, 174.56, 128.68, 86.49, 69.19, 61.12, 56.83, 54.16, 52.44, 50.31, 48.91, 48.56, 47.87, 48.12, 47.75},
    {0, 174.56, 128.68, 86.49, 69.19, 61.12, 56.83, 54.16, 52.44, 50.31, 48.91, 48.56, 47.87, 48.12, 47.75},
    {0, 174.56, 128.68, 86.49, 69.19, 61.12, 56.83, 54.16, 52.44, 50.31, 48.91, 48.56, 47.87, 48.12, 47.75}
  },
  {
    {0, 149.12, 114.84, 84.52, 69.78, 61.61, 56.23, 54.25, 51.30, 50.19, 50.19, 49.11, 48.73, 48.37, 47.89},
    {0, 149.12, 114.84, 84.52, 69.78, 61.61, 56.23, 54.25, 51.30, 50.19, 50.19, 49.11, 48.73, 48.37, 47.89},
    {0, 149.12, 114.84, 84.52, 69.78, 61.61, 56.23, 54.25, 51.30, 50.19, 50.19, 49.11, 48.73, 48.37, 47.89}
  }
};


#endif 



static const char *PartName[] = {"ha","mu"}; 
static const int nParts = 2;

static const char *ChargeName[] = {"p","m"};
static const int nCharges = 2;

static const char *BeamName[] = {"B","Y"};
static const int nBeams = 2;

static const char *FbName[] = {"f","b"}; //forward or backward
static const int nFbs = 2;



//-- Single Muon pT
saModuleSngmuonFeng::saModuleSngmuonFeng(const std::string &name) :
  saModuleBase(name), _nevts(0),
  _nHistAnBins(12000), _HistAnLow(-1.2), _HistAnHigh(1.2),
  _PhiSpinUp(TMath::PiOver2()), _PhiSpinDown(-TMath::PiOver2()),
  _verbosity(0)
{
  static double PtBins[] = {1.0, 1.5, 2.0, 2.5, 3.0, 3.85, 5.0};
  int nPtBins = sizeof(PtBins)/sizeof(double)-1;

  _HistPtBins = PtBins;
  _nHistPtBins = nPtBins;

  static double PzBins[] = {1.0, 5.0, 20.0};
  int nPzBins = sizeof(PzBins)/sizeof(double)-1;

  _HistPzBins = PzBins;
  _nHistPzBins = nPzBins;

}

saModuleSngmuonFeng::~saModuleSngmuonFeng()
{

}

//! global initialization
int
saModuleSngmuonFeng::init(PHCompositeNode *topNode, sa_hist_mangager_ptr hm)
{
  //add a new node to output DST
  // make the flag node which is syncronized with cut on picodst_object
  PHCompositeNode* dstNode = NULL;
  PHNodeIterator nodeItr(topNode);
  dstNode = static_cast<PHCompositeNode*>(nodeItr.findFirst("PHCompositeNode",
      "DST"));
  if (!dstNode)
    {
      dstNode = new PHCompositeNode("DST");
      topNode->addNode(dstNode);
    }

  saFlagC *flags = new saFlagC();

  if (flags)
    {      
      // make a new flag node called DST/SimpleDimuonFlag
      PHIODataNode<PHObject> * node = new PHIODataNode<PHObject>(flags,
			      "SngmuonFengFlag", "PHObject");
      if(!node)
	{
	  cout << "saModuleSngmuonFeng::Init failed to create saFlag Node" << endl;
	  return ABORTRUN;
	}
      
      dstNode->addNode(node);
      cout << "saFlag Node is added with version " << flags->ClassName()
	   << " as " << node->getName() << endl;
    }
  else
    {
      cout << "saModuleSngmuonFeng::Init failed to create saFlag" << endl;
      return ABORTRUN;
    }


  // --- define user hisotgrams ---

  //
  // --- histograms for spin asymmetry study ----
  //

  if(_nHistPtBins==0 || _nHistPzBins==0 || _nHistAnBins==0)
    { 
      cout << "saModuleSngmuonFeng::init - ERROR - no valid PtBins or PzBins or AnBins" << endl;
      return ABORTRUN;
    }
  else
    {
      cout << "saModuleSngmuonFeng::init - nHistPtBins = " << _nHistPtBins << ", nHistPzBins = " << _nHistPzBins << ", nHistAnBins = " << _nHistAnBins << endl;
    }

  TH2D *h_temp;
  saHist *sah_temp;
  char histname[200];

  for(int ipart=0; ipart<2; ipart++)
    for(int icharge=0; icharge<2; icharge++)
      for(int ibeam=0; ibeam<2; ibeam++)
	for(int ifb=0; ifb<2; ifb++)
	  {
	    // For Max Likelihood Method
	    sprintf(histname, "h2_lh_pt_an_%s_%s_%s_%s", PartName[ipart], ChargeName[icharge], BeamName[ibeam], FbName[ifb]);
	    h_temp = new TH2D(histname, histname, _nHistPtBins, _HistPtBins, _nHistAnBins, _HistAnLow, _HistAnHigh);	    
	    sah_temp = new saHist(h_temp, saHist::RUN_PROPERTY | saHist::FILL_PROPERTY, Verbosity());
	    hm.get()->registerHisto(sah_temp);
	    
	    sprintf(histname, "h2_lh_pt_anerr2_%s_%s_%s_%s", PartName[ipart], ChargeName[icharge], BeamName[ibeam], FbName[ifb]);
	    h_temp = new TH2D(histname, histname, _nHistPtBins, _HistPtBins, _nHistAnBins, _HistAnLow, _HistAnHigh);
	    sah_temp = new saHist(h_temp, saHist::RUN_PROPERTY | saHist::FILL_PROPERTY, Verbosity());
	    hm.get()->registerHisto(sah_temp);

	    sprintf(histname, "h2_lh_pz_an_%s_%s_%s_%s", PartName[ipart], ChargeName[icharge], BeamName[ibeam], FbName[ifb]);
	    h_temp = new TH2D(histname, histname, _nHistPzBins, _HistPzBins, _nHistAnBins, _HistAnLow, _HistAnHigh);
	    sah_temp = new saHist(h_temp, saHist::RUN_PROPERTY | saHist::FILL_PROPERTY, Verbosity());
	    hm.get()->registerHisto(sah_temp);
	    
	    sprintf(histname, "h2_lh_pz_anerr2_%s_%s_%s_%s", PartName[ipart], ChargeName[icharge], BeamName[ibeam], FbName[ifb]);
	    h_temp = new TH2D(histname, histname, _nHistPzBins, _HistPzBins, _nHistAnBins, _HistAnLow, _HistAnHigh);
	    sah_temp = new saHist(h_temp, saHist::RUN_PROPERTY | saHist::FILL_PROPERTY, Verbosity());
	    hm.get()->registerHisto(sah_temp);

	    // For SQRT Formula
	    sprintf(histname, "h2_sqrt_pt_count_%s_%s_%s_%s", PartName[ipart], ChargeName[icharge], BeamName[ibeam], FbName[ifb]);
	    h_temp = new TH2D(histname, histname, _nHistPtBins, _HistPtBins, SUM, UPLEFT-0.5, SUM+0.5);
	    sah_temp = new saHist(h_temp, saHist::RUN_PROPERTY | saHist::FILL_PROPERTY, Verbosity());
	    hm.get()->registerHisto(sah_temp);

	    sprintf(histname, "h2_sqrt_pz_count_%s_%s_%s_%s", PartName[ipart], ChargeName[icharge], BeamName[ibeam], FbName[ifb]);
	    h_temp = new TH2D(histname, histname, _nHistPzBins, _HistPzBins, SUM, UPLEFT-0.5, SUM+0.5);
	    sah_temp = new saHist(h_temp, saHist::RUN_PROPERTY | saHist::FILL_PROPERTY, Verbosity());
	    hm.get()->registerHisto(sah_temp);

	    sprintf(histname, "h2_sqrt_pt_cosphi_%s_%s_%s_%s", PartName[ipart], ChargeName[icharge], BeamName[ibeam], FbName[ifb]);
	    h_temp = new TH2D(histname, histname, _nHistPtBins, _HistPtBins, SUM, UPLEFT-0.5, SUM+0.5);
	    sah_temp = new saHist(h_temp, saHist::RUN_PROPERTY | saHist::FILL_PROPERTY, Verbosity());
	    hm.get()->registerHisto(sah_temp);

	    sprintf(histname, "h2_sqrt_pz_cosphi_%s_%s_%s_%s", PartName[ipart], ChargeName[icharge], BeamName[ibeam], FbName[ifb]);
	    h_temp = new TH2D(histname, histname, _nHistPzBins, _HistPzBins, SUM, UPLEFT-0.5, SUM+0.5);
	    sah_temp = new saHist(h_temp, saHist::RUN_PROPERTY | saHist::FILL_PROPERTY, Verbosity());
	    hm.get()->registerHisto(sah_temp);

	    sprintf(histname, "h2_sqrt_pt_varmean_%s_%s_%s_%s", PartName[ipart], ChargeName[icharge], BeamName[ibeam], FbName[ifb]);
	    h_temp = new TH2D(histname, histname, _nHistPtBins, _HistPtBins, SUM, UPLEFT-0.5, SUM+0.5);
	    sah_temp = new saHist(h_temp, saHist::RUN_PROPERTY | saHist::FILL_PROPERTY, Verbosity());
	    hm.get()->registerHisto(sah_temp);

	    sprintf(histname, "h2_sqrt_pz_varmean_%s_%s_%s_%s", PartName[ipart], ChargeName[icharge], BeamName[ibeam], FbName[ifb]);
	    h_temp = new TH2D(histname, histname, _nHistPzBins, _HistPzBins, SUM, UPLEFT-0.5, SUM+0.5);
	    sah_temp = new saHist(h_temp, saHist::RUN_PROPERTY | saHist::FILL_PROPERTY, Verbosity());
	    hm.get()->registerHisto(sah_temp);

	    sprintf(histname, "h2_sqrt_pt_phi_%s_%s_%s_%s", PartName[ipart], ChargeName[icharge], BeamName[ibeam], FbName[ifb]);
	    h_temp = new TH2D(histname, histname, _nHistPtBins, _HistPtBins, 8, -TMath::Pi(), TMath::Pi());
	    sah_temp = new saHist(h_temp, saHist::EVENT_PROPERTY | saHist::SPIN_DEPENDENT, Verbosity());
	    hm.get()->registerHisto(sah_temp);

	    sprintf(histname, "h2_sqrt_pz_phi_%s_%s_%s_%s", PartName[ipart], ChargeName[icharge], BeamName[ibeam], FbName[ifb]);
	    h_temp = new TH2D(histname, histname, _nHistPzBins, _HistPzBins, 8, -TMath::Pi(), TMath::Pi());
	    sah_temp = new saHist(h_temp, saHist::EVENT_PROPERTY | saHist::SPIN_DEPENDENT, Verbosity());
	    hm.get()->registerHisto(sah_temp);
	  }

  return EVENT_OK;
}

//! Run initialization
int
saModuleSngmuonFeng::init_run(PHCompositeNode *topNode,
    sa_hist_mangager_ptr hm)
{
  _nevts = 0;

  return EVENT_OK;
}

//! event method
int
saModuleSngmuonFeng::event(PHCompositeNode *topNode, sa_hist_mangager_ptr hm)
{
  SingleMuonContainer* sngmuons = findNode::getClass<SingleMuonContainer>(topNode,
      "SingleMuonContainer");
  if (!sngmuons)
    {
      cout << "saModuleSngmuonFeng:: SngMuonContainer - ERROR - not in Node Tree"
          << endl;
      return ABORTRUN;
    }

  RpcSingleMuonContainer* rpcsngmuons = findNode::getClass<RpcSingleMuonContainer>(topNode, "RpcSingleMuonContainer");
  if (!rpcsngmuons)
    {
      cout << "saModuleSngmuonFeng:: RpcSingelMuonContainer - ERROR - not in Node Tree"	<< endl;
    }

  MutrRefitSingleMuonContainer* mutrrefitsngmuons = findNode::getClass<MutrRefitSingleMuonContainer>(topNode, "MutrRefitSingleMuonContainer");
  if (!mutrrefitsngmuons)
    {
      cout << "saModuleSngmuonFeng:: MutrRefitSingelMuonContainer - ERROR - not in Node Tree"	<< endl;
    }

  saFlagC *flag = findNode::getClass<saFlagC>(topNode, "SngmuonFengFlag");
  if(!flag)
    {
      cout << "saModuleSngmuonFeng::event - SngmuonFengFlag not in Node Tree" << endl;
      return ABORTRUN;
    }

  flag->flags.Set(sngmuons->get_nSingleMuons());
  flag->flags.Reset(0);

  const saEventProperty * ep = findNode::getClass<saEventProperty>(topNode,
      "saEventProperty");
  if (!ep)
    {
      cout << "saModuleSngmuonFeng::event - ERROR  - Cannot find EventProperty node in Top Node"
          << endl;
      return ABORTRUN;
    }

  if (Verbosity() >= 2)
    {
      cout << "saModuleSngmuonFeng::event - INFO  - ";
      ep->identify();
    }  
  
  int run = ep->get_run_number();

  double Pol[2]; // 0 - blue, 1 - yellow
  Pol[0] = ep->get_polarization_blue();
  Pol[1] = ep->get_polarization_yellow();


  _nevts++;
  if(_nevts%1000==0) cout << "process event : " << _nevts << endl;
  
  unsigned int spin = ep->get_spin_state(); //spin: PP - 0, PM - 1, MP - 2, MM - 3, UNKOWN - 4, SUM - 5
  if(spin>3) return ABORTEVENT;

  double PhiSpin[2] = {0, 0}; // 0 - blue, 1 - yellow

  // spin: PP - 0, PM - 1, MP - 2, MM - 3, UNKOWN - 4, SUM - 5
  if(spin==0 || spin==1) //blue spin up
    {
      PhiSpin[0] = _PhiSpinUp;
    }
  else if(spin==2 || spin==3) // blue spin down
    {
      PhiSpin[0] = _PhiSpinDown;
    }

  if(spin==0 || spin==2) //yellow spin up
    { 
      PhiSpin[1] = _PhiSpinUp;
    }
  else if(spin==1 || spin==3) //yellow spin down
    {
      PhiSpin[1] = _PhiSpinDown;
    }


  TrigLvl1* trig_lvl1 = findNode::getClass<TrigLvl1>(topNode,"TrigLvl1");
  if(!trig_lvl1)
    {
      cout <<  PHWHERE << "saModuleSngmuonFeng::event - ERROR -  TrigLvl1 does not exist"<<endl;
      return ABORTRUN;
    }

  bool realtrigSN_1D = false;
  bool realtrigS_SG3_1DH = false;
  bool realtrigN_SG3_1DH = false;

  if (trig_lvl1->get_lvl1_trigscaled()&0x00020000) realtrigSN_1D=true;
  if (trig_lvl1->get_lvl1_trigscaled()&0x00100000) realtrigS_SG3_1DH=true;
  if (trig_lvl1->get_lvl1_trigscaled()&0x00200000) realtrigN_SG3_1DH=true;



    
  bool fill_ok = false;  // output this event to ooutput pdst file if "true"

  double Evt_bbcZ = sngmuons->get_Evt_bbcZ();
  if(fabs(Evt_bbcZ)>30)  return ABORTEVENT;

  for (unsigned int imuon = 0; imuon < sngmuons->get_nSingleMuons(); imuon++)
    {
      SingleMuon* sngmuon = sngmuons->get_SingleMuon(imuon);

      //      int charge = 2*sngmuon->get_charge() -1; // charge in pdst is  "1 for Q+" and "0 for Q-"
      int arm = (sngmuon->get_pz() > 0) ? 1 : 0;
      
      double px = sngmuon->get_px();
      double py = sngmuon->get_py();
      double pz = sngmuon->get_pz();
      
      double phi = atan2(py, px);
      
      double pT = sqrt(px*px +  py*py); 
      double p = sqrt(px*px +  py*py + pz*pz);

      double eta     = sngmuon->get_rapidity();
      int lastGap = sngmuon->get_lastgap();
      double DG0 = sngmuon->get_DG0();
      double DDG0 = sngmuon->get_DDG0();
 
      double ntrhits = sngmuon->get_ntrhits();
      double nidhits = sngmuon->get_nidhits();
      double trchi2 = sngmuon->get_trchi2();
      //      double idchi2 = sngmuon->get_trchi2();
      
      double st1px = sngmuon->get_st1px();
      double st1py = sngmuon->get_st1py();
      double st1pz = sngmuon->get_st1pz();
      double st1p = sqrt(st1px*st1px+st1py*st1py+st1pz*st1pz);
      double costheta = (px*st1px+py*st1py+pz*st1pz)/(p*st1p);
      double dAngle = 0.0;
      if(fabs(costheta)<1.0) dAngle = acos(costheta);
      dAngle = dAngle*0.5*(p+st1p);

      double z_ref = sngmuon->get_z0();
      double x_ref = sngmuon->get_x0() + (Evt_bbcZ-z_ref)*(px/pz);
      double y_ref = sngmuon->get_y0() + (Evt_bbcZ-z_ref)*(py/pz);
      double ref_vtx_r = sqrt(x_ref*x_ref + y_ref*y_ref + z_ref*z_ref);

      // RpcSingleMuon* rpcsngmuon = rpcsngmuons->get_RpcSingleMuon(imuon);
      // double Rpc3DCA = rpcsngmuon->get_Rpc3DCA();

      MutrRefitSingleMuon* mutrrefitsngmuon = mutrrefitsngmuons->get_MutrRefitSingleMuon(imuon);
      double refit_z = mutrrefitsngmuon->get_RfZVtx();
      double delta_z = refit_z - Evt_bbcZ;
      double vtx_chi2 = mutrrefitsngmuon->get_RfVtxChi2();
      //      double vtx_chi2pdf = mutrrefitsngmuon->get_RfVtxChi2PDF();

      //make cuts
 
      if(pz>0) // north arm
	{
	  if(run==360370) continue;
	}
      else // south arm
	{
	  if((run==360372) || (run==361004) || (run==361760) || (run==362807) || (run==362849) || (run==362951)) continue;
	}
      
      if(lastGap<2) continue;
      if(p<2) continue;
      if(fabs(pz) > 20) continue;
      if(pT < 1 || pT > 5) continue;
      if(fabs(eta)<1.4 || fabs(eta)>1.9) continue;
      //      if(fabs(Evt_bbcZ - refit_z) > 2.0) continue;
      if(fabs(ntrhits) < 13.0) continue; 
      //      if (fabs(nidhits) < 6.0) continue; 
      if(pz<0)
	{ 
	  if(fabs(trchi2) > 25.0) continue; 
	}
      else
	{
	  if(fabs(trchi2) > 35.0) continue; 
	}
      // if (fabs(idchi2) > 5.0) continue; 
      
      //      int pt_bin = findBin(pT, _n_pt_cut_bins, _pt_cut_bins);      
      int p_bin = findBin(p, _n_p_cut_bins, _p_cut_bins);
      
      //    cout << pT << " -> " << pt_bin << ", " << p << " -> " << p_bin << endl; 
	  
      // if(DDG0 > _ddg0_cuts[arm][lastGap-2][pt_bin]) continue;
      // if(DG0 > _dg0_cuts[arm][lastGap-2][pt_bin]) continue;
      // if(ref_vtx_r > _VtxRad_cuts[arm][lastGap-2][pt_bin]) continue;
      // if(vtx_chi2 > _VtxChi2_cuts[arm][lastGap-2][p_bin]) continue;

      if(lastGap==4)
	{
	  if(DG0 > _dg0_cuts[arm][lastGap-2][p_bin]) continue;
	  if(DDG0>_ddg0_cuts[arm][lastGap-2][p_bin]) continue;
	  if(ref_vtx_r > _VtxRad_cuts[arm][lastGap-2][p_bin]) continue;    
	  if(vtx_chi2 > _VtxChi2_cuts[arm][lastGap-2][p_bin]) continue;  
	  if(delta_z< _deltaZ_lowcuts[arm][lastGap-2][p_bin] || delta_z> _deltaZ_highcuts[arm][lastGap-2][p_bin]) continue;
	}
      else
	{
	  if(pz<0)	    
	    {
	      if(DG0>12 || DDG0>9) continue;
	    }
	  else
	    {
	      if(DG0>10 || DDG0>9) continue;
	    }
	  if(pT<1.5)
	    {
	      if(vtx_chi2>15) continue;
	    }
	  else
	    {
	      if(vtx_chi2>10) continue;
	    }
	}

      if(lastGap == 2)
	{
	  if(pz<0)
	    {
	      if(fabs(pz)<2.9) continue;
	      if(dAngle>0.35) continue;
	    }
	  else
	    {
	      if(fabs(pz)<3.2) continue;
	      if(dAngle>0.29) continue;
	    }
	}
      else if(lastGap == 3)
	{
	  if(dAngle>0.35) continue;
	  if(pz<0)
	    {
	      if(fabs(pz)<3.1) continue;
	    }
	  else
	    {
	      if(fabs(pz)<3.4) continue;
	    }
	}  
      else if(lastGap == 4)
	{
	  if(!realtrigSN_1D) continue;
	  if(dAngle>0.25) continue;
	  if(fabs(nidhits) < 8.0) continue; 
	}
      
      fill_ok = true;  
      flag->flags[imuon] = 1;  //good track;
      
      //Fill Histogram

      int ipart = lastGap<4 ? 0 : 1;
      int icharge = (sngmuon->get_charge()+1)%2; // charge in pdst is  "1 for Q+" and "0 for Q-"
      int fb[2]; // 0 - blue, 1 - yellow
      fb[0] = pz>0 ? 0 : 1;
      fb[1] = pz>0 ? 1 : 0;

      double modu[2]; // 0 - blue, 1 - yellow
      modu[0] = sin(PhiSpin[0] - phi);
      modu[1] = sin(PhiSpin[1] - phi);

      int sign[2];
      sign[0] = 1;
      sign[1] = -1;

      int sqrt_bin[2];  // 0 - blue, 1 - yellow      
      if(spin==0 || spin==1) //blue spin up
	{
	  if(px>0) // left
	    {
	      sqrt_bin[0] = UPLEFT;
	    }
	  else //right
	    {
	      sqrt_bin[0] = UPRIGHT;
	    }
	}
      else if(spin==2 || spin==3) // blue spin down
	{
	  if(px>0) // left
	    {
	      sqrt_bin[0] = DOWNLEFT;
	    }
	  else //right
	    {
	      sqrt_bin[0] = DOWNRIGHT;
	    }
	}
      
      if(spin==0 || spin==2) //yellow spin up
	{ 
	  if(px>0) // right
	    {
	      sqrt_bin[1] = UPRIGHT;
	    }
	  else //left
	    {
	      sqrt_bin[1] = UPLEFT;
	    }
	}
      else if(spin==1 || spin==3) //yellow spin down
	{
	  if(px>0) // right
	    {
	      sqrt_bin[1] = DOWNRIGHT;
	    }
	  else //left
	    {
	      sqrt_bin[1] = DOWNLEFT;
	    }
	}

      char histname[200];
      saHist *sah_pt_an, *sah_pt_anerr2, *sah_pz_an, *sah_pz_anerr2;

      for(int ibeam=0; ibeam<2; ibeam++)
	{
	  // fill likelihood histogram
	  sprintf(histname, "h2_lh_pt_an_%s_%s_%s_%s", PartName[ipart], ChargeName[icharge], BeamName[ibeam], FbName[fb[ibeam]]);
	  sah_pt_an = hm.get()->getHisto(histname);

	  sprintf(histname, "h2_lh_pt_anerr2_%s_%s_%s_%s", PartName[ipart], ChargeName[icharge], BeamName[ibeam], FbName[fb[ibeam]]);
	  sah_pt_anerr2 = hm.get()->getHisto(histname);

	  sprintf(histname, "h2_lh_pz_an_%s_%s_%s_%s", PartName[ipart], ChargeName[icharge], BeamName[ibeam], FbName[fb[ibeam]]);
	  sah_pz_an = hm.get()->getHisto(histname);

	  sprintf(histname, "h2_lh_pz_anerr2_%s_%s_%s_%s", PartName[ipart], ChargeName[icharge], BeamName[ibeam], FbName[fb[ibeam]]);
	  sah_pz_anerr2 = hm.get()->getHisto(histname);
	  
	  double AnWidth = (_HistAnHigh - _HistAnLow)/_nHistAnBins;
	  //	  double limit = 1.0/Pol[ibeam];
	  //	  int low_bin = 2+int((-limit-_HistAnLow)*_nHistAnBins/(_HistAnHigh - _HistAnLow));
	  //	  int high_bin = int((limit-_HistAnLow)*_nHistAnBins/(_HistAnHigh - _HistAnLow));
	  //	  for(int i = low_bin; i<high_bin; i++)
	  for(int i=0; i<_nHistAnBins; i++)
	    {
	      double AN = _HistAnLow+AnWidth*(i+0.5); //use bin center val
	      //	      double L = numeric_limits<double>::infinity(); 
	      //	      double temp = 1+Pol[ibeam]*AN*modu[ibeam]*sign[ibeam];
	      //	      if(temp>0) L = -log(temp);
	      double L = -log(1+Pol[ibeam]*AN*modu[ibeam]*sign[ibeam]);
	      double L_err2 = pow(Pol[ibeam]*modu[ibeam]/(1+Pol[ibeam]*AN*modu[ibeam]*sign[ibeam]),2);
	      sah_pt_an->Fill(ep, pT, AN, L);
	      sah_pt_anerr2->Fill(ep, pT, AN, L_err2);
	      sah_pz_an->Fill(ep, fabs(pz), AN, L);
	      sah_pz_anerr2->Fill(ep, fabs(pz), AN, L_err2);
	    }
	  
	  // fill sqrt histogram
	  sprintf(histname, "h2_sqrt_pt_count_%s_%s_%s_%s", PartName[ipart], ChargeName[icharge], BeamName[ibeam], FbName[fb[ibeam]]);
	  hm.get()->getHisto(histname)->Fill(ep, pT, sqrt_bin[ibeam]);
	  hm.get()->getHisto(histname)->Fill(ep, pT, SUM);

	  sprintf(histname, "h2_sqrt_pz_count_%s_%s_%s_%s", PartName[ipart], ChargeName[icharge], BeamName[ibeam], FbName[fb[ibeam]]);
	  hm.get()->getHisto(histname)->Fill(ep, fabs(pz), sqrt_bin[ibeam]);
	  hm.get()->getHisto(histname)->Fill(ep, fabs(pz), SUM);

	  sprintf(histname, "h2_sqrt_pt_cosphi_%s_%s_%s_%s", PartName[ipart], ChargeName[icharge], BeamName[ibeam], FbName[fb[ibeam]]);
	  hm.get()->getHisto(histname)->Fill(ep, pT, sqrt_bin[ibeam], fabs(cos(phi)));
	  hm.get()->getHisto(histname)->Fill(ep, pT, SUM, fabs(cos(phi)));

	  sprintf(histname, "h2_sqrt_pz_cosphi_%s_%s_%s_%s", PartName[ipart], ChargeName[icharge], BeamName[ibeam], FbName[fb[ibeam]]);
	  hm.get()->getHisto(histname)->Fill(ep, fabs(pz), sqrt_bin[ibeam], fabs(cos(phi)));
	  hm.get()->getHisto(histname)->Fill(ep, fabs(pz), SUM, fabs(cos(phi)));

	  sprintf(histname, "h2_sqrt_pt_varmean_%s_%s_%s_%s", PartName[ipart], ChargeName[icharge], BeamName[ibeam], FbName[fb[ibeam]]);
	  hm.get()->getHisto(histname)->Fill(ep, pT, sqrt_bin[ibeam], pT);
	  hm.get()->getHisto(histname)->Fill(ep, pT, SUM, pT);

	  sprintf(histname, "h2_sqrt_pz_varmean_%s_%s_%s_%s", PartName[ipart], ChargeName[icharge], BeamName[ibeam], FbName[fb[ibeam]]);
	  hm.get()->getHisto(histname)->Fill(ep, fabs(pz), sqrt_bin[ibeam], fabs(pz));
	  hm.get()->getHisto(histname)->Fill(ep, fabs(pz), SUM, fabs(pz));

	  sprintf(histname, "h2_sqrt_pt_phi_%s_%s_%s_%s", PartName[ipart], ChargeName[icharge], BeamName[ibeam], FbName[fb[ibeam]]);
	  hm.get()->getHisto(histname)->Fill(ep, pT, phi);

	  sprintf(histname, "h2_sqrt_pz_phi_%s_%s_%s_%s", PartName[ipart], ChargeName[icharge], BeamName[ibeam], FbName[fb[ibeam]]);
	  hm.get()->getHisto(histname)->Fill(ep, fabs(pz), phi);

	}//ibeam
            
      
    }//for imuon
  
  return fill_ok? EVENT_OK : ABORTEVENT;
}

//! global termination
int
saModuleSngmuonFeng::end(PHCompositeNode *topNode, sa_hist_mangager_ptr hm)
{

  return EVENT_OK;
}

void
saModuleSngmuonFeng::setHistPtBins(Int_t nbins, Double_t *bins)
{
  _nHistPtBins = nbins;
  _HistPtBins = bins;
}

void
saModuleSngmuonFeng::setHistPzBins(Int_t nbins, Double_t *bins)
{
  _nHistPzBins = nbins;
  _HistPzBins = bins;
}

void
saModuleSngmuonFeng::setHistAnBins(Int_t nbins, Double_t low, Double_t high)
{
  _nHistAnBins = nbins;
  _HistAnLow = low;
  _HistAnHigh = high;
}

int 
saModuleSngmuonFeng::findBin(float val, const int nbins, const float *cut_bins) const
{
  int start = 0;
  int end = nbins-1;

  if(val < cut_bins[start]) return start;

  while(start<end)
    {
      int mid = (start + end)/2;
      if(val<cut_bins[mid]) end = mid;
      else if(start<mid) start = mid;
      else
	{ 
	  start = end;
	  break;	
	}
    }

  return start;
}


int saModuleSngmuonFeng::get_MaxLikelihood(const TH1D *h1) const
{
  int left = 1;
  int right = h1->GetNbinsX();
  
  int result = find_MinVal_Bin(h1, left, right);

  return result;
}

int saModuleSngmuonFeng::find_MinVal_Bin(const TH1D *h1, int left, int right) const
{
  if(left>right) return find_MinVal_Bin(h1, right, left);
  if(left==right) return left;
  if(left+1==right) return h1->GetBinContent(left) < h1->GetBinContent(right) ? left : right;

  int mid = (left+right)/2;

  double v = h1->GetBinContent(mid);
  double lv = h1->GetBinContent(mid-1);
  double rv = h1->GetBinContent(mid+1);
  
  if(lv < v ) return find_MinVal_Bin(h1, left, mid-1);
  else if(rv < v) return find_MinVal_Bin(h1, right, mid+1);
  
  return mid;
  
}


const TH2D* saModuleSngmuonFeng::calcFillAN_lh(saHist* sah_an, saHist* sah_anerr2, saHist* sah_lumi)
{  
  assert(sah_an);
  assert(sah_anerr2);
  assert(sah_lumi);

  vector<int> v_fill = sah_lumi->get_fill_list();
  
  TH2D *h2_base = (TH2D*) sah_an->getHisto("BASE");
  int nbinxs = h2_base->GetNbinsX();
  double *binx = new double[nbinxs+1];
  for(int i=0; i<=nbinxs; i++) 
    binx[i] = h2_base->GetXaxis()->GetBinUpEdge(i);


  string histname = sah_an->get_name();
  unsigned int pos = histname.find("an");
  if(pos == string::npos)
    {
      cout << "saModuleSngmuonFeng::calcAN_lh_FBF - ERROR: invalid sah_an name" << endl;

      delete[] binx;//otherwise it leads to memory leaks
      return NULL;
    }
  histname.replace(pos, 2, "FillAN");

  pos = histname.find("_saHist");
  if(pos!=string::npos) histname.erase(pos, 7);

  TH2D *h2_AN = new TH2D(histname.c_str(), histname.c_str(), nbinxs, binx, v_fill.back()-v_fill.front()+1, v_fill.front()-0.5, v_fill.back()+0.5);

  bool is_blue = true;
  if(histname.find("_B_")==string::npos) is_blue = false;
 
  for(unsigned int ifill=0; ifill<v_fill.size(); ifill++)
    {
      unsigned int fill_number = static_cast<unsigned int>(v_fill[ifill]);
      
      saHist *sah_an_fill = sah_an->getsaHist_fill(fill_number);
      saHist *sah_anerr2_fill = sah_anerr2->getsaHist_fill(fill_number);     
      saHist *sah_lumi_fill = sah_lumi->getsaHist_fill(fill_number);     

      TH2D* h2 = (TH2D*)sah_an_fill->getHisto("SUM");
      TH2D* h2_anerr2 = (TH2D*)sah_anerr2_fill->getHisto("SUM");
      TH1D* h1_lumi = (TH1D*)sah_lumi_fill->getHisto("SUM");
      TH1D *h1_pol = NULL, *h1_polerr = NULL;
      if(is_blue)
	{
	  h1_pol = (TH1D*)sah_lumi_fill->getHisto("Pol_Blue");
	  h1_polerr = (TH1D*)sah_lumi_fill->getHisto("PolStatErr_Blue");
	}
      else
	{
	  h1_pol = (TH1D*)sah_lumi_fill->getHisto("Pol_Yellow");
	  h1_polerr = (TH1D*)sah_lumi_fill->getHisto("PolStatErr_Yellow");
	}
      double lumi = h1_lumi->GetBinContent(1);
      double pol = h1_pol->GetBinContent(1)/lumi;
      double polerr = h1_polerr->GetBinContent(1)/lumi;

      if(Verbosity()>1)
	{
	  if(is_blue) cout << "Fill " << fill_number << ": Blue Pol = " << pol << " +/- " << polerr << endl;
	  else cout << "Fill " << fill_number << ": Yellow Pol = " << pol << " +/- " << polerr << endl;
	}

      for(int ix=1; ix<=nbinxs; ix++)
	{
	  TH1D *h1 = h2->ProjectionY("lh_an",ix,ix);
	  int ibin = get_MaxLikelihood(h1);    //time complexity O( log(N) )
	  if(Verbosity())
	    {
	      int ibin2 = h1->GetMinimumBin();   // time complexity O(N)
	      cout << "ibin = " << ibin << ", ibin2 = " << ibin2 << endl;			  
	    }
	  double AN = h1->GetBinCenter(ibin);
	  double AN_err2 = 1.0/h2_anerr2->GetBinContent(ix, ibin);
	  delete h1;
	  
	  // ------------apply Statistical Uncertainty Correction for Likelihood AN Histograms----------
	  // ------------delta_AN_true^2 = delta_AN_fill^2 + (Pol*AN*delta_Pol)^2-----------------------            
	  AN_err2 += pow(AN*pol*polerr, 2);	
	  
	  if(Verbosity()>1)
	    {
	      cout << "Likelihood error correction " << AN_err2 << " (original) + " << pow(AN*pol*polerr, 2) << " (correction)" << endl;
	    }
	  int iy = h2_AN->GetYaxis()->FindBin(fill_number);
	  h2_AN->SetBinContent(ix, iy, AN);
	  h2_AN->SetBinError(ix, iy, sqrt(AN_err2));
	  
	}// ix
  
      
    }// ifill	
  delete[] binx; // detelet an array should be delete[]!
  return h2_AN;
}


const TH2D* saModuleSngmuonFeng::calcFillAN_sqrt(saHist* sah_count, saHist* sah_cosphi, saHist* sah_lumi)
{

  assert(sah_count);
  assert(sah_cosphi);
  assert(sah_lumi);

  vector<int> v_fill = sah_lumi->get_fill_list();
  
  TH2D *h2 = (TH2D*) sah_count->getHisto("BASE");
  int nbinxs = h2->GetNbinsX();
  double *binx = new double[nbinxs+1];
  for(int i=0; i<=nbinxs; i++) 
    binx[i] = h2->GetXaxis()->GetBinUpEdge(i);

  string histname = sah_count->get_name();
  unsigned int pos = histname.find("count");
  if(pos == string::npos)
    {
      cout << "saModuleSngmuonFeng::calcAN_sqrt_FBF - ERROR: invalid sah_count name" << endl;

      delete[] binx;//otherwise it leads to memory leaks
      return NULL;
    }
  histname.replace(pos, 5, "FillAN");

  pos = histname.find("_saHist");
  if(pos!=string::npos) histname.erase(pos, 7);

  TH2D *h2_AN = new TH2D(histname.c_str(), histname.c_str(), nbinxs, binx, v_fill.back()-v_fill.front()+1, v_fill.front()-0.5, v_fill.back()+0.5);
  
  bool is_blue = true;
  if(histname.find("_B_")==string::npos) is_blue = false;

  for(unsigned int ifill=0; ifill<v_fill.size(); ifill++)
    {
      unsigned int fill_number = static_cast<unsigned int>(v_fill[ifill]);
      
      saHist *sah_count_fill = sah_count->getsaHist_fill(fill_number);
      saHist *sah_cosphi_fill = sah_cosphi->getsaHist_fill(fill_number);
      saHist *sah_lumi_fill = sah_lumi->getsaHist_fill(fill_number);  
      
      TH2D* h2_count_fill = (TH2D *)sah_count_fill->getHisto("SUM");
      TH2D* h2_cosphi_fill = (TH2D *)sah_cosphi_fill->getHisto("SUM");
      TH1D* h1_lumi = (TH1D*)sah_lumi_fill->getHisto("SUM");
      TH1D *h1_pol = NULL, *h1_polerr = NULL;
      if(is_blue)
	{
	  h1_pol = (TH1D*)sah_lumi_fill->getHisto("Pol_Blue");
	  h1_polerr = (TH1D*)sah_lumi_fill->getHisto("PolStatErr_Blue");
	}
      else
	{
	  h1_pol = (TH1D*)sah_lumi_fill->getHisto("Pol_Yellow");
	  h1_polerr = (TH1D*)sah_lumi_fill->getHisto("PolStatErr_Yellow");
	}
      double lumi = h1_lumi->GetBinContent(1);
      double pol = h1_pol->GetBinContent(1)/lumi;
      double polerr = h1_polerr->GetBinContent(1)/lumi;

      if(Verbosity()>1)
	{
	  if(is_blue) cout << "Fill " << fill_number << ": Blue Pol = " << pol << " +/- " << polerr << endl;
	  else cout << "Fill " << fill_number << ": Yellow Pol = " << pol << " +/- " << polerr << endl;
	}
      for(int ix=1; ix<=nbinxs; ix++)
	{
	  double NUL = h2_count_fill->GetBinContent(ix, UPLEFT);
	  double NUR = h2_count_fill->GetBinContent(ix, UPRIGHT);
	  double NDL = h2_count_fill->GetBinContent(ix, DOWNLEFT);
	  double NDR = h2_count_fill->GetBinContent(ix, DOWNRIGHT);
	  
	  double cosphi_sum = h2_cosphi_fill->GetBinContent(ix, SUM);
	  double count_sum = h2_count_fill->GetBinContent(ix, SUM);
	  double cosphi_mean = cosphi_sum/count_sum;
	  
	  if(NUL*NDR<0.1 || NUR*NDL<0.1 || pol*cosphi_mean<10e-8) continue;
	  double AN = 1.0/(pol*cosphi_mean)*(sqrt(NUL*NDR)-sqrt(NUR*NDL))/(sqrt(NUL*NDR)+sqrt(NUR*NDL));
	  double AN_err = sqrt((NUL*NDL*NDR+NUR*NDL*NDR+NUL*NUR*NDL+NUL*NUR*NDR)
			       /(pow(sqrt(NUL*NDR)+sqrt(NUR*NDL),4)*pow(pol*cosphi_mean, 2))
			       + pow(AN*polerr/pol,2));
	  
	  int iy = h2_AN->GetYaxis()->FindBin(fill_number);
	  h2_AN->SetBinContent(ix, iy, AN);
	  h2_AN->SetBinError(ix, iy, AN_err);
	}//ix
      
    }//ifill
  
  delete[] binx;
  return h2_AN;
}

const TH1D* saModuleSngmuonFeng::calcAN_sqrt(saHist* sah_count, saHist* sah_cosphi, saHist* sah_lumi)
{  
  assert(sah_count);
  assert(sah_cosphi);
  assert(sah_lumi);

  TH2D *h2_base = (TH2D*) sah_count->getHisto("BASE");
  int nbinxs = h2_base->GetNbinsX();
  double *binx = new double[nbinxs+1];
  for(int i=0; i<=nbinxs; i++) 
    binx[i] = h2_base->GetXaxis()->GetBinUpEdge(i);


  string histname = sah_count->get_name();
  histname.replace(0, 2, "h1"); // change "h2" to "h1"

  unsigned int pos = histname.find("count");
  if(pos == string::npos)
    {
      cout << "saModuleSngmuonFeng::calcAN_sqrt - ERROR: invalid sah_count name" << endl;

      delete[] binx;//otherwise it leads to memory leaks
      return NULL;
    }
  histname.replace(pos, 5, "AN");

  pos = histname.find("_saHist");
  if(pos!=string::npos) histname.erase(pos, 7);

  TH1D *h1_AN = new TH1D(histname.c_str(), histname.c_str(), nbinxs, binx);

  bool is_blue = true;
  if(histname.find("_B_")==string::npos) is_blue = false;

  //  TH1D *h1_lumi = (TH1D*)sah_lumi->getHisto("SUM");  
  //  double lumi = h1_lumi->GetBinContent(1);

  
  const TArrayD *arr = getAvgPol(sah_lumi);
  
  double pol = arr->At(0);
  double polerr = arr->At(1);
  if(!is_blue) 
    {
      pol = arr->At(2);
      polerr = arr->At(3);
    }
  delete arr;

  if(Verbosity()>1)
    {
      if(is_blue) cout << "Global Run12: Blue Pol = " << pol << " +/- " << polerr << endl;
      else cout << "Global Run12: Yellow Pol = " << pol << " +/- " << polerr << endl;
    }
  
  TH2D* h2_count = (TH2D *)sah_count->getHisto("SUM");
  TH2D* h2_cosphi = (TH2D *)sah_cosphi->getHisto("SUM");

  for(int ix=1; ix<=nbinxs; ix++)
    {
      double NUL = h2_count->GetBinContent(ix, UPLEFT);
      double NUR = h2_count->GetBinContent(ix, UPRIGHT);
      double NDL = h2_count->GetBinContent(ix, DOWNLEFT);
      double NDR = h2_count->GetBinContent(ix, DOWNRIGHT);
      
      double cosphi_sum = h2_cosphi->GetBinContent(ix, SUM);
      double count_sum = h2_count->GetBinContent(ix, SUM);
      double cosphi_mean = cosphi_sum/count_sum;
      
      if(NUL*NDR<0.1 || NUR*NDL<0.1 || pol*cosphi_mean<10e-8) continue;
      double AN = 1.0/(pol*cosphi_mean)*(sqrt(NUL*NDR)-sqrt(NUR*NDL))/(sqrt(NUL*NDR)+sqrt(NUR*NDL));
      double AN_err = sqrt((NUL*NDL*NDR+NUR*NDL*NDR+NUL*NUR*NDL+NUL*NUR*NDR)
			   /(pow(sqrt(NUL*NDR)+sqrt(NUR*NDL),4)*pow(pol*cosphi_mean, 2))
			   + pow(AN*polerr/pol,2));
      
      h1_AN->SetBinContent(ix, AN);
      h1_AN->SetBinError(ix, AN_err);
    }//ix

  delete[] binx;
  return h1_AN;
  
}

const TH1D* saModuleSngmuonFeng::calcAN_lh(saHist* sah_an, saHist* sah_anerr2, saHist* sah_lumi)
{
  assert(sah_an);
  assert(sah_anerr2);
  assert(sah_lumi);
  
  TH2D *h2_base = (TH2D*) sah_an->getHisto("BASE");
  int nbinxs = h2_base->GetNbinsX();
  double *binx = new double[nbinxs+1];
  for(int i=0; i<=nbinxs; i++) 
    binx[i] = h2_base->GetXaxis()->GetBinUpEdge(i);

  string histname = sah_an->get_name();
  histname.replace(0, 2, "h1"); // change "h2" to "h1"

  unsigned int pos = histname.find("an");
  if(pos == string::npos)
    {
      cout << "saModuleSngmuonFeng::calcAN_lh - ERROR: invalid sah_an name" << endl;

      delete[] binx;//otherwise it leads to memory leaks
      return NULL;
    }
  histname.replace(pos, 2, "AN");

  pos = histname.find("_saHist");
  if(pos!=string::npos) histname.erase(pos, 7);

  TH1D *h1_AN = new TH1D(histname.c_str(), histname.c_str(), nbinxs, binx);

  bool is_blue = true;
  if(histname.find("_B_")==string::npos) is_blue = false;

  const TArrayD *arr = getAvgPol(sah_lumi);
  
  double pol = arr->At(0);
  double polerr = arr->At(1);
  if(!is_blue) 
    {
      pol = arr->At(2);
      polerr = arr->At(3);
    }
  delete arr;
   
  if(Verbosity()>1)
    {
      if(is_blue) cout << "Global Run12: Blue Pol = " << pol << " +/- " << polerr << endl;
      else cout << "Global Run12: Yellow Pol = " << pol << " +/- " << polerr << endl;
    }
  TH2D* h2 = (TH2D*)sah_an->getHisto("SUM");
  TH2D* h2_anerr2 = (TH2D*)sah_anerr2->getHisto("SUM");
  for(int ix=1; ix<=nbinxs; ix++)
    {
      TH1D *h1 = h2->ProjectionY("lh_an",ix,ix);
      int ibin = get_MaxLikelihood(h1);    //time complexity O( log(N) )
      if(Verbosity())
	{
	  int ibin2 = h1->GetMinimumBin();   // time complexity O(N)
	  cout << "ibin = " << ibin << ", ibin2 = " << ibin2 << endl;			  
	}
      double AN = h1->GetBinCenter(ibin);
      double AN_err2 = 1.0/h2_anerr2->GetBinContent(ix, ibin);
      delete h1;
      
      // ------------apply Statistical Uncertainty Correction for Likelihood AN Histograms----------
      // ------------delta_AN_true^2 = delta_AN_fill^2 + (Pol*AN*delta_Pol)^2-----------------------            
      AN_err2 += pow(AN*pol*polerr, 2);	

      if(Verbosity()>1)
	{
	  cout << "Likelihood error correction " << AN_err2 << " (original) + " << pow(AN*pol*polerr, 2) << " (correction)" << endl;
	}      
      h1_AN->SetBinContent(ix, AN);
      h1_AN->SetBinError(ix, sqrt(AN_err2));
      
    }// ix

  delete[] binx;
  return h1_AN;
  
}

const TObjArray* saModuleSngmuonFeng::FitFillAN(TH2D *h2_AN)
{
  assert(h2_AN);

  int nbins = h2_AN->GetNbinsY();
  double low = h2_AN->GetYaxis()->GetBinLowEdge(1);
  double high = h2_AN->GetYaxis()->GetBinUpEdge(nbins);

  TF1 *f1 = new TF1("f1","[0]", low-5, high+5);

  string histname = h2_AN->GetName();
  histname.replace(0,2,"h1");
  int pos = histname.find("FillAN");
  histname.replace(pos, 6, "FitFillAN");
  TH1D *h1_AN = (TH1D*)h2_AN->ProjectionX(histname.c_str(), 0,0);
  histname.replace(pos, 9, "FitFillChi2");
  TH1D *h1_Chi2 = (TH1D*)h2_AN->ProjectionX(histname.c_str(), 0,0);

  nbins = h2_AN->GetNbinsX();

  for(int i=1; i<=nbins; i++)
    {
      TH1D* h1 = (TH1D*)h2_AN->ProjectionY("_py",i,i);
      h1->Fit(f1,"QNR");
      h1_AN->SetBinContent(i, f1->GetParameter(0));
      h1_AN->SetBinError(i, f1->GetParError(0));
      h1_Chi2->SetBinContent(i, f1->GetChisquare()/f1->GetNDF());
      delete h1;
    }

  TObjArray *AN_slices = new TObjArray(2, 0);
  AN_slices->AddAt(h1_AN, 0);
  AN_slices->AddAt(h1_Chi2, 1);
  
  delete f1;

  return AN_slices;  
}

const TH1D* saModuleSngmuonFeng::combineBYAN(const TH1D *h1_blue, const TH1D *h1_yellow)
{
  assert(h1_blue);
  assert(h1_yellow);

  string histname = h1_blue->GetName();
  int pos = histname.find("_B_");
  histname.replace(pos, 3, "_C_");

  TH1D* h1_combined = (TH1D*)h1_blue->Clone();
  h1_combined->Reset();
  h1_combined->SetNameTitle(histname.c_str(), histname.c_str());

  int nbins = h1_blue->GetNbinsX();
  
  for(int i=1; i<=nbins; i++)
    {
      double bc_blue = h1_blue->GetBinContent(i);
      double be_blue = h1_blue->GetBinError(i);
      double bc_yellow = h1_yellow->GetBinContent(i);
      double be_yellow = h1_yellow->GetBinError(i);
      
      double bc_combined = 0.0;	      
      double be_combined = 0.0;
      if(be_blue>10e-20 || be_yellow>10e-20)
	{
	  bc_combined = (bc_blue*be_yellow*be_yellow + bc_yellow*be_blue*be_blue)/(be_blue*be_blue + be_yellow*be_yellow);
	  be_combined = be_blue*be_yellow/sqrt(be_blue*be_blue + be_yellow*be_yellow);
	}
      h1_combined->SetBinContent(i, bc_combined);
      h1_combined->SetBinError(i, be_combined);
    } 

  return h1_combined;
}

const TGraphErrors* saModuleSngmuonFeng::makeTGE_pt(const TH1D *h1_AN, const TArrayD *pvx, const double offset)
{
  assert(h1_AN);
  assert(pvx);

  int n = pvx->GetSize();

  if(h1_AN->GetNbinsX()!=n) 
    {
      cout << "saModuleSngmuonFeng::makeTGE_pt - ERROR : number of points for x and y are not same" << endl;
      return NULL;
    }
  
  double *x = new double[n];
  double *ex = new double[n];
  double *y = new double[n];
  double *ey = new double[n];

  for(int i=0; i<n; i++)
    {
      x[i] = pvx->At(i) + offset;
      ex[i] = 0.0;
      y[i] = h1_AN->GetBinContent(i+1);
      ey[i] = h1_AN->GetBinError(i+1);
    }

  TGraphErrors *tge = new TGraphErrors(n, x, y, ex, ey);

  
  delete[] x;
  delete[] ex;
  delete[] y;
  delete[] ey;

  return tge;

}


const TGraphErrors* saModuleSngmuonFeng::makeTGE_pz(const TH1D *h1_AN_forward, const TH1D *h1_AN_backward, 
					      const TArrayD *pvx, const double offset)
{
  assert(h1_AN_forward);
  assert(h1_AN_backward);
  assert(pvx);

  int n = pvx->GetSize();
  if(h1_AN_forward->GetNbinsX()!=n || h1_AN_backward->GetNbinsX()!=n) 
    {
      cout << "saModuleSngmuonFeng::makeTGE_pz - ERROR : number of points for x and y are not same" << endl;
      return NULL;
    }
  
  double *x = new double[n*2];
  double *ex = new double[n*2];
  double *y = new double[n*2];
  double *ey = new double[n*2];

  for(int i=0; i<n; i++)
    {      
      x[n+i] = pvx->At(i) + offset;
      ex[n+i] = 0.0;
      y[n+i] = h1_AN_forward->GetBinContent(i+1);
      ey[n+i] = h1_AN_forward->GetBinError(i+1);

      x[n-1-i] = (-1.0) * pvx->At(i) + offset;
      ex[n-1-i] = 0.0;
      y[n-1-i] = h1_AN_backward->GetBinContent(i+1);
      ey[n-1-i] = h1_AN_backward->GetBinError(i+1);
    }

  
  TGraphErrors *tge = new TGraphErrors(n*2, x, y, ex, ey);
  
  delete[] x;
  delete[] ex;
  delete[] y;
  delete[] ey;

  return tge;


}


const TArrayD* saModuleSngmuonFeng::getVarMean(TH2D *h2_varmean, TH2D *h2_count)
{
  assert(h2_varmean);
  assert(h2_count);
  
  h2_varmean->Divide(h2_count);
  
  int nbins = h2_varmean->GetNbinsX();  
  TArrayD *pvx = new TArrayD(nbins);
  for(int ix=0; ix<nbins; ix++)
    {
      double x = h2_varmean->GetBinContent(ix+1, SUM);
      (*pvx)[ix] = x;
    }

  return pvx;

}

const TArrayD* saModuleSngmuonFeng::getAvgPol(saHist* sah_lumi)
{
  assert(sah_lumi);
  
  int lumi_id = 1;

  TH1D *h1_lumi = (TH1D*)sah_lumi->getHisto("SUM");
  TH1D *h1_pol_blue = (TH1D*)sah_lumi->getHisto("Pol_Blue");
  TH1D *h1_pol_yellow = (TH1D*)sah_lumi->getHisto("Pol_Yellow");

  double lumi = h1_lumi->GetBinContent(lumi_id);
  double pol_blue = h1_pol_blue->GetBinContent(lumi_id)/lumi;
  double pol_yellow = h1_pol_yellow->GetBinContent(lumi_id)/lumi;

  double polerr_blue = 0.0;
  double polerr_yellow = 0.0;

  vector<int> v_fill = sah_lumi->get_fill_list();
  
  for(unsigned int i=0; i<v_fill.size(); i++)
    {
      unsigned int fill_number = v_fill[i];
      saHist *sah_lumi_fill = sah_lumi->getsaHist_fill(fill_number);

      TH1D *h1_polerr_blue = (TH1D*)sah_lumi_fill->getHisto("PolStatErr_Blue"); 
      TH1D *h1_polerr_yellow = (TH1D*)sah_lumi_fill->getHisto("PolStatErr_Yellow"); 

      polerr_blue += pow(h1_polerr_blue->GetBinContent(lumi_id), 2);
      polerr_yellow += pow(h1_polerr_yellow->GetBinContent(lumi_id), 2);
    }

  polerr_blue = sqrt(polerr_blue)/lumi;
  polerr_yellow = sqrt(polerr_yellow)/lumi;
  

  TArrayD *arr = new TArrayD(4);
  (*arr)[0] = pol_blue;
  (*arr)[1] = polerr_blue;
  (*arr)[2] = pol_yellow;
  (*arr)[3] = polerr_yellow;

  return arr;
 
}





