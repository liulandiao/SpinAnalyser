// $Id: saModuleJPsi.C,v 1.4 2014/04/15 05:20:43 jinhuang Exp $                                                                                             

/*!
 * \file saModuleJPsi.C
 * \brief 
 * \author Aaron Key <keyaaron@unm.edu>
 * \version $Revision: 1.4 $
 * \date $Date: 2014/04/15 05:20:43 $
 */

#include <Fun4AllReturnCodes.h>
#include <PHCompositeNode.h>
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
#include <DiMuonContainer.h>
#include <DiMuon.h>
#include <MCDiMuonContainer.h>
#include <MCDiMuon.h>

#include <cmath>
#include <algorithm>
#include <boost/bind.hpp>

#include "saFlag.h"

#include <Fun4AllReturnCodes.h>
#include <PHNodeIterator.h>
#include <getClass.h>

#include "saModuleJPsi.h"

using namespace std;

using std::sqrt;
using boost::bind;


saModuleJPsi::saModuleJPsi(const std::string &name) :
  saModuleBase(name), 
  _h_mass_pT(NULL), _h_pT(NULL), _h_pT_m(NULL), 
  _h_pT_p(NULL), _h_pT_side(NULL), _h_pT_wide(NULL), 
  _h_mass_pT_FVTX(NULL), _h_pT_FVTX(NULL), _h_pT_m_FVTX(NULL), 
  _h_pT_p_FVTX(NULL), _h_pT_side_FVTX(NULL), _h_pT_wide_FVTX(NULL), 
  _h_mass_pT_sel(NULL), _h_pT_sel(NULL), _h_pT_m_sel(NULL), 
  _h_pT_p_sel(NULL), _h_pT_side_sel(NULL), _h_pT_wide_sel(NULL),
  _dimuon(NULL), _ep(NULL), _histos()
{
}

saModuleJPsi::~saModuleJPsi()
{
}

saHist*
saModuleJPsi::make_saHist(sa_hist_mangager_ptr hm, TH1* h)
{
  saHist* newHist = new saHist(
       //
       h, //the template histogram, any TH1 and derivatives is accepted
       saHist::DEFAULT_FLAG | saHist::RUN_PROPERTY | saHist::FILL_PROPERTY
       | saHist::CROSSING_CHECKS, // flags
       Verbosity() // verbosity
       );
  hm.get()->registerHisto(newHist);
  _histos.push_back(newHist);

  return newHist;
}

//! global initialization
int
saModuleJPsi::init(PHCompositeNode *topNode, sa_hist_mangager_ptr hm)
{
  const int NMASSBINS = 60;
  const int NPTBINS = 3;
  float mass_bins[NMASSBINS+10];//other wise its boundary overflows in the next for loop
  for (int i=0; i < NMASSBINS+1 ; i++)
    mass_bins[i] = 1.0 + i*(10.0-1.0)/NMASSBINS;

  float pT_bins [] = { 0.0, 2.0, 4.0, 10.0 };

  // make a simple histogram
  TH2F * h_mass_pT = new TH2F("InvMass_pT", 
                              "Invariant Mass;Invariant Mass (GeV);"
                              "Transverse Momentum (GeV)", 
                              NMASSBINS,mass_bins,NPTBINS,pT_bins);
  _h_mass_pT = make_saHist(hm, h_mass_pT); 
  
  TH1F * h_pT = new TH1F("pT", 
                         "Transverse Momentum (+-);"
                         "Transverse Momentum (GeV)",
                         NPTBINS,pT_bins);
  _h_pT = make_saHist(hm, h_pT);

  TH1F * h_pT_like = new TH1F("pT_like", 
                              "Transverse Momentum (+-);"
                              "Transverse Momentum (GeV)",
                              NPTBINS,pT_bins);
  _h_pT_like = make_saHist(hm, h_pT_like);

  TH1F * h_pT_wide = new TH1F("pT_wide", 
                              "Wide Transverse Momentum (+-);"
                              "Transverse Momentum (GeV)",
                              NPTBINS, pT_bins);
  _h_pT_wide = make_saHist(hm, h_pT_wide);

  TH1F * h_pT_m = new TH1F("pT_m", 
                           "Transverse Momentum (--);"
                           "Transverse Momentum (GeV)",
                           NPTBINS, pT_bins);
  _h_pT_m = make_saHist(hm, h_pT_m);

  TH1F * h_pT_p = new TH1F("pT_p",
                           "Transverse Momentum (++);"
                           "Transverse Momentum (GeV)",
                           NPTBINS, pT_bins);
  _h_pT_p = make_saHist(hm, h_pT_p);

  TH1F * h_pT_side = new TH1F("pT_side", 
                              "Sideband Transverse Momentum (+-);"
                              "Transverse Momentum (GeV)",
                              NPTBINS, pT_bins);
  _h_pT_side = make_saHist(hm, h_pT_side);

  // Histograms for events with FVTX information
  TH2F * h_mass_pT_FVTX = new TH2F("InvMass_pT_FVTX", 
                                   "Invariant Mass;Invariant Mass (GeV);"
                                   "Transverse Momentum (GeV)", 
                                   60,1,10,20,0,10);
  _h_mass_pT_FVTX = make_saHist(hm, h_mass_pT_FVTX);

  TH1F * h_pT_FVTX = new TH1F("pT_FVTX",
                              "Transverse Momentum (+-);"
                              "Transverse Momentum (GeV)",
                              NPTBINS, pT_bins);
  _h_pT_FVTX = make_saHist(hm, h_pT_FVTX);

  TH1F * h_pT_like_FVTX = new TH1F("pT_like_FVTX",
                                   "Transverse Momentum (+-);"
                                   "Transverse Momentum (GeV)",
                                   NPTBINS, pT_bins);
  _h_pT_like_FVTX = make_saHist(hm, h_pT_like_FVTX);

  TH1F * h_pT_wide_FVTX = new TH1F("pT_wide_FVTX", 
                                   "Wide Transverse Momentum (+-);"
                                   "Transverse Momentum (GeV)",
                                   NPTBINS, pT_bins);
  _h_pT_wide_FVTX = make_saHist(hm, h_pT_wide_FVTX);

  TH1F * h_pT_m_FVTX = new TH1F("pT_m_FVTX", 
                                "Transverse Momentum (--);"
                                "Transverse Momentum (GeV)",
                                NPTBINS, pT_bins);
  _h_pT_m_FVTX = make_saHist(hm, h_pT_m_FVTX);

  TH1F * h_pT_p_FVTX = new TH1F("pT_p_FVTX", 
                                "Transverse Momentum (++);"
                                "Transverse Momentum (GeV)",
                                NPTBINS, pT_bins);
  _h_pT_p_FVTX = make_saHist(hm, h_pT_p_FVTX);

  TH1F * h_pT_side_FVTX = new TH1F("pT_side_FVTX", 
                                   "Sideband Transverse Momentum (+-);"
                                   "Transverse Momentum (GeV)",
                                   NPTBINS, pT_bins);
  _h_pT_side_FVTX = make_saHist(hm, h_pT_side_FVTX);

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

      // make a new flag node called DST/JPsiFlag
      PHIODataNode<PHObject> * node = new PHIODataNode<PHObject>(flags,
          "JPsiFlag", "PHObject");

      if (!node)
        {

          cout
              << "saModuleJPsi::Init failed to create saEventProperty Node"
              << endl;
          return ABORTRUN;

        }
      else
        {

          dstNode->addNode(node);
          cout << "saFlag Node is added with version " << flags->ClassName()
              << " as " << node->getName() << endl;
        }
    }
  else
    {
      cout << "saModuleJPsi::Init failed to create saEventProperty"
          << endl;
      return ABORTRUN;
    }

  return EVENT_OK;
}

//! Run initialization
int
saModuleJPsi::init_run(PHCompositeNode *topNode,
    sa_hist_mangager_ptr hm)
{

  return EVENT_OK;
}

//! event method
int
saModuleJPsi::event(PHCompositeNode *topNode, sa_hist_mangager_ptr hm)
{

  DiMuonContainer* dimuons = findNode::getClass<DiMuonContainer>(topNode,
      "DiMuonContainer");
  if (!dimuons)
    {
      cout
          << "saModuleJPsi:: DiMuonContainer - ERROR - not in Node Tree"
          << endl;
      return ABORTRUN;
    }

  _ep = findNode::getClass<saEventProperty>(topNode,
      "saEventProperty");
  if (!_ep)
    {
      cout
          << "saModuleJPsi::event - ERROR  - Cannot find EventProperty node in Top Node"
          << endl;
      return ABORTRUN;
    }

  if (Verbosity() >= 2)
    {
      cout << "saModuleJPsi::event - INFO  - ";
      _ep->identify();
    }

  saFlagC *flag = findNode::getClass<saFlagC>(topNode, "JPsiFlag");
  if (!flag)
    {
      cout << "saModuleJPsi::event - JPsiFlag not in Node Tree"
          << endl;
      return ABORTRUN;
    }

  // syncronizely set the cut flag with dimuon container
  flag->flags.Set(dimuons->get_nDiMuons());
  flag->flags.Reset(0);

  // save the event for DST files?
  bool save_event = false;

  // determine this event's primary VTX
      
  // double Evt_fvtxZ = dimuons->get_Evt_fvtxZ();
  // double Evt_fvtxZ2 = dimuons->get_Evt_fvtxZ2();
  double Evt_bbcZ = dimuons->get_Evt_bbcZ();

  if ( fabs(Evt_bbcZ) > 100 ) return DISCARDEVENT;
    
  for (unsigned int imuon = 0; imuon < dimuons->get_nDiMuons(); imuon++)
    {
      _dimuon = dimuons->get_DiMuon(imuon);

      if ( !passesCuts() ) continue;
      
      // Check to see if FVTX information is good
      bool has_fvtx = passesFVTXCuts();
      
      selectOutputHistos(has_fvtx);
      fillHistos();

      save_event = true;
      flag->flags[imuon] = 1;
    }

  return save_event ? EVENT_OK : DISCARDEVENT;
}

//! global termination
int
saModuleJPsi::end(PHCompositeNode *topNode, sa_hist_mangager_ptr hm)
{
  std::for_each(_histos.begin(),_histos.end(),
                bind(&saHist::CalcAsymmetry,_1,
                     hm.get()->getHisto_DefaultLumi(),
                     saHist::BbcNoCutPHENIX));

  std::for_each(_histos.begin(),_histos.end(),
                bind(&saHist::CalcAsymmetry_Chi2Fit,_1,false));
  
  return EVENT_OK;
}

bool
saModuleJPsi::passesCuts() const
{
  // Cuts
  if ( !_dimuon->get_same_event() ) 
    return false;

  int arm0 = (_dimuon->get_Tr0_pz() > 0) ? 1 : 0;
  int arm1 = (_dimuon->get_Tr1_pz() > 0) ? 1 : 0;

  if ( arm0 != arm1 ) 
    return false;
      
  if ( _dimuon->get_Evt_vtxchi2() > 5 ) 
    return false;

  if ( _dimuon->get_Evt_vtxoor() > 1 ) 
    return false;

  if ( sqrt(std::pow(_dimuon->get_X0(),2)
            + std::pow(_dimuon->get_Y0(),2)) > 2 )
    return false;

  if ( _dimuon->get_Tr0_DG0() > 15 
       || _dimuon->get_Tr1_DG0() > 15 )
    return false;
  if ( _dimuon->get_Tr0_DDG0() > 10 
       || _dimuon->get_Tr1_DDG0() > 10 )
    return false;

  if ( _dimuon->get_Tr0_ntrhits() < 10 
       || _dimuon->get_Tr1_ntrhits() < 10 )
    return false;
  if ( _dimuon->get_Tr0_nidhits() < 6 
       || _dimuon->get_Tr1_nidhits() < 6 )
    return false;
      
  if ( _dimuon->get_Tr0_lastgap() < 3 
       || _dimuon->get_Tr1_lastgap() < 3 ) 
    return false;

  if ( _dimuon->get_Tr0_dca_r() > 5 
       || _dimuon->get_Tr1_dca_r() > 5 )
    return false;
     
  double pT = sqrt(std::pow(_dimuon->get_px(),2) 
                   + std::pow(_dimuon->get_py(),2));
     
  if ( pT > 10 ) 
    return false;
  if ( _dimuon->get_pz() > 100 ) 
    return false;

  return true;
}

bool
saModuleJPsi::passesFVTXCuts() const
{
  if ( _dimuon->get_mass_fvtxmutr() == 0 ) 
    return false;

  if ( _dimuon->get_dca_r() > 0.2 )
    return false;

  if ( _dimuon->get_Tr0_chi2_fvtxmutr() > 10 
       || _dimuon->get_Tr1_chi2_fvtxmutr() > 10 )
    return false;

  if ( _dimuon->get_Tr0_dr_fvtx() > 2 
       || _dimuon->get_Tr1_dr_fvtx() > 2 )
    return false;

  if ( fabs(_dimuon->get_Tr0_dtheta_fvtx()) > 0.5 
       || fabs(_dimuon->get_Tr1_dtheta_fvtx()) > 0.5 )
    return false;

  if ( fabs(_dimuon->get_Tr0_dphi_fvtx()) > 0.5 
       || fabs(_dimuon->get_Tr1_dphi_fvtx()) > 0.5 )
    return false;

  return true;
}

void
saModuleJPsi::fillHistos()
{  
  double pT = sqrt(std::pow(_dimuon->get_px(),2) 
                   + std::pow(_dimuon->get_py(),2));

  if (dimuonType() == "+-")
    {
      _h_mass_pT_sel->Fill(_ep, _dimuon->get_mass(), pT);
      if ( _dimuon->get_mass() > 2.6 && _dimuon->get_mass() < 3.6 )
        _h_pT_wide_sel->Fill(_ep, pT);
      if ( _dimuon->get_mass() > 2.8 && _dimuon->get_mass() < 3.44 )
        _h_pT_sel->Fill(_ep, pT);
      else if ( _dimuon->get_mass() > 2.0 && _dimuon->get_mass() < 2.5 )
        _h_pT_side_sel->Fill(_ep, pT);
    }
  else if (dimuonType() == "++")
    {
      if ( _dimuon->get_mass() > 2.6 && _dimuon->get_mass() < 3.6 )
        _h_pT_wide_sel->Fill(_ep, pT, -1);
      if ( _dimuon->get_mass() > 2.8 && _dimuon->get_mass() < 3.44 )
        {
          _h_pT_p_sel->Fill(_ep, pT);              
          _h_pT_sel->Fill(_ep, pT, -1);
          _h_pT_like_sel->Fill(_ep, pT);
        }
      else if ( _dimuon->get_mass() > 2.0 && _dimuon->get_mass() < 2.5 )
        _h_pT_side_sel->Fill(_ep, pT, -1);      
    }
  else if (dimuonType() == "--")
    {
      if ( _dimuon->get_mass() > 2.6 && _dimuon->get_mass() < 3.6 )
        _h_pT_wide_sel->Fill(_ep, pT, -1);
      if ( _dimuon->get_mass() > 2.8 && _dimuon->get_mass() < 3.44 )
        {
          _h_pT_m_sel->Fill(_ep, pT);              
          _h_pT_sel->Fill(_ep, pT, -1);
          _h_pT_like_sel->Fill(_ep, pT);
        }
      else if ( _dimuon->get_mass() > 2.0 && _dimuon->get_mass() < 2.5 )
        _h_pT_side_sel->Fill(_ep, pT, -1);      
    } 
}

void 
saModuleJPsi::selectOutputHistos(bool has_fvtx)
{
  if ( has_fvtx )
    {
      _h_mass_pT_sel = _h_mass_pT_FVTX;
      _h_pT_sel = _h_pT_FVTX;
      _h_pT_like_sel = _h_pT_like_FVTX;
      _h_pT_m_sel = _h_pT_m_FVTX;
      _h_pT_p_sel = _h_pT_p_FVTX;
      _h_pT_side_sel = _h_pT_side_FVTX;
      _h_pT_wide_sel = _h_pT_wide_FVTX;      
    }
  else
    {
      _h_mass_pT_sel = _h_mass_pT;
      _h_pT_sel = _h_pT;
      _h_pT_like_sel = _h_pT_like;
      _h_pT_m_sel = _h_pT_m;;
      _h_pT_p_sel = _h_pT_p;
      _h_pT_side_sel = _h_pT_side;
      _h_pT_wide_sel = _h_pT_wide;      
    }
}

std::string
saModuleJPsi::dimuonType() const
{
  std::string type;
  if (_dimuon->get_charge() == 0 && _dimuon->get_same_event())
    type = "+-";
  if (_dimuon->get_charge() < 0 && _dimuon->get_same_event())
    type = "--";
  if (_dimuon->get_charge() > 0 && _dimuon->get_same_event())
    type = "++";

  return type;
}

