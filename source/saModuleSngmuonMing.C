// $Id: saModuleSngmuonMing.C,v 1.5 2013/08/17 05:41:42 mxliu Exp $                                                                                             
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
#include <DiMuonContainer.h>
#include <DiMuon.h>
#include <SingleMuonContainer.h>
#include <SingleMuon.h>
#include <RpcSingleMuonContainer.h>
#include <RpcSingleMuon.h>
#include <MCDiMuonContainer.h>
#include <MCDiMuon.h>

#include "saFlag.h"
#include "saModuleSngmuonMing.h"

using namespace std;

//-- Single Muon pT
saModuleSngmuonMing::saModuleSngmuonMing(const std::string &name) :
  saModuleBase(name),                                         // muons
  _h_pT(NULL),   _h_pT_n(NULL),   _h_pT_s(NULL),              //
  _h_pT_p(NULL), _h_pT_p_n(NULL), _h_pT_p_s(NULL),            //
  _h_pT_m(NULL), _h_pT_m_n(NULL), _h_pT_m_s(NULL),            //
  _h_pTH(NULL),  _h_pTH_n(NULL),  _h_pTH_s(NULL),             // hadrons
  _h_pTH_p(NULL),_h_pTH_p_n(NULL),_h_pTH_p_s(NULL),           //
  _h_pTH_m(NULL),_h_pTH_m_n(NULL),_h_pTH_m_s(NULL),           //
  _h_eta(NULL),    _h_eta_p(NULL),    _h_eta_m(NULL),         // heavy flavor
  _h_etaH(NULL),   _h_etaH_p(NULL),   _h_etaH_m(NULL),         // heavy flavor
  _h_etaW10(NULL), _h_etaW10_p(NULL), _h_etaW10_m(NULL),      // W->muon
  _h_etaW15(NULL), _h_etaW15_p(NULL), _h_etaW15_m(NULL),      // 
  _h_etaW20(NULL), _h_etaW20_p(NULL), _h_etaW20_m(NULL)
{

}

saModuleSngmuonMing::~saModuleSngmuonMing()
{

}

//! global initialization
int
saModuleSngmuonMing::init(PHCompositeNode *topNode, sa_hist_mangager_ptr hm)
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
          "SngmuonFlag", "PHObject");

      if (!node)
        {
          cout
              << "saModuleSngmuon::Init failed to create saEventProperty Node"
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
      cout << "saModuleSngmuonMing::Init failed to create saEventProperty"
          << endl;
      return ABORTRUN;
    }


  // --- define user hisotgrams ---

  //
  // -- general event and/or track quality histograms 
  //

  TH1F * h_ST1Z_N = new TH1F("SngMuon_ST1Z_N", "Single muon St1-Z; ", 400, 180, 200);
  _h_ST1Z_N = new saHist(
      //
      h_ST1Z_N, //the template histogram, any TH1 and derivatives is accepted
      saHist::EVENT_PROPERTY, // flags
      Verbosity() // verbosity
      );
  hm.get()->registerHisto(_h_ST1Z_N);

  TH1F * h_ST1Z_S = new TH1F("SngMuon_ST1Z_S", "Single muon St1-Z; ", 400, -200, -180);
  _h_ST1Z_S = new saHist(
      //
      h_ST1Z_S, //the template histogram, any TH1 and derivatives is accepted
      saHist::EVENT_PROPERTY, // flags
      Verbosity() // verbosity
      );
  hm.get()->registerHisto(_h_ST1Z_S);



  // ST-1 multiple scatering angle 
  TH1F * h_dA = new TH1F("SngMuon_dA", "Single muon dA; ", 100, -1, 1);
  _h_dA = new saHist(
      //
      h_dA, //the template histogram, any TH1 and derivatives is accepted
      saHist::EVENT_PROPERTY | saHist::FILL_PROPERTY, // flags
      Verbosity() // verbosity
      );
  hm.get()->registerHisto(_h_dA);

  TH2F * h2_dA = new TH2F("SngMuon_dA2", "Single muon dA vs pz; ", 200, -50, 50, 100, 0.0, 0.5);
  _h2_dA = new saHist(
      //
      h2_dA, //the template histogram, any TH1 and derivatives is accepted
      saHist::EVENT_PROPERTY | saHist::FILL_PROPERTY, // flags
      Verbosity() // verbosity
      );
  hm.get()->registerHisto(_h2_dA);

  // ST-1 <phi-1>
  TH1F * h_phi_1N = new TH1F("SngMuon_phi_1N", "North Single muon ST-1 Phi; ", 200, -4, 4);
  _h_phi_1N = new saHist(
      //
      h_phi_1N, //the template histogram, any TH1 and derivatives is accepted
      saHist::EVENT_PROPERTY | saHist::FILL_PROPERTY, // flags
      Verbosity() // verbosity
      );
  hm.get()->registerHisto(_h_phi_1N);

  TH1F * h_phi_1S = new TH1F("SngMuon_phi_1S", "South Single muon ST-1 Phi; ", 200, -4, 4);
  _h_phi_1S = new saHist(
      //
      h_phi_1S, //the template histogram, any TH1 and derivatives is accepted
      saHist::EVENT_PROPERTY | saHist::FILL_PROPERTY, // flags
      Verbosity() // verbosity
      );
  hm.get()->registerHisto(_h_phi_1S);



  // ST-2 <phi-2>
  TH1F * h_phi_2N = new TH1F("SngMuon_phi_2N", "North Single muon ST-2 Phi; ", 200, -4, 4);
  _h_phi_2N = new saHist(
      //
      h_phi_2N, //the template histogram, any TH1 and derivatives is accepted
      saHist::EVENT_PROPERTY | saHist::FILL_PROPERTY, // flags
      Verbosity() // verbosity
      );
  hm.get()->registerHisto(_h_phi_2N);

  TH1F * h_phi_2S = new TH1F("SngMuon_phi_2S", "South Single muon ST-2 Phi; ", 200, -4, 4);
  _h_phi_2S = new saHist(
      //
      h_phi_2S, //the template histogram, any TH1 and derivatives is accepted
      saHist::EVENT_PROPERTY | saHist::FILL_PROPERTY, // flags
      Verbosity() // verbosity
      );
  hm.get()->registerHisto(_h_phi_2S);


  // ST-3 <phi-3>
  TH1F * h_phi_3N = new TH1F("SngMuon_phi_3N", "North Single muon ST-3 Phi; ", 200, -4, 4);
  _h_phi_3N = new saHist(
      //
      h_phi_3N, //the template histogram, any TH1 and derivatives is accepted
      saHist::EVENT_PROPERTY | saHist::FILL_PROPERTY, // flags
      Verbosity() // verbosity
      );
  hm.get()->registerHisto(_h_phi_3N);

  TH1F * h_phi_3S = new TH1F("SngMuon_phi_3S", "North Single muon ST-3 Phi; ", 200, -4, 4);
  _h_phi_3S = new saHist(
      //
      h_phi_3S, //the template histogram, any TH1 and derivatives is accepted
      saHist::EVENT_PROPERTY | saHist::FILL_PROPERTY, // flags
      Verbosity() // verbosity
      );
  hm.get()->registerHisto(_h_phi_3S);



  // ST-1-3 <phi-1-3>
  TH1F * h_phi_13N = new TH1F("SngMuon_phi_13N", "North Single muon ST-13 dPhi; ", 100, -0.5, 0.5);
  _h_phi_13N = new saHist(
      //
      h_phi_13N, //the template histogram, any TH1 and derivatives is accepted
      saHist::EVENT_PROPERTY | saHist::FILL_PROPERTY, // flags
      Verbosity() // verbosity
      );
  hm.get()->registerHisto(_h_phi_13N);

  TH1F * h_phi_13S = new TH1F("SngMuon_phi_13S", "South Single muon ST-13 dPhi; ", 100, -0.5, 0.5);
  _h_phi_13S = new saHist(
      //
      h_phi_13S, //the template histogram, any TH1 and derivatives is accepted
      saHist::EVENT_PROPERTY | saHist::FILL_PROPERTY, // flags
      Verbosity() // verbosity
      );
  hm.get()->registerHisto(_h_phi_13S);


  // ST-2-3 <phi-2-3>
  TH1F * h_phi_23N = new TH1F("SngMuon_phi_23N", "North Single muon ST-23 dPhi; ", 100, -0.5, 0.5);
  _h_phi_23N = new saHist(
      //
      h_phi_23N, //the template histogram, any TH1 and derivatives is accepted
      saHist::EVENT_PROPERTY | saHist::FILL_PROPERTY, // flags
      Verbosity() // verbosity
      );
  hm.get()->registerHisto(_h_phi_23N);

  TH1F * h_phi_23S = new TH1F("SngMuon_phi_23S", "South Single muon ST-23 dPhi; ", 100, -0.5, 0.5);
  _h_phi_23S = new saHist(
      //
      h_phi_23S, //the template histogram, any TH1 and derivatives is accepted
      saHist::EVENT_PROPERTY | saHist::FILL_PROPERTY, // flags
      Verbosity() // verbosity
      );
  hm.get()->registerHisto(_h_phi_23S);


  // -- 2D plots vs pZ ST-1-3
  TH2F * h2_phi_13N = new TH2F("SngMuon_phi_13N2", "North Single muon pZ*dPhi-13 vs pZ; ", 200, -50, 50,  100, 0, 5);
  _h2_phi_13N = new saHist(
      //
      h2_phi_13N, //the template histogram, any TH1 and derivatives is accepted
      saHist::EVENT_PROPERTY | saHist::FILL_PROPERTY, // flags
      Verbosity() // verbosity
      );
  hm.get()->registerHisto(_h2_phi_13N);

  TH2F * h2_phi_13S = new TH2F("SngMuon_phi_13S2", "South Single muon pZ*dPhi-13  vs pZ; ", 200, -50, 50,  100, 0, 5);
  _h2_phi_13S = new saHist(
      //
      h2_phi_13S, //the template histogram, any TH1 and derivatives is accepted
      saHist::EVENT_PROPERTY | saHist::FILL_PROPERTY, // flags
      Verbosity() // verbosity
      );
  hm.get()->registerHisto(_h2_phi_13S);

  // -- 2D plots vs pZ ST-2-3
  TH2F * h2_phi_23N = new TH2F("SngMuon_phi_23N2", "North Single muon pZ*dPhi-23 vs pZ; ", 200, -50, 50,  100, 0, 5);
  _h2_phi_23N = new saHist(
      //
      h2_phi_23N, //the template histogram, any TH1 and derivatives is accepted
      saHist::EVENT_PROPERTY | saHist::FILL_PROPERTY, // flags
      Verbosity() // verbosity
      );
  hm.get()->registerHisto(_h2_phi_23N);

  TH2F * h2_phi_23S = new TH2F("SngMuon_phi_23S2", "South Single muon pZ*dPhi-23 vs pZ; ", 200, -50, 50,  100, 0, 5);
  _h2_phi_23S = new saHist(
      //
      h2_phi_23S, //the template histogram, any TH1 and derivatives is accepted
      saHist::EVENT_PROPERTY | saHist::FILL_PROPERTY, // flags
      Verbosity() // verbosity
      );
  hm.get()->registerHisto(_h2_phi_23S);




  // -- RPC1
  TH1F * h_Rpc1DCA_N = new TH1F("SngMuonRPC1_DCA_N", "Single muon North RPC1 DCA; ", 150, -5, 10);
  _h_Rpc1DCA_N = new saHist(
      //
      h_Rpc1DCA_N, //the template histogram, any TH1 and derivatives is accepted
      saHist::EVENT_PROPERTY | saHist::FILL_PROPERTY, // flags
      Verbosity() // verbosity
      );
  hm.get()->registerHisto(_h_Rpc1DCA_N);

  TH1F * h_Rpc1DCA_S = new TH1F("SngMuonRPC1_DCA_S", "Single muon South RPC1 DCA; ", 150, -5, 10);
  _h_Rpc1DCA_S = new saHist(
      //
      h_Rpc1DCA_S, //the template histogram, any TH1 and derivatives is accepted
      saHist::EVENT_PROPERTY | saHist::FILL_PROPERTY, // flags
      Verbosity() // verbosity
      );
  hm.get()->registerHisto(_h_Rpc1DCA_S);

  // --- RPC3 
  TH1F * h_Rpc3DCA_N = new TH1F("SngMuonRPC3_DCA_N", "Single muon North RPC3 DCA; ", 100, -10, 100);
  _h_Rpc3DCA_N = new saHist(
      //
      h_Rpc3DCA_N, //the template histogram, any TH1 and derivatives is accepted
      saHist::EVENT_PROPERTY | saHist::FILL_PROPERTY, // flags
      Verbosity() // verbosity
      );
  hm.get()->registerHisto(_h_Rpc3DCA_N);

  TH1F * h_Rpc3DCA_S = new TH1F("SngMuonRPC3_DCA_S", "Single muon South RPC3 DCA; ", 100, -10, 100);
  _h_Rpc3DCA_S = new saHist(
      //
      h_Rpc3DCA_S, //the template histogram, any TH1 and derivatives is accepted
      saHist::EVENT_PROPERTY | saHist::FILL_PROPERTY, // flags
      Verbosity() // verbosity
      );
  hm.get()->registerHisto(_h_Rpc3DCA_S);


  //
  // --- for comparision with run11/12 style W->muon analysis 
  //

  TH1F * h_dw23N_p = new TH1F("SngMuon_dw23N_p", "North (+) SingleMuon dw23", 100, -1, 1);
  _h_dw23N_p = new saHist(
      //
      h_dw23N_p, //the template histogram, any TH1 and derivatives is accepted
      saHist::EVENT_PROPERTY | saHist::FILL_PROPERTY, // flags
      Verbosity() // verbosity
      );
  hm.get()->registerHisto(_h_dw23N_p);

  TH1F * h_dw23N_m = new TH1F("SngMuon_dw23N_m", "North (-) SingleMuon dw23", 100, -1, 1);
  _h_dw23N_m = new saHist(
      //
      h_dw23N_m, //the template histogram, any TH1 and derivatives is accepted
      saHist::EVENT_PROPERTY | saHist::FILL_PROPERTY, // flags
      Verbosity() // verbosity
      );
  hm.get()->registerHisto(_h_dw23N_m);



  TH1F * h_dw23S_p = new TH1F("SngMuon_dw23S_p", "South (+) SingleMuon dw23", 100, -1, 1);
  _h_dw23S_p = new saHist(
      //
      h_dw23S_p, //the template histogram, any TH1 and derivatives is accepted
      saHist::EVENT_PROPERTY | saHist::FILL_PROPERTY, // flags
      Verbosity() // verbosity
      );
  hm.get()->registerHisto(_h_dw23S_p);

  TH1F * h_dw23S_m = new TH1F("SngMuon_dw23S_m", "South (-) SingleMuon dw23", 100, -1, 1);
  _h_dw23S_m = new saHist(
      //
      h_dw23S_m, //the template histogram, any TH1 and derivatives is accepted
      saHist::EVENT_PROPERTY | saHist::FILL_PROPERTY, // flags
      Verbosity() // verbosity
      );
  hm.get()->registerHisto(_h_dw23S_m);

  // pT > 10 

  TH1F * h_dw23N_PT10_p = new TH1F("SngMuon_dw23N_PT10_p", "North (+) pT> 10 SingleMuon dw23", 100, -1, 1);
  _h_dw23N_PT10_p = new saHist(
      //
      h_dw23N_PT10_p, //the template histogram, any TH1 and derivatives is accepted
      saHist::EVENT_PROPERTY | saHist::FILL_PROPERTY, // flags
      Verbosity() // verbosity
      );
  hm.get()->registerHisto(_h_dw23N_PT10_p);

  TH1F * h_dw23N_PT10_m = new TH1F("SngMuon_dw23N_PT10_m", "North (-) pT> 10 SingleMuon dw23", 100, -1, 1);
  _h_dw23N_PT10_m = new saHist(
      //
      h_dw23N_PT10_m, //the template histogram, any TH1 and derivatives is accepted
      saHist::EVENT_PROPERTY | saHist::FILL_PROPERTY, // flags
      Verbosity() // verbosity
      );
  hm.get()->registerHisto(_h_dw23N_PT10_m);



  TH1F * h_dw23S_PT10_p = new TH1F("SngMuon_dw23S_PT10_p", "South (+) pT> 10 SingleMuon dw23", 100, -1, 1);
  _h_dw23S_PT10_p = new saHist(
      //
      h_dw23S_PT10_p, //the template histogram, any TH1 and derivatives is accepted
      saHist::EVENT_PROPERTY | saHist::FILL_PROPERTY, // flags
      Verbosity() // verbosity
      );
  hm.get()->registerHisto(_h_dw23S_PT10_p);

  TH1F * h_dw23S_PT10_m = new TH1F("SngMuon_dw23S_PT10_m", "South (-) pT>10 SingleMuon dw23", 100, -1, 1);
  _h_dw23S_PT10_m = new saHist(
      //
      h_dw23S_PT10_m, //the template histogram, any TH1 and derivatives is accepted
      saHist::EVENT_PROPERTY | saHist::FILL_PROPERTY, // flags
      Verbosity() // verbosity
      );
  hm.get()->registerHisto(_h_dw23S_PT10_m);

  // pT > 15 

  TH1F * h_dw23N_PT15_p = new TH1F("SngMuon_dw23N_PT15_p", "North (+) pT> 10 SingleMuon dw23", 100, -1, 1);
  _h_dw23N_PT15_p = new saHist(
      //
      h_dw23N_PT15_p, //the template histogram, any TH1 and derivatives is accepted
      saHist::EVENT_PROPERTY | saHist::FILL_PROPERTY, // flags
      Verbosity() // verbosity
      );
  hm.get()->registerHisto(_h_dw23N_PT15_p);

  TH1F * h_dw23N_PT15_m = new TH1F("SngMuon_dw23N_PT15_m", "North (-) pT> 10 SingleMuon dw23", 100, -1, 1);
  _h_dw23N_PT15_m = new saHist(
      //
      h_dw23N_PT15_m, //the template histogram, any TH1 and derivatives is accepted
      saHist::EVENT_PROPERTY | saHist::FILL_PROPERTY, // flags
      Verbosity() // verbosity
      );
  hm.get()->registerHisto(_h_dw23N_PT15_m);



  TH1F * h_dw23S_PT15_p = new TH1F("SngMuon_dw23S_PT15_p", "South (+) pT> 10 SingleMuon dw23", 100, -1, 1);
  _h_dw23S_PT15_p = new saHist(
      //
      h_dw23S_PT15_p, //the template histogram, any TH1 and derivatives is accepted
      saHist::EVENT_PROPERTY | saHist::FILL_PROPERTY, // flags
      Verbosity() // verbosity
      );
  hm.get()->registerHisto(_h_dw23S_PT15_p);

  TH1F * h_dw23S_PT15_m = new TH1F("SngMuon_dw23S_PT15_m", "South (-) pT>10 SingleMuon dw23", 100, -1, 1);
  _h_dw23S_PT15_m = new saHist(
      //
      h_dw23S_PT15_m, //the template histogram, any TH1 and derivatives is accepted
      saHist::EVENT_PROPERTY | saHist::FILL_PROPERTY, // flags
      Verbosity() // verbosity
      );
  hm.get()->registerHisto(_h_dw23S_PT15_m);

  // pT > 20 

  TH1F * h_dw23N_PT20_p = new TH1F("SngMuon_dw23N_PT20_p", "North (+) pT> 10 SingleMuon dw23", 100, -1, 1);
  _h_dw23N_PT20_p = new saHist(
      //
      h_dw23N_PT20_p, //the template histogram, any TH1 and derivatives is accepted
      saHist::EVENT_PROPERTY | saHist::FILL_PROPERTY, // flags
      Verbosity() // verbosity
      );
  hm.get()->registerHisto(_h_dw23N_PT20_p);

  TH1F * h_dw23N_PT20_m = new TH1F("SngMuon_dw23N_PT20_m", "North (-) pT> 10 SingleMuon dw23", 100, -1, 1);
  _h_dw23N_PT20_m = new saHist(
      //
      h_dw23N_PT20_m, //the template histogram, any TH1 and derivatives is accepted
      saHist::EVENT_PROPERTY | saHist::FILL_PROPERTY, // flags
      Verbosity() // verbosity
      );
  hm.get()->registerHisto(_h_dw23N_PT20_m);



  TH1F * h_dw23S_PT20_p = new TH1F("SngMuon_dw23S_PT20_p", "South (+) pT> 10 SingleMuon dw23", 100, -1, 1);
  _h_dw23S_PT20_p = new saHist(
      //
      h_dw23S_PT20_p, //the template histogram, any TH1 and derivatives is accepted
      saHist::EVENT_PROPERTY | saHist::FILL_PROPERTY, // flags
      Verbosity() // verbosity
      );
  hm.get()->registerHisto(_h_dw23S_PT20_p);

  TH1F * h_dw23S_PT20_m = new TH1F("SngMuon_dw23S_PT20_m", "South (-) pT>10 SingleMuon dw23", 100, -1, 1);
  _h_dw23S_PT20_m = new saHist(
      //
      h_dw23S_PT20_m, //the template histogram, any TH1 and derivatives is accepted
      saHist::EVENT_PROPERTY | saHist::FILL_PROPERTY, // flags
      Verbosity() // verbosity
      );
  hm.get()->registerHisto(_h_dw23S_PT20_m);




  //
  // --- histograms for spin asymmetry study ----
  //

  //-- Single Muon pT (+,-) 
  TH1F * h_pT = new TH1F("SngMuonPT", "Single muon pT; pT (GeV)", 20, 1, 10);
  _h_pT = new saHist(
      //
      h_pT, //the template histogram, any TH1 and derivatives is accepted
      saHist::DEFAULT_FLAG | saHist::RUN_PROPERTY | saHist::FILL_PROPERTY
          | saHist::CROSSING_CHECKS, // flags
      Verbosity() // verbosity
      );
  hm.get()->registerHisto(_h_pT);

  // --- north ---
  TH1F * h_pT_n = new TH1F("SngMuonPTN", "North Single Muon pT; pT (GeV)", 20, 1, 10);
  _h_pT_n = new saHist(
      //
      h_pT_n, //the template histogram, any TH1 and derivatives is accepted
      saHist::DEFAULT_FLAG | saHist::RUN_PROPERTY | saHist::FILL_PROPERTY
          | saHist::CROSSING_CHECKS, // flags
      Verbosity() // verbosity
      );
  hm.get()->registerHisto(_h_pT_n);

  // --- south ---
  TH1F * h_pT_s = new TH1F("SngMuonPTS", "South Single Muon pT; pT (GeV)", 20, 1, 10);
  _h_pT_s = new saHist(
      //
      h_pT_s, //the template histogram, any TH1 and derivatives is accepted
      saHist::DEFAULT_FLAG | saHist::RUN_PROPERTY | saHist::FILL_PROPERTY
          | saHist::CROSSING_CHECKS, // flags
      Verbosity() // verbosity
      );
  hm.get()->registerHisto(_h_pT_s);


  //-- Single Muon pT (+)
  TH1F * h_pT_p = new TH1F("SngMuonPTp", "Single muon (+) pT; pT (GeV)", 20, 1, 10);
  _h_pT_p = new saHist(
      //
      h_pT_p, //the template histogram, any TH1 and derivatives is accepted
      saHist::DEFAULT_FLAG | saHist::RUN_PROPERTY | saHist::FILL_PROPERTY
          | saHist::CROSSING_CHECKS, // flags
      Verbosity() // verbosity
      );
  hm.get()->registerHisto(_h_pT_p);

  // --- north ---
  TH1F * h_pT_p_n = new TH1F("SngMuonPTpN", "North Single Muon (+)pT; pT (GeV)", 20, 1, 10);
  _h_pT_p_n = new saHist(
      //
      h_pT_p_n, //the template histogram, any TH1 and derivatives is accepted
      saHist::DEFAULT_FLAG | saHist::RUN_PROPERTY | saHist::FILL_PROPERTY
          | saHist::CROSSING_CHECKS, // flags
      Verbosity() // verbosity
      );
  hm.get()->registerHisto(_h_pT_p_n);

  // --- south ---
  TH1F * h_pT_p_s = new TH1F("SngMuonPTpS", "South Single Muon (+)pT; pT (GeV)", 20, 1, 10);
  _h_pT_p_s = new saHist(
      //
      h_pT_p_s, //the template histogram, any TH1 and derivatives is accepted
      saHist::DEFAULT_FLAG | saHist::RUN_PROPERTY | saHist::FILL_PROPERTY
          | saHist::CROSSING_CHECKS, // flags
      Verbosity() // verbosity
      );
  hm.get()->registerHisto(_h_pT_p_s);


  //-- Single Muon pT (-)
  TH1F * h_pT_m = new TH1F("SngMuonPTm", "Single muon (-) pT; pT (GeV)", 20, 1, 10);
  _h_pT_m = new saHist(
      //
      h_pT_m, //the template histogram, any TH1 and derivatives is accepted
      saHist::DEFAULT_FLAG | saHist::RUN_PROPERTY | saHist::FILL_PROPERTY
          | saHist::CROSSING_CHECKS, // flags
      Verbosity() // verbosity
      );
  hm.get()->registerHisto(_h_pT_m);

  // --- north ---
  TH1F * h_pT_m_n = new TH1F("SngMuonPTmN", "North Single Muon (-)pT; pT (GeV)", 20, 1, 10);
  _h_pT_m_n = new saHist(
      //
      h_pT_m_n, //the template histogram, any TH1 and derivatives is accepted
      saHist::DEFAULT_FLAG | saHist::RUN_PROPERTY | saHist::FILL_PROPERTY
          | saHist::CROSSING_CHECKS, // flags
      Verbosity() // verbosity
      );
  hm.get()->registerHisto(_h_pT_m_n);

  // --- south ---
  TH1F * h_pT_m_s = new TH1F("SngMuonPTmS", "South Single Muon (-)pT; pT (GeV)", 20, 1, 10);
  _h_pT_m_s = new saHist(
      //
      h_pT_m_s, //the template histogram, any TH1 and derivatives is accepted
      saHist::DEFAULT_FLAG | saHist::RUN_PROPERTY | saHist::FILL_PROPERTY
          | saHist::CROSSING_CHECKS, // flags
      Verbosity() // verbosity
      );
  hm.get()->registerHisto(_h_pT_m_s);


  //
  //--- stopped hadrons -------------------------------
  //


  //-- Stoped Hadron pT (+,-) 
  TH1F * h_pTH = new TH1F("SngHadronPT", "Single hadron pT; pT (GeV)", 20, 1, 10);
  _h_pTH = new saHist(
      //
      h_pTH, //the template histogram, any TH1 and derivatives is accepted
      saHist::DEFAULT_FLAG | saHist::RUN_PROPERTY | saHist::FILL_PROPERTY
          | saHist::CROSSING_CHECKS, // flags
      Verbosity() // verbosity
      );
  hm.get()->registerHisto(_h_pTH);

  // --- north ---
  TH1F * h_pTH_n = new TH1F("SngHadronPTN", "North Single Hadron pT; pT (GeV)", 20, 1, 10);
  _h_pTH_n = new saHist(
      //
      h_pTH_n, //the template histogram, any TH1 and derivatives is accepted
      saHist::DEFAULT_FLAG | saHist::RUN_PROPERTY | saHist::FILL_PROPERTY
          | saHist::CROSSING_CHECKS, // flags
      Verbosity() // verbosity
      );
  hm.get()->registerHisto(_h_pTH_n);

  // --- south ---
  TH1F * h_pTH_s = new TH1F("SngHadronPTS", "South Single Hadron pT; pT (GeV)", 20, 1, 10);
  _h_pTH_s = new saHist(
      //
      h_pTH_s, //the template histogram, any TH1 and derivatives is accepted
      saHist::DEFAULT_FLAG | saHist::RUN_PROPERTY | saHist::FILL_PROPERTY
          | saHist::CROSSING_CHECKS, // flags
      Verbosity() // verbosity
      );
  hm.get()->registerHisto(_h_pTH_s);


  //-- Single Hadron pT (+)
  TH1F * h_pTH_p = new TH1F("SngHadronPTp", "Single Hadron (+) pT; pT (GeV)", 20, 1, 10);
  _h_pTH_p = new saHist(
      //
      h_pTH_p, //the template histogram, any TH1 and derivatives is accepted
      saHist::DEFAULT_FLAG | saHist::RUN_PROPERTY | saHist::FILL_PROPERTY
          | saHist::CROSSING_CHECKS, // flags
      Verbosity() // verbosity
      );
  hm.get()->registerHisto(_h_pTH_p);

  // --- north ---
  TH1F * h_pTH_p_n = new TH1F("SngHadronPTpN", "North Single Hadron (+)pT; pT (GeV)", 20, 1, 10);
  _h_pTH_p_n = new saHist(
      //
      h_pTH_p_n, //the template histogram, any TH1 and derivatives is accepted
      saHist::DEFAULT_FLAG | saHist::RUN_PROPERTY | saHist::FILL_PROPERTY
          | saHist::CROSSING_CHECKS, // flags
      Verbosity() // verbosity
      );
  hm.get()->registerHisto(_h_pTH_p_n);

  // --- south ---
  TH1F * h_pTH_p_s = new TH1F("SngHadronPTpS", "South Single Hadron (+)pT; pT (GeV)", 20, 1, 10);
  _h_pTH_p_s = new saHist(
      //
      h_pTH_p_s, //the template histogram, any TH1 and derivatives is accepted
      saHist::DEFAULT_FLAG | saHist::RUN_PROPERTY | saHist::FILL_PROPERTY
          | saHist::CROSSING_CHECKS, // flags
      Verbosity() // verbosity
      );
  hm.get()->registerHisto(_h_pTH_p_s);


  //-- Single Hadron pT (-)
  TH1F * h_pTH_m = new TH1F("SngHadronPTm", "Single Hadron (-) pT; pT (GeV)", 20, 1, 10);
  _h_pTH_m = new saHist(
      //
      h_pTH_m, //the template histogram, any TH1 and derivatives is accepted
      saHist::DEFAULT_FLAG | saHist::RUN_PROPERTY | saHist::FILL_PROPERTY
          | saHist::CROSSING_CHECKS, // flags
      Verbosity() // verbosity
      );
  hm.get()->registerHisto(_h_pTH_m);

  // --- north ---
  TH1F * h_pTH_m_n = new TH1F("SngHadronPTmN", "North Single Hadron (-)pT; pT (GeV)", 20, 1, 10);
  _h_pTH_m_n = new saHist(
      //
      h_pTH_m_n, //the template histogram, any TH1 and derivatives is accepted
      saHist::DEFAULT_FLAG | saHist::RUN_PROPERTY | saHist::FILL_PROPERTY
          | saHist::CROSSING_CHECKS, // flags
      Verbosity() // verbosity
      );
  hm.get()->registerHisto(_h_pTH_m_n);

  // --- south ---
  TH1F * h_pTH_m_s = new TH1F("SngHadronPTmS", "South Single Hadron (-)pT; pT (GeV)", 20, 1, 10);
  _h_pTH_m_s = new saHist(
      //
      h_pTH_m_s, //the template histogram, any TH1 and derivatives is accepted
      saHist::DEFAULT_FLAG | saHist::RUN_PROPERTY | saHist::FILL_PROPERTY
          | saHist::CROSSING_CHECKS, // flags
      Verbosity() // verbosity
      );
  hm.get()->registerHisto(_h_pTH_m_s);


  //
  // -------  single muon  rapidity --------
  //

  TH1F * h_eta = new TH1F("SngMuonEta", "Single Muon (+,-) Eta;  ", 12, -3, 3);
  _h_eta = new saHist(
      //
      h_eta, //the template histogram, any TH1 and derivatives is accepted
      saHist::DEFAULT_FLAG | saHist::RUN_PROPERTY | saHist::FILL_PROPERTY
          | saHist::CROSSING_CHECKS, // flags
      Verbosity() // verbosity
      );
  hm.get()->registerHisto(_h_eta);


  TH1F * h_eta_p = new TH1F("SngMuonEtap", "Single Muon (+) Eta; ", 12, -3, 3);
  _h_eta_p = new saHist(
      //
      h_eta_p, //the template histogram, any TH1 and derivatives is accepted
      saHist::DEFAULT_FLAG | saHist::RUN_PROPERTY | saHist::FILL_PROPERTY
          | saHist::CROSSING_CHECKS, // flags
      Verbosity() // verbosity
      );
  hm.get()->registerHisto(_h_eta_p);

  TH1F * h_eta_m = new TH1F("SngMuonEtam", "Single Muon (-) Eta;  ", 12, -3, 3);
  _h_eta_m = new saHist(
      //
      h_eta_m, //the template histogram, any TH1 and derivatives is accepted
      saHist::DEFAULT_FLAG | saHist::RUN_PROPERTY | saHist::FILL_PROPERTY
          | saHist::CROSSING_CHECKS, // flags
      Verbosity() // verbosity
      );
  hm.get()->registerHisto(_h_eta_m);

  //
  // ----  rapidity for W->muon pT > 10 GeV ---
  //

  TH1F * h_etaW10 = new TH1F("SngMuonEtaW10", "Single Muon (+,-) Eta, pT>10;  ", 12, -3, 3);
  _h_etaW10 = new saHist(
      //
      h_etaW10, //the template histogram, any TH1 and derivatives is accepted
      saHist::DEFAULT_FLAG | saHist::RUN_PROPERTY | saHist::FILL_PROPERTY
          | saHist::CROSSING_CHECKS, // flags
      Verbosity() // verbosity
      );
  hm.get()->registerHisto(_h_etaW10);


  TH1F * h_etaW10_p = new TH1F("SngMuonEtaW10p", "Single Muon (+) Eta; pT>10", 12, -3, 3);
  _h_etaW10_p = new saHist(
      //
      h_etaW10_p, //the template histogram, any TH1 and derivatives is accepted
      saHist::DEFAULT_FLAG | saHist::RUN_PROPERTY | saHist::FILL_PROPERTY
          | saHist::CROSSING_CHECKS, // flags
      Verbosity() // verbosity
      );
  hm.get()->registerHisto(_h_etaW10_p);

  TH1F * h_etaW10_m = new TH1F("SngMuonEtaW10m", "Single Muon (-) Eta; pT>10 ", 12, -3, 3);
  _h_etaW10_m = new saHist(
      //
      h_etaW10_m, //the template histogram, any TH1 and derivatives is accepted
      saHist::DEFAULT_FLAG | saHist::RUN_PROPERTY | saHist::FILL_PROPERTY
          | saHist::CROSSING_CHECKS, // flags
      Verbosity() // verbosity
      );
  hm.get()->registerHisto(_h_etaW10_m);

  // --- W->muon rapidity , pT>15 
  TH1F * h_etaW15 = new TH1F("SngMuonEtaW15", "Single Muon (+,-) Eta, pT>15;  ", 12, -3, 3);
  _h_etaW15 = new saHist(
      //
      h_etaW15, //the template histogram, any TH1 and derivatives is accepted
      saHist::DEFAULT_FLAG | saHist::RUN_PROPERTY | saHist::FILL_PROPERTY
          | saHist::CROSSING_CHECKS, // flags
      Verbosity() // verbosity
      );
  hm.get()->registerHisto(_h_etaW15);


  TH1F * h_etaW15_p = new TH1F("SngMuonEtaW15p", "Single Muon (+) Eta; pT>15", 12, -3, 3);
  _h_etaW15_p = new saHist(
      //
      h_etaW15_p, //the template histogram, any TH1 and derivatives is accepted
      saHist::DEFAULT_FLAG | saHist::RUN_PROPERTY | saHist::FILL_PROPERTY
          | saHist::CROSSING_CHECKS, // flags
      Verbosity() // verbosity
      );
  hm.get()->registerHisto(_h_etaW15_p);

  TH1F * h_etaW15_m = new TH1F("SngMuonEtaW15m", "Single Muon (-) Eta; pT>15 ", 12, -3, 3);
  _h_etaW15_m = new saHist(
      //
      h_etaW15_m, //the template histogram, any TH1 and derivatives is accepted
      saHist::DEFAULT_FLAG | saHist::RUN_PROPERTY | saHist::FILL_PROPERTY
          | saHist::CROSSING_CHECKS, // flags
      Verbosity() // verbosity
      );
  hm.get()->registerHisto(_h_etaW15_m);

  // --- W->muon rapidity , pT>20 
  TH1F * h_etaW20 = new TH1F("SngMuonEtaW20", "Single Muon (+,-) Eta, pT>20;  ", 12, -3, 3);
  _h_etaW20 = new saHist(
      //
      h_etaW20, //the template histogram, any TH1 and derivatives is accepted
      saHist::DEFAULT_FLAG | saHist::RUN_PROPERTY | saHist::FILL_PROPERTY
          | saHist::CROSSING_CHECKS, // flags
      Verbosity() // verbosity
      );
  hm.get()->registerHisto(_h_etaW20);


  TH1F * h_etaW20_p = new TH1F("SngMuonEtaW20p", "Single Muon (+) Eta; pT>20", 12, -3, 3);
  _h_etaW20_p = new saHist(
      //
      h_etaW20_p, //the template histogram, any TH1 and derivatives is accepted
      saHist::DEFAULT_FLAG | saHist::RUN_PROPERTY | saHist::FILL_PROPERTY
          | saHist::CROSSING_CHECKS, // flags
      Verbosity() // verbosity
      );
  hm.get()->registerHisto(_h_etaW20_p);

  TH1F * h_etaW20_m = new TH1F("SngMuonEtaW20m", "Single Muon (-) Eta; pT>20 ", 12, -3, 3);
  _h_etaW20_m = new saHist(
      //
      h_etaW20_m, //the template histogram, any TH1 and derivatives is accepted
      saHist::DEFAULT_FLAG | saHist::RUN_PROPERTY | saHist::FILL_PROPERTY
          | saHist::CROSSING_CHECKS, // flags
      Verbosity() // verbosity
      );
  hm.get()->registerHisto(_h_etaW20_m);

  //
  // --- stopped hadron rapidity ---
  //

  TH1F * h_etaH = new TH1F("SngHadronEta", "Single Hadron (+,-) Eta;  ", 12, -3, 3);
  _h_etaH = new saHist(
      //
      h_etaH, //the template histogram, any TH1 and derivatives is accepted
      saHist::DEFAULT_FLAG | saHist::RUN_PROPERTY | saHist::FILL_PROPERTY
          | saHist::CROSSING_CHECKS, // flags
      Verbosity() // verbosity
      );
  hm.get()->registerHisto(_h_etaH);


  TH1F * h_etaH_p = new TH1F("SngHadronEtap", "Single Hadron (+) Eta; ", 12, -3, 3);
  _h_etaH_p = new saHist(
      //
      h_etaH_p, //the template histogram, any TH1 and derivatives is accepted
      saHist::DEFAULT_FLAG | saHist::RUN_PROPERTY | saHist::FILL_PROPERTY
          | saHist::CROSSING_CHECKS, // flags
      Verbosity() // verbosity
      );
  hm.get()->registerHisto(_h_etaH_p);

  TH1F * h_etaH_m = new TH1F("SngHadronEtam", "Single Hadron (-) Eta;  ", 12, -3, 3);
  _h_etaH_m = new saHist(
      //
      h_etaH_m, //the template histogram, any TH1 and derivatives is accepted
      saHist::DEFAULT_FLAG | saHist::RUN_PROPERTY | saHist::FILL_PROPERTY 
          | saHist::CROSSING_CHECKS, // flags
      Verbosity() // verbosity
      );
  hm.get()->registerHisto(_h_etaH_m);



  return EVENT_OK;
}

//! Run initialization
int
saModuleSngmuonMing::init_run(PHCompositeNode *topNode,
    sa_hist_mangager_ptr hm)
{

  return EVENT_OK;
}

//! event method
int
saModuleSngmuonMing::event(PHCompositeNode *topNode, sa_hist_mangager_ptr hm)
{

  SingleMuonContainer* sngmuons = findNode::getClass<SingleMuonContainer>(topNode,
      "SingleMuonContainer");
  if (!sngmuons)
    {
      cout << "saModuleSngmuonMing:: SngMuonContainer - ERROR - not in Node Tree"
          << endl;
      return ABORTRUN;
    }

  const saEventProperty * ep = findNode::getClass<saEventProperty>(topNode,
      "saEventProperty");
  if (!ep)
    {
      cout << "saModuleSngmuonMing::event - ERROR  - Cannot find EventProperty node in Top Node"
          << endl;
      return ABORTRUN;
    }

  //  int fill_number = ep->get_fill_number();
  int run_number  = ep->get_run_number();

  // good run12pp510 run # = 366283 - 366059
  // good run13pp510 run # = 386773 - 397990

  bool run12pp510 = false;
  bool run13pp510 = false;

  if (run_number > 360000 && run_number < 370000) run12pp510 = true;
  if (run_number > 380000 && run_number < 400000) run13pp510 = true;


  RpcSingleMuonContainer* rpcsngmuons = findNode::getClass<RpcSingleMuonContainer>(topNode,
      "RpcSingleMuonContainer");

  if (!rpcsngmuons)
    {
      if(run13pp510){
	cout << "saModuleSngmuonMing:: RpcSingelMuonContainer - ERROR - not in Node Tree"
	     << endl;
	return ABORTRUN;
      }
      else {
	cout << "this is <run12 pp510>, NO RpcSingleMuonContainer !" << endl;
      }
    }

  if (Verbosity() >= 2)
    {
      cout << "saModuleSngmuonMing::event - INFO  - ";
      ep->identify();
    }


  // --- 
  saFlagC *flag = findNode::getClass<saFlagC>(topNode, "SngmuonFlag");
  if (!flag)
    {
      cout << "saModuleSngmuonMing::event - SngmuonFlag not in Node Tree"
          << endl;
      return ABORTRUN;
    }

  // syncronizely set the cut flag with the SingleMuonContainer
  flag->flags.Set(sngmuons->get_nSingleMuons());
  flag->flags.Reset(0);

  // save the event for DST files?
  //  int return_code = dimuons->get_nDiMuons() >0 ? EVENT_OK : DISCARDEVENT;

  // determine this event's primary VTX      
  double Evt_fvtxX = sngmuons->get_Evt_fvtxX();
  double Evt_fvtxY = sngmuons->get_Evt_fvtxY();
  double Evt_fvtxZ = sngmuons->get_Evt_fvtxZ();

  double Evt_fvtxX2 = sngmuons->get_Evt_fvtxX2();
  double Evt_fvtxY2 = sngmuons->get_Evt_fvtxY2();
  double Evt_fvtxZ2 = sngmuons->get_Evt_fvtxZ2();

  double Evt_bbcZ = sngmuons->get_Evt_fvtxZ();

  double Evt_X0 =  -999;   // the primary z-vtx associated with muon candidate
  double Evt_Y0 =  -999;   // the primary z-vtx associated with muon candidate
  double Evt_Z0 =  -999;   // the primary z-vtx associated with muon candidate
    
  if ( fabs(Evt_bbcZ) > 100) {

    return DISCARDEVENT;
  }

  if ( fabs(Evt_fvtxZ) > 100 && fabs(Evt_fvtxZ2) > 100 ) {

    return DISCARDEVENT;
  }

  if ( sqrt(Evt_fvtxX*Evt_fvtxX + Evt_fvtxY*Evt_fvtxY) > 5 && sqrt(Evt_fvtxX2*Evt_fvtxX2 + Evt_fvtxY2*Evt_fvtxY2) > 5 ) {
    return DISCARDEVENT;
  }
 

  bool fill_ok = false;  // output this event to ooutput pdst file if "true"

  for (unsigned int imuon = 0; imuon < sngmuons->get_nSingleMuons(); imuon++)
    {
      SingleMuon* sngmuon = sngmuons->get_SingleMuon(imuon);

      // get RPC information for Run13 pp510
      RpcSingleMuon* rpcsngmuon = NULL;
      if (run13pp510) {
        rpcsngmuon = rpcsngmuons->get_RpcSingleMuon(imuon);
      }

      int charge = 0;

      charge = 2*sngmuon->get_charge() -1 ;   // charge in pdst = 1 for Q+ and 0 for Q-

      int arm = (sngmuon->get_pz() > 0) ? 1 : 0;
     
      double pX = sngmuon->get_px();
      double pY = sngmuon->get_py();
      double pZ = sngmuon->get_pz();

      double pT = sqrt(pX*pX +  pY*pY); 

      double pTot =  sqrt(pT*pT + pZ*pZ);

      double pX1 = sngmuon->get_st1px();
      double pY1 = sngmuon->get_st1py();
      double pZ1 = sngmuon->get_st1pz();

      double X1 = sngmuon->get_xst1();
      double Y1 = sngmuon->get_yst1();
      //      double Z1 = sngmuon->get_zst1();
      double Z1 = 0;

      if (run13pp510){
	Z1 = rpcsngmuon->get_zst1();
      }
      else {
	if (pZ >0){
	  Z1 = 189.1;   // <Z1> = 189.1 cm in run13pp510 data
	}
	else {
	  Z1 =  -189.1;
	}
      }

      double X2 = sngmuon->get_xst2();
      double Y2 = sngmuon->get_yst2();

      double X3 = sngmuon->get_xst3();
      double Y3 = sngmuon->get_yst3();

      double phi_1 = atan2(Y1,X1);
      double phi_2 = atan2(Y2,X2);
      double phi_3 = atan2(Y3,X3);

      double phi_13 = phi_3 - phi_1;
      double phi_23 = phi_3 - phi_2;

      //      double phi_13 = acos( (X1*X3 + Y1*Y3)/sqrt(X1*X1+Y1*Y1)/sqrt(X3*X3 + Y3*Y3) );    
      //      double phi_23 = acos( (X2*X3 + Y2*Y3)/sqrt(X2*X2+Y2*Y2)/sqrt(X3*X3 + Y3*Y3) );    

      double dA01 = -99;  // normalized multiple scattering angle VTX-ST1 vector x St1P 

      double eta     = sngmuon->get_rapidity();
      double lastGap = sngmuon->get_lastgap();
      double DG0 = sngmuon->get_DG0();
      double DDG0 = sngmuon->get_DDG0();

      double ntrhits = sngmuon->get_ntrhits();
      double nidhits = sngmuon->get_nidhits();

      double trchi2 = sngmuon->get_trchi2();
      double idchi2 = sngmuon->get_trchi2();

      double dca_z = sngmuon->get_dca_z();
      double dca_r = sngmuon->get_dca_r();

      double dphi_fvtx = sngmuon->get_dphi_fvtx();
      double dtheta_fvtx = sngmuon->get_dtheta_fvtx();
      double dr_fvtx = sngmuon->get_dr_fvtx();

      double Rpc3DCA =  -999;
      double Rpc1DCA =  -999;

      if (run12pp510){
	//      double Rpc3time = sngmuon->get_Rpctime();
	Rpc3DCA = sngmuon->get_RpcDCA();

	//      double Rpc1time = sngmuon->get_Rpc1time();
	Rpc1DCA = sngmuon->get_Rpc1DCA();
      }
      if (run13pp510){
	Rpc3DCA = rpcsngmuon->get_Rpc3St3DCA();
	Rpc1DCA = rpcsngmuon->get_Rpc1St1DCA();
      }


      //
      //--- basic event and track selections ---
      //

      if (fabs(pZ)  < 3.0 or fabs(pZ) > 500) continue;

      if (pT < 1 || pT > 100) continue;

      if (fabs(DG0) > 10.0) continue; 
      if (fabs(DDG0) > 10.0) continue; 

      if (fabs(ntrhits) < 10.0) continue; 
      if (fabs(nidhits) < 6.0) continue; 

      if (fabs(trchi2) > 10.0) continue; 
      if (fabs(idchi2) > 5.0) continue; 

      if (fabs(dca_z) > 100.0) continue;   // not cut 
      if (fabs(dca_r) > 5.0) continue; 

      if (fabs(dphi_fvtx) > 20.0) continue;   // D= -99
      if (fabs(dtheta_fvtx) > 20.0) continue; // D= -99
      if (fabs(dr_fvtx) > 20.0) continue;     // D= -99


      // --- RPC timing cut --- effectively seltected "deep muons" ------ 1H and 1D are different 
      //      if (fabs(Rpc3DCA) > 50.0 && fabs(Rpc1DCA) > 10.0 ) continue;   // D = -9999


      //
      // using the 1st absorber multiple-scattrering to reject fake soft "high pT" tracks
      //


      // determine the primary vtx as the one closest to the Evt_bbcZ (or better muon track candidate)
      // at least one of them are close to the muon track z
      if (fabs(Evt_fvtxZ - Evt_bbcZ) > 100 && fabs(Evt_fvtxZ2 - Evt_bbcZ) > 100) continue ; 


      //select the primary vtx with FVTX, need to use "muon's DCA or Chi2 fit"
      if (fabs(Evt_fvtxZ - Evt_bbcZ) < fabs(Evt_fvtxZ2 - Evt_bbcZ) ){
	Evt_X0 = Evt_fvtxX;
	Evt_Y0 = Evt_fvtxY;
	Evt_Z0 = Evt_fvtxZ;
      }
      if (fabs(Evt_fvtxZ - Evt_bbcZ) > fabs(Evt_fvtxZ2 - Evt_bbcZ) ){
	Evt_X0 = Evt_fvtxX2;
	Evt_Y0 = Evt_fvtxY2;
	Evt_Z0 = Evt_fvtxZ2;
      }

      //calcualte VTX-ST1: dTheta and dPhi; need th eprimary vtx infor

      double v1  = (X1 - Evt_X0)*pX1 + (Y1 - Evt_Y0)*pY1 + (Z1 - Evt_Z0)*pZ1;
      double v2a = sqrt( (X1 - Evt_X0)*(X1 - Evt_X0)  + (Y1 - Evt_Y0)*(Y1 - Evt_Y0) + (Z1 - Evt_Z0)*(Z1 - Evt_Z0) );
      double v2b = sqrt( pX1*pX1 + pY1*pY1 + pZ1*pZ1 );

      dA01 = pZ*acos(v1/v2a/v2b)*sqrt(pTot/abs(pZ));     // singed according to arms 
      if (fabs(dA01)>0.5) continue ;


      //calcualte VTX-ST1-ST3: dTheta and dPhi

      // -- comparision ---
      double dw23 = sin(phi_23)*pT*(pT/pTot) ;


      //
      // this is an event used to fill the histogram
      // 
      flag->flags[imuon] = 1;   

      //
      //--- end of event and track selections ---
      //


      fill_ok = true;   // output this event for later analysis 

      //
      //---- fill general event and track properties ---
      //

      _h_dA->Fill(ep, dA01);

      _h2_dA->Fill(ep, pZ, fabs(dA01) );

      if (pZ > 0) {
	_h_Rpc1DCA_N->Fill(ep,Rpc1DCA);
	_h_Rpc3DCA_N->Fill(ep,Rpc3DCA);
      }
      else {
	_h_Rpc1DCA_S->Fill(ep, Rpc1DCA);
	_h_Rpc3DCA_S->Fill(ep, Rpc3DCA);
      }

      // MuTr phi angles 
      if (pZ > 0) {
	_h_phi_1N->Fill(ep,phi_1);
	_h_phi_2N->Fill(ep,phi_2);
	_h_phi_3N->Fill(ep,phi_3);

	_h_phi_13N->Fill(ep,phi_13);
	_h_phi_23N->Fill(ep,phi_23);

	_h2_phi_13N->Fill(ep,pZ,fabs(pZ)*phi_13);
	_h2_phi_23N->Fill(ep,pZ,fabs(pZ)*phi_23);

	_h_ST1Z_N->Fill(ep,Z1);

      }
      else if (pZ < 0){
	_h_phi_1S->Fill(ep,phi_1);
	_h_phi_2S->Fill(ep,phi_2);
	_h_phi_3S->Fill(ep,phi_3);

	_h_phi_13S->Fill(ep,phi_13);
	_h_phi_23S->Fill(ep,phi_23);

	_h2_phi_13S->Fill(ep,pZ,fabs(pZ)*phi_13);
	_h2_phi_23S->Fill(ep,pZ,fabs(pZ)*phi_23);

	_h_ST1Z_S->Fill(ep,Z1);

      }

      //
      //  --- now for spin asymmetry stuff ---
      //

      //
      // --- deep muon pT dependent asymetry  ----
      //

      if ( (fabs(Rpc3DCA) < 50.0 || fabs(Rpc1DCA) < 10.0) && lastGap > 3 ) {   // D = -9999

	// all sngmuons Q(+,-)      
	_h_pT->Fill(ep, pT);
	
	if (arm ==1){      // north
	  _h_pT_n->Fill(ep, pT);
	}
	if (arm ==0){      // south
	  _h_pT_s->Fill(ep, pT);
	}
	
	// sngmuons Q(+)
	if (charge >0 and arm ==1){      // north
	  _h_pT_p->Fill(ep, pT);
	  _h_pT_p_n->Fill(ep, pT);
	}
	if (charge <0 and arm ==0){      // south
	  _h_pT_m->Fill(ep, pT);
	  _h_pT_m_s->Fill(ep, pT);
	}
	// sngmuons Q(-)
	if (charge >0 and arm ==1){      // north
	  _h_pT_p->Fill(ep, pT);
	  _h_pT_p_n->Fill(ep, pT);
	}
	if (charge <0 and arm ==0){      // south
	  _h_pT_m->Fill(ep, pT);
	  _h_pT_m_s->Fill(ep, pT);
	}
	
      } //  --- en dof deep muon + RPCs


      //
      // --- stopped hadron pT dependent asymetry  ----
      //

      if ( (fabs(Rpc1DCA) < 10.0) && lastGap < 4 ) {   // D = -9999

	// all sngmuons Q(+,-)      
	_h_pTH->Fill(ep, pT);
	
	if (arm ==1){      // north
	  _h_pTH_n->Fill(ep, pT);
	}
	if (arm ==0){      // south
	  _h_pTH_s->Fill(ep, pT);
	}
	
	// sngmuons Q(+)
	if (charge >0 and arm ==1){      // north
	  _h_pTH_p->Fill(ep, pT);
	  _h_pTH_p_n->Fill(ep, pT);
	}
	if (charge <0 and arm ==0){      // south
	  _h_pTH_p->Fill(ep, pT);
	  _h_pTH_p_s->Fill(ep, pT);
	}
	// sngmuons Q(-)
	if (charge <0 and arm ==1){      // north
	  _h_pTH_m->Fill(ep, pT);
	  _h_pTH_m_n->Fill(ep, pT);
	}
	if (charge <0 and arm ==0){      // south
	  _h_pTH_m->Fill(ep, pT);
	  _h_pTH_m_s->Fill(ep, pT);
	}
     
      } // --- end of hadron selection lastGap < 4 & RPC1 


      //
      // -- repidity dependance --- W analysis "deep muon" + RPCs
      //

      // -- hadron background ---
      if ( (fabs(Rpc1DCA) < 10.0) && lastGap < 4 ) {   // D = -9999

	// sngmuons stopped hadrons -> mu +/-
	if (charge >0 and lastGap < 4){     
	  _h_etaH->Fill(ep, eta);
	  _h_etaH_p->Fill(ep, eta);
	}
	if (charge < 0 and lastGap < 4 ){
	  _h_etaH->Fill(ep, eta);
	  _h_etaH_m->Fill(ep, eta);
	}
      }  // -- hadrons -- 


      // -- deep muons ---
      if ( (fabs(Rpc3DCA) < 50.0 || fabs(Rpc1DCA) < 10.0) && lastGap > 3 ) {   // D = -9999
      
	// sngmuons heavy -> mu +/-
	if (charge >0 and lastGap > 3){     
	  _h_eta->Fill(ep, eta);
	  _h_eta_p->Fill(ep, eta);

	  if (pZ > 0) {
	    _h_dw23N_p->Fill(ep, dw23);
	  }
	  else if (pZ < 0) {
	    _h_dw23S_p->Fill(ep, dw23);
	  } //

	}
	if (charge < 0 and lastGap >3 ){
	  _h_eta->Fill(ep, eta);
	  _h_eta_m->Fill(ep, eta);

	  if (pZ > 0) {
	    _h_dw23N_m->Fill(ep, dw23);
	  }
	  else if (pZ < 0) {
	    _h_dw23S_m->Fill(ep, dw23);
	  } //

	}

	// sngmuons W -> mu +/-, pT>10
	if (charge >0 and lastGap >3 and pT > 10){            // W->mu+
	  _h_etaW10->Fill(ep, eta);
	  _h_etaW10_p->Fill(ep, eta);

	  if (pZ > 0) {
	    _h_dw23N_PT10_p->Fill(ep, dw23);
	  }
	  else if (pZ < 0) {
	    _h_dw23S_PT10_p->Fill(ep, dw23);
	  } //

	}

	if  (charge <0 and lastGap >3 and pT > 10){      // W->mu-
	  _h_etaW10->Fill(ep, eta);
	  _h_etaW10_m->Fill(ep, eta);

	  if (pZ > 0) {
	    _h_dw23N_PT10_m->Fill(ep, dw23);
	  }
	  else if (pZ < 0) {
	    _h_dw23S_PT10_m->Fill(ep, dw23);
	  } //

	}

	
	// sngmuons W -> mu +/-, pT>15
	if (charge >0 and lastGap >3 and pT > 15){            // W->mu+
	  _h_etaW15->Fill(ep, eta);
	  _h_etaW15_p->Fill(ep, eta);

	  if (pZ > 0) {
	    _h_dw23N_PT15_p->Fill(ep, dw23);
	  }
	  else if (pZ < 0) {
	    _h_dw23S_PT15_p->Fill(ep, dw23);
	  } //

	}
	if  (charge <0 and lastGap >3 and pT > 15){      // W->mu-
	  _h_etaW15->Fill(ep, eta);
	  _h_etaW15_m->Fill(ep, eta);

	  if (pZ > 0) {
	    _h_dw23N_PT15_m->Fill(ep, dw23);
	  }
	  else if (pZ < 0) {
	    _h_dw23S_PT15_m->Fill(ep, dw23);
	  } //

	}


	// sngmuons W -> mu +/-, pT>20
	if (charge >0 and lastGap >3 and pT > 20){            // W->mu+
	  _h_etaW20->Fill(ep, eta);
	  _h_etaW20_p->Fill(ep, eta);

	  if (pZ > 0) {
	    _h_dw23N_PT20_p->Fill(ep, dw23);
	  }
	  else if (pZ < 0) {
	    _h_dw23S_PT20_p->Fill(ep, dw23);
	  } //

	}
	if  (charge <0 and lastGap >3 and pT > 20){      // W->mu-
	  _h_etaW20->Fill(ep, eta);
	  _h_etaW20_m->Fill(ep, eta);

	  if (pZ > 0) {
	    _h_dw23N_PT20_m->Fill(ep, dw23);
	  }
	  else if (pZ < 0) {
	    _h_dw23S_PT20_m->Fill(ep, dw23);
	  } //

	}

      } // -- end of deep muon + RPCs 
      
      //      fill_ok = true;

    }

  return fill_ok? EVENT_OK : DISCARDEVENT;
}

//! global termination
int
saModuleSngmuonMing::end(PHCompositeNode *topNode, sa_hist_mangager_ptr hm)
{
  //
  //calculate asymmetry with relative lumi of BbcNoCutPHENIX
  //

  // --- pT: Muons --- 
  _h_pT->CalcAsymmetry(hm.get()->getHisto_DefaultLumi(),
      saHist::BbcNoCutPHENIX);
  _h_pT_n->CalcAsymmetry(hm.get()->getHisto_DefaultLumi(),
      saHist::BbcNoCutPHENIX);
  _h_pT_s->CalcAsymmetry(hm.get()->getHisto_DefaultLumi(),
      saHist::BbcNoCutPHENIX);

  _h_pT_p->CalcAsymmetry(hm.get()->getHisto_DefaultLumi(),
      saHist::BbcNoCutPHENIX);
  _h_pT_p_n->CalcAsymmetry(hm.get()->getHisto_DefaultLumi(),
      saHist::BbcNoCutPHENIX);
  _h_pT_p_s->CalcAsymmetry(hm.get()->getHisto_DefaultLumi(),
      saHist::BbcNoCutPHENIX);

  _h_pT_m->CalcAsymmetry(hm.get()->getHisto_DefaultLumi(),
      saHist::BbcNoCutPHENIX);
  _h_pT_m_n->CalcAsymmetry(hm.get()->getHisto_DefaultLumi(),
      saHist::BbcNoCutPHENIX);
  _h_pT_m_s->CalcAsymmetry(hm.get()->getHisto_DefaultLumi(),
      saHist::BbcNoCutPHENIX);

  // --- pT: Hadron --- 

  _h_pTH->CalcAsymmetry(hm.get()->getHisto_DefaultLumi(),
      saHist::BbcNoCutPHENIX);
  _h_pTH_n->CalcAsymmetry(hm.get()->getHisto_DefaultLumi(),
      saHist::BbcNoCutPHENIX);
  _h_pTH_s->CalcAsymmetry(hm.get()->getHisto_DefaultLumi(),
      saHist::BbcNoCutPHENIX);

  _h_pTH_p->CalcAsymmetry(hm.get()->getHisto_DefaultLumi(),
      saHist::BbcNoCutPHENIX);
  _h_pTH_p_n->CalcAsymmetry(hm.get()->getHisto_DefaultLumi(),
      saHist::BbcNoCutPHENIX);
  _h_pTH_p_s->CalcAsymmetry(hm.get()->getHisto_DefaultLumi(),
      saHist::BbcNoCutPHENIX);

  _h_pTH_m->CalcAsymmetry(hm.get()->getHisto_DefaultLumi(),
      saHist::BbcNoCutPHENIX);
  _h_pTH_m_n->CalcAsymmetry(hm.get()->getHisto_DefaultLumi(),
      saHist::BbcNoCutPHENIX);
  _h_pTH_m_s->CalcAsymmetry(hm.get()->getHisto_DefaultLumi(),
      saHist::BbcNoCutPHENIX);


  // -- eta: Muons and Hadrons ---

  _h_eta->CalcAsymmetry(hm.get()->getHisto_DefaultLumi(),
      saHist::BbcNoCutPHENIX);
  _h_eta_p->CalcAsymmetry(hm.get()->getHisto_DefaultLumi(),
      saHist::BbcNoCutPHENIX);
  _h_eta_m->CalcAsymmetry(hm.get()->getHisto_DefaultLumi(),
      saHist::BbcNoCutPHENIX);

  _h_etaH->CalcAsymmetry(hm.get()->getHisto_DefaultLumi(),
      saHist::BbcNoCutPHENIX);
  _h_etaH_p->CalcAsymmetry(hm.get()->getHisto_DefaultLumi(),
      saHist::BbcNoCutPHENIX);
  _h_etaH_m->CalcAsymmetry(hm.get()->getHisto_DefaultLumi(),
      saHist::BbcNoCutPHENIX);

  // -- eta: W ->Muons ---

  _h_etaW10->CalcAsymmetry(hm.get()->getHisto_DefaultLumi(),
      saHist::BbcNoCutPHENIX);
  _h_etaW10_p->CalcAsymmetry(hm.get()->getHisto_DefaultLumi(),
      saHist::BbcNoCutPHENIX);
  _h_etaW10_m->CalcAsymmetry(hm.get()->getHisto_DefaultLumi(),
      saHist::BbcNoCutPHENIX);

  _h_etaW15->CalcAsymmetry(hm.get()->getHisto_DefaultLumi(),
      saHist::BbcNoCutPHENIX);
  _h_etaW15_p->CalcAsymmetry(hm.get()->getHisto_DefaultLumi(),
      saHist::BbcNoCutPHENIX);
  _h_etaW15_m->CalcAsymmetry(hm.get()->getHisto_DefaultLumi(),
      saHist::BbcNoCutPHENIX);

  _h_etaW20->CalcAsymmetry(hm.get()->getHisto_DefaultLumi(),
      saHist::BbcNoCutPHENIX);
  _h_etaW20_p->CalcAsymmetry(hm.get()->getHisto_DefaultLumi(),
      saHist::BbcNoCutPHENIX);
  _h_etaW20_m->CalcAsymmetry(hm.get()->getHisto_DefaultLumi(),
      saHist::BbcNoCutPHENIX);


  return EVENT_OK;
}

