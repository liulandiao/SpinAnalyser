// $Id: saModuleDimuonMing.C,v 1.5 2013/08/17 05:42:47 mxliu Exp $                                                                                             
//
// modified from Jin's saModuleSimpleDimuon  07/12/2013
// for Run12+ dimuon analyses: J/Psi, DY and B2B
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
#include <MCDiMuonContainer.h>
#include <MCDiMuon.h>

#include "saFlag.h"

#include "saModuleDimuonMing.h"


using namespace std;

//-- DiMuon Mass
saModuleDimuonMing::saModuleDimuonMing(const std::string &name) :
  saModuleBase(name),                                     //
  _h_mass(NULL), _h_mass_n(NULL), _h_mass_s(NULL),        //
  _h_mass_p(NULL), _h_mass_p_n(NULL), _h_mass_p_s(NULL),  //
  _h_mass_m(NULL), _h_mass_m_n(NULL), _h_mass_m_s(NULL),  //
  _h_pT(NULL), _h_pT_n(NULL), _h_pT_s(NULL),              //
  _h_pT_p(NULL), _h_pT_p_n(NULL), _h_pT_p_s(NULL),              //
  _h_pT_m(NULL), _h_pT_m_n(NULL), _h_pT_m_s(NULL),              //
  _h_b2b(NULL), _h_b2b_p(NULL), _h_b2b_m(NULL)
{

}



saModuleDimuonMing::~saModuleDimuonMing()
{

}

//! global initialization
int
saModuleDimuonMing::init(PHCompositeNode *topNode, sa_hist_mangager_ptr hm)
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
          "DimuonFlag", "PHObject");

      if (!node)
        {
          cout
              << "saModuleSimpleDimuon::Init failed to create saEventProperty Node"
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
      cout << "saModuleDimuonMing::Init failed to create saEventProperty"
          << endl;
      return ABORTRUN;
    }

  //
  // --- define user histograms ----
  //

  //-- same arm DiMuons Mass (+-)Q
  TH1F * h_mass = new TH1F("DiMuonMass", "Same arm Invariant Mass (+,-) (GeV)", 20, 1, 10);
  _h_mass = new saHist(
      //
      h_mass, //the template histogram, any TH1 and derivatives is accepted
      saHist::DEFAULT_FLAG | saHist::RUN_PROPERTY | saHist::FILL_PROPERTY
          | saHist::CROSSING_CHECKS, // flags
      Verbosity() // verbosity
      );
  hm.get()->registerHisto(_h_mass);

  // --- north ---
  TH1F * h_mass_n = new TH1F("DiMuonMassN", "North Invariant Mass (+,-) (GeV)", 20, 1, 10);
  _h_mass_n = new saHist(
      //
      h_mass_n, //the template histogram, any TH1 and derivatives is accepted
      saHist::DEFAULT_FLAG | saHist::RUN_PROPERTY | saHist::FILL_PROPERTY
          | saHist::CROSSING_CHECKS, // flags
      Verbosity() // verbosity
      );
  hm.get()->registerHisto(_h_mass_n);

  // --- south ---
  TH1F * h_mass_s = new TH1F("DiMuonMassS", "South Invariant Mass (+,-) (GeV)", 20, 1, 10);
  _h_mass_s = new saHist(
      //
      h_mass_s, //the template histogram, any TH1 and derivatives is accepted
      saHist::DEFAULT_FLAG | saHist::RUN_PROPERTY | saHist::FILL_PROPERTY
          | saHist::CROSSING_CHECKS, // flags
      Verbosity() // verbosity
      );
  hm.get()->registerHisto(_h_mass_s);


  //-- same arm DiMuons Mass Q(++)
  TH1F * h_mass_p = new TH1F("DiMuonMassQp", "Q(++)Invariant Mass;Invariant Mass (GeV)", 20, 1, 10);
  _h_mass_p = new saHist(
      //
      h_mass_p, //the template histogram, any TH1 and derivatives is accepted
      saHist::DEFAULT_FLAG | saHist::RUN_PROPERTY | saHist::FILL_PROPERTY
          | saHist::CROSSING_CHECKS, // flags
      Verbosity() // verbosity
      );
  hm.get()->registerHisto(_h_mass_p);


  // --- north muons ---
  TH1F * h_mass_p_n = new TH1F("DiMuonMassQpN", "North Q(++)Invariant Mass;Invariant Mass (GeV)", 20, 1, 10);
  _h_mass_p_n = new saHist(
      //
      h_mass_p_n, //the template histogram, any TH1 and derivatives is accepted
      saHist::DEFAULT_FLAG | saHist::RUN_PROPERTY | saHist::FILL_PROPERTY
          | saHist::CROSSING_CHECKS, // flags
      Verbosity() // verbosity
      );
  hm.get()->registerHisto(_h_mass_p_n);


  // --- south muons ---
  TH1F * h_mass_p_s = new TH1F("DiMuonMassQpS", "South Q(++)Invariant Mass;Invariant Mass (GeV)", 20, 1, 10);
  _h_mass_p_s = new saHist(
      //
      h_mass_p_s, //the template histogram, any TH1 and derivatives is accepted
      saHist::DEFAULT_FLAG | saHist::RUN_PROPERTY | saHist::FILL_PROPERTY
          | saHist::CROSSING_CHECKS, // flags
      Verbosity() // verbosity
      );
  hm.get()->registerHisto(_h_mass_p_s);



  //-- same arm DiMuons Mass Q(--)
  TH1F * h_mass_m = new TH1F("DiMuonMassQm", "Q(--)Invariant Mass;Invariant Mass (GeV)", 20,
      1, 10);
  _h_mass_m = new saHist(
      //
      h_mass_m, //the template histogram, any TH1 and derivatives is accepted
      saHist::DEFAULT_FLAG | saHist::RUN_PROPERTY | saHist::FILL_PROPERTY
          | saHist::CROSSING_CHECKS, // flags
      Verbosity() // verbosity
      );
  hm.get()->registerHisto(_h_mass_m);


  // --- north ---
  TH1F * h_mass_m_n = new TH1F("DiMuonMassQmN", "North Q(--)Invariant Mass;Invariant Mass (GeV)", 20,
      1, 10);
  _h_mass_m_n = new saHist(
      //
      h_mass_m_n, //the template histogram, any TH1 and derivatives is accepted
      saHist::DEFAULT_FLAG | saHist::RUN_PROPERTY | saHist::FILL_PROPERTY
          | saHist::CROSSING_CHECKS, // flags
      Verbosity() // verbosity
      );
  hm.get()->registerHisto(_h_mass_m_n);


  // --- south ---
  TH1F * h_mass_m_s = new TH1F("DiMuonMassQmS", "South Q(--)Invariant Mass;Invariant Mass (GeV)", 20,
      1, 10);
  _h_mass_m_s = new saHist(
      //
      h_mass_m_s, //the template histogram, any TH1 and derivatives is accepted
      saHist::DEFAULT_FLAG | saHist::RUN_PROPERTY | saHist::FILL_PROPERTY
          | saHist::CROSSING_CHECKS, // flags
      Verbosity() // verbosity
      );
  hm.get()->registerHisto(_h_mass_m_s);


  //-- same arm DiMuons pT
  TH1F * h_pT = new TH1F("DiMuonPT", "DiMuon pT; DiMuon pT (GeV)", 20,
      0, 10);
  _h_pT = new saHist(
      //
      h_pT, //the template histogram, any TH1 and derivatives is accepted
      saHist::DEFAULT_FLAG | saHist::RUN_PROPERTY | saHist::FILL_PROPERTY
          | saHist::CROSSING_CHECKS, // flags
      Verbosity() // verbosity
      );
  hm.get()->registerHisto(_h_pT);

  // ---- north ---
  TH1F * h_pT_n = new TH1F("DiMuonPTN", "North DiMuon pT; DiMuon pT (GeV)", 20,
      0, 10);
  _h_pT_n = new saHist(
      //
      h_pT_n, //the template histogram, any TH1 and derivatives is accepted
      saHist::DEFAULT_FLAG | saHist::RUN_PROPERTY | saHist::FILL_PROPERTY
          | saHist::CROSSING_CHECKS, // flags
      Verbosity() // verbosity
      );
  hm.get()->registerHisto(_h_pT_n);

  // ---- south ---
  TH1F * h_pT_s = new TH1F("DiMuonPTS", "South DiMuon pT; DiMuon pT (GeV)", 20,
      0, 10);
  _h_pT_s = new saHist(
      //
      h_pT_s, //the template histogram, any TH1 and derivatives is accepted
      saHist::DEFAULT_FLAG | saHist::RUN_PROPERTY | saHist::FILL_PROPERTY
          | saHist::CROSSING_CHECKS, // flags
      Verbosity() // verbosity
      );
  hm.get()->registerHisto(_h_pT_s);

  // -- Q(++)
  TH1F * h_pT_p = new TH1F("DiMuonPTp", "DiMuon(++) pT; DiMuon pT (GeV)", 20,
      0, 10);
  _h_pT_p = new saHist(
      //
      h_pT_p, //the template histogram, any TH1 and derivatives is accepted
      saHist::DEFAULT_FLAG | saHist::RUN_PROPERTY | saHist::FILL_PROPERTY
          | saHist::CROSSING_CHECKS, // flags
      Verbosity() // verbosity
      );
  hm.get()->registerHisto(_h_pT_p);

  // ---- north ---
  TH1F * h_pT_p_n = new TH1F("DiMuonPTpN", "North DiMuon (++) pT; DiMuon pT (GeV)", 20,
      0, 10);
  _h_pT_p_n = new saHist(
      //
      h_pT_p_n, //the template histogram, any TH1 and derivatives is accepted
      saHist::DEFAULT_FLAG | saHist::RUN_PROPERTY | saHist::FILL_PROPERTY
          | saHist::CROSSING_CHECKS, // flags
      Verbosity() // verbosity
      );
  hm.get()->registerHisto(_h_pT_p_n);

  // ---- south ---
  TH1F * h_pT_p_s = new TH1F("DiMuonPTpS", "South DiMuon (++) pT; DiMuon pT (GeV)", 20,
      0, 10);
  _h_pT_p_s = new saHist(
      //
      h_pT_p_s, //the template histogram, any TH1 and derivatives is accepted
      saHist::DEFAULT_FLAG | saHist::RUN_PROPERTY | saHist::FILL_PROPERTY
          | saHist::CROSSING_CHECKS, // flags
      Verbosity() // verbosity
      );
  hm.get()->registerHisto(_h_pT_p_s);

  // -- Q(--)
  TH1F * h_pT_m = new TH1F("DiMuonPTm", "DiMuon(--) pT; DiMuon pT (GeV)", 20,
      0, 10);
  _h_pT_m = new saHist(
      //
      h_pT_m, //the template histogram, any TH1 and derivatives is accepted
      saHist::DEFAULT_FLAG | saHist::RUN_PROPERTY | saHist::FILL_PROPERTY
          | saHist::CROSSING_CHECKS, // flags
      Verbosity() // verbosity
      );
  hm.get()->registerHisto(_h_pT_m);

  // ---- north ---
  TH1F * h_pT_m_n = new TH1F("DiMuonPTmN", "North DiMuon (--) pT; DiMuon pT (GeV)", 20,
      0, 10);
  _h_pT_m_n = new saHist(
      //
      h_pT_m_n, //the template histogram, any TH1 and derivatives is accepted
      saHist::DEFAULT_FLAG | saHist::RUN_PROPERTY | saHist::FILL_PROPERTY
          | saHist::CROSSING_CHECKS, // flags
      Verbosity() // verbosity
      );
  hm.get()->registerHisto(_h_pT_m_n);

  // ---- south ---
  TH1F * h_pT_m_s = new TH1F("DiMuonPTmS", "South DiMuon(--) pT; DiMuon pT (GeV)", 20,
      0, 10);
  _h_pT_m_s = new saHist(
      //
      h_pT_m_s, //the template histogram, any TH1 and derivatives is accepted
      saHist::DEFAULT_FLAG | saHist::RUN_PROPERTY | saHist::FILL_PROPERTY
          | saHist::CROSSING_CHECKS, // flags
      Verbosity() // verbosity
      );
  hm.get()->registerHisto(_h_pT_m_s);



  //-- B2B DiMuons
  TH1F * h_b2b = new TH1F("DiMuonB2B", "B2B DiMuon Mass; DiMuon Mass (GeV)", 20,
      5, 15);
  _h_b2b = new saHist(
      //
      h_b2b, //the template histogram, any TH1 and derivatives is accepted
      saHist::DEFAULT_FLAG | saHist::RUN_PROPERTY | saHist::FILL_PROPERTY
          | saHist::CROSSING_CHECKS, // flags
      Verbosity() // verbosity
      );
  hm.get()->registerHisto(_h_b2b);


  //-- B2B DiMuons Q(++)
  TH1F * h_b2b_p = new TH1F("DiMuonB2BQp", "Q(++) B2B DiMuon Mass; DiMuon Mass (GeV)", 20,
      5, 15);
  _h_b2b_p = new saHist(
      //
      h_b2b_p, //the template histogram, any TH1 and derivatives is accepted
      saHist::DEFAULT_FLAG | saHist::RUN_PROPERTY | saHist::FILL_PROPERTY
          | saHist::CROSSING_CHECKS, // flags
      Verbosity() // verbosity
      );
  hm.get()->registerHisto(_h_b2b_p);


  //-- B2B DiMuons Q(--)
  TH1F * h_b2b_m = new TH1F("DiMuonB2BQm", "Q(--) B2B DiMuon Mass; DiMuon Mass (GeV)", 20,
      5, 15);
  _h_b2b_m = new saHist(
      //
      h_b2b_m, //the template histogram, any TH1 and derivatives is accepted
      saHist::DEFAULT_FLAG | saHist::RUN_PROPERTY | saHist::FILL_PROPERTY
          | saHist::CROSSING_CHECKS, // flags
      Verbosity() // verbosity
      );
  hm.get()->registerHisto(_h_b2b_m);



  return EVENT_OK;
}

//! Run initialization
int
saModuleDimuonMing::init_run(PHCompositeNode *topNode,
    sa_hist_mangager_ptr hm)
{

  return EVENT_OK;
}

//! event method
int
saModuleDimuonMing::event(PHCompositeNode *topNode, sa_hist_mangager_ptr hm)
{

  DiMuonContainer* dimuons = findNode::getClass<DiMuonContainer>(topNode,
      "DiMuonContainer");
  if (!dimuons)
    {
      cout << "saModuleDimuonMing:: DiMuonContainer - ERROR - not in Node Tree"
          << endl;
      return ABORTRUN;
    }

  const saEventProperty * ep = findNode::getClass<saEventProperty>(topNode,
      "saEventProperty");
  if (!ep)
    {
      cout << "saModuleDimuonMing::event - ERROR  - Cannot find EventProperty node in Top Node"
          << endl;
      return ABORTRUN;
    }

  if (Verbosity() >= 2)
    {
      cout << "saModuleDimuonMing::event - INFO  - ";
      ep->identify();
    }

  saFlagC *flag = findNode::getClass<saFlagC>(topNode, "DimuonFlag");
  if (!flag)
    {
      cout << "saModuleDimuonMing::event - DimuonFlag not in Node Tree"
          << endl;
      return ABORTRUN;
    }

  // syncronizely set the cut flag with dimuon container
  flag->flags.Set(dimuons->get_nDiMuons());
  flag->flags.Reset(0);

  // save the event for DST files?
  //  int return_code = dimuons->get_nDiMuons() >0 ? EVENT_OK : DISCARDEVENT;

  bool fill_ok = false;   // save the event for DST output 

  // determine this event's primary VTX
      
  double Evt_fvtxZ = dimuons->get_Evt_fvtxZ();
  double Evt_fvtxZ2 = dimuons->get_Evt_fvtxZ2();
  double Evt_bbcZ = dimuons->get_Evt_fvtxZ();
    
  //  int fill_number = ep->get_fill_number();
  int run_number = ep->get_run_number();

  // good run12pp510 run # = 366283 - 366059
  // good run13pp510 run # = 386773 - 397990

  bool run12pp510 = false;
  bool run13pp510 = false;

  if (run_number > 360000 && run_number < 370000) run12pp510 = true;
  if (run_number > 380000 && run_number < 400000) run13pp510 = true;
 

  if ( fabs(Evt_bbcZ) > 100) {

    return DISCARDEVENT;
  }

  if ( fabs(Evt_fvtxZ) > 100) {

    return DISCARDEVENT;
  }


  for (unsigned int imuon = 0; imuon < dimuons->get_nDiMuons(); imuon++)
    {
      DiMuon* dimuon = dimuons->get_DiMuon(imuon);


      int type = 99;
      if (dimuon->get_charge() == 0 && dimuon->get_same_event())
        type = 0;
      if (dimuon->get_charge() < 0 && dimuon->get_same_event())
        type = 1;
      if (dimuon->get_charge() > 0 && dimuon->get_same_event())
        type = 2;

      if (dimuon->get_charge() == 0 && !dimuon->get_same_event())
        type = 3;
      if (dimuon->get_charge() < 0 && !dimuon->get_same_event())
        type = 4;
      if (dimuon->get_charge() > 0 && !dimuon->get_same_event())
        type = 5;

      int arm0 = (dimuon->get_Tr0_pz() > 0) ? 1 : 0;
      int arm1 = (dimuon->get_Tr1_pz() > 0) ? 1 : 0;

      double pX_DiMuon = dimuon->get_Tr0_px()+dimuon->get_Tr1_px();
      double pY_DiMuon = dimuon->get_Tr0_py()+dimuon->get_Tr1_py();
      double pZ_DiMuon = dimuon->get_Tr0_pz()+dimuon->get_Tr1_pz();

      double pT_DiMuon = sqrt(pX_DiMuon*pX_DiMuon +  pY_DiMuon*pY_DiMuon); 
      //      int charge = dimuon->get_charge(); // Q= -1, 0, +1

      double Evt_vtxchi2 = dimuon->get_Evt_vtxchi2();
      double Evt_vtxoor = dimuon->get_Evt_vtxoor();

      double X0 = dimuon->get_X0();
      double Y0 = dimuon->get_Y0();
      double Z0 = dimuon->get_Z0();

      double dca_r = dimuon->get_dca_r();
      double dca_z = dimuon->get_dca_z();

      double Tr0_DG0 = dimuon->get_Tr0_DG0();
      double Tr1_DG0 = dimuon->get_Tr1_DG0();

      double Tr0_DDG0 = dimuon->get_Tr0_DDG0();
      double Tr1_DDG0 = dimuon->get_Tr1_DDG0();

      double Tr0_ntrhits = dimuon->get_Tr0_ntrhits();
      double Tr1_ntrhits = dimuon->get_Tr1_ntrhits();

      double Tr0_lastgap = dimuon->get_Tr0_lastgap();
      double Tr1_lastgap = dimuon->get_Tr1_lastgap();

      double Tr0_dca_r = dimuon->get_Tr0_dca_r();;
      double Tr1_dca_r = dimuon->get_Tr1_dca_r();;


      double Tr0_Rpc3DCA = -999;
      double Tr1_Rpc3DCA = -999;

      double Tr0_Rpc1DCA = -999;
      double Tr1_Rpc1DCA = -999;

      //RPC information, run12 and run13 are different 

      if (run12pp510){
	Tr0_Rpc3DCA = dimuon->get_Tr0_RpcDCA();
        Tr1_Rpc3DCA = dimuon->get_Tr1_RpcDCA();

	Tr0_Rpc1DCA = dimuon->get_Tr0_Rpc1DCA();
	Tr1_Rpc1DCA = dimuon->get_Tr1_Rpc1DCA();
      }
      if (run13pp510){
	Tr0_Rpc3DCA = dimuon->get_Tr0_Rpc3St3DCA();
	Tr1_Rpc3DCA = dimuon->get_Tr1_Rpc3St3DCA();

	Tr0_Rpc1DCA = dimuon->get_Tr0_Rpc1St1DCA();
	Tr1_Rpc1DCA = dimuon->get_Tr1_Rpc1St1DCA();
      }

      //FVTX information 
      int Tr0_nfvtx_clusters_cone = dimuon->get_Tr0_nfvtx_clusters_cone(7); // 0-7
      int Tr0_nfvtx_tracklets_cone = dimuon->get_Tr0_nfvtx_tracklets_cone(7);

      int Tr1_nfvtx_clusters_cone = dimuon->get_Tr0_nfvtx_clusters_cone(7);
      int Tr1_nfvtx_tracklets_cone = dimuon->get_Tr0_nfvtx_tracklets_cone(7);

      //
      // --- event and track selections ----
      //
      // --- same event dimuon selection 
      if (type > 2)  continue;   // dimuon tracks from the same event only: type = 0, 1, 2

      if (fabs(pZ_DiMuon)>100) continue ;
      if (fabs(pT_DiMuon)>10) continue ;

      if (fabs(dca_r)>10) continue ;
      if (fabs(dca_z)>100) continue ;   // not used 
      if (Evt_vtxchi2 > 5) continue ;
      if (Evt_vtxoor > 1) continue ;

      if (sqrt(X0*X0 + Y0*Y0) > 2.0 ) continue;  // must be close to the bean axis
      
      if (fabs(Evt_bbcZ) > 100) continue ;
      if (fabs(Evt_fvtxZ - Z0) > 10 && fabs(Evt_fvtxZ2 - Z0) > 10) continue ; // at least one of them are close 
      
      if (sqrt(Tr0_DG0*Tr0_DG0 + Tr1_DG0*Tr1_DG0) > 15) continue;
      if (sqrt(Tr0_DDG0*Tr0_DDG0 + Tr1_DDG0*Tr1_DDG0) > 10) continue;

      if (Tr0_ntrhits < 10) continue;
      if (Tr1_ntrhits < 10) continue;

      if (Tr0_lastgap < 3) continue;
      if (Tr1_lastgap < 3) continue;

      if (Tr0_dca_r > 5) continue;
      if (Tr1_dca_r > 5) continue;

      //at least one RPC1 or RPC2 hits for both track at least
      if ( (Tr0_Rpc3DCA <-10 || Tr0_Rpc3DCA > 50 ) and (Tr0_Rpc1DCA <-10 || Tr0_Rpc1DCA > 10 ) ) continue;
      if ( (Tr1_Rpc3DCA <-10 || Tr1_Rpc3DCA > 50 ) and (Tr1_Rpc1DCA <-10 || Tr1_Rpc1DCA > 10 ) ) continue;


      //FVTX cone cuts - just require with FVTX cone-cut variables
      if (Tr0_nfvtx_tracklets_cone < -1 && Tr1_nfvtx_tracklets_cone < -1 ) continue;
      if (Tr0_nfvtx_clusters_cone < -1 && Tr1_nfvtx_clusters_cone < -1 ) continue;


      // --- end of event/track selection ---


      // this is an event used to fill the histogram 
      flag->flags[imuon] = 1;   


      // same arm dimuons Q(+-)
      if (type == 0 and arm0 == arm1)            // same arm (+,-)
	{
        _h_mass->Fill(ep, dimuon->get_mass());
	_h_pT->Fill(ep, pT_DiMuon);
	}
      if (type ==0 and arm0==1 and arm1==1){      // north
        _h_mass_n->Fill(ep, dimuon->get_mass());
	_h_pT_n->Fill(ep, pT_DiMuon);
      }
      if (type ==0 and arm0==0 and arm1==0){    // south
        _h_mass_s->Fill(ep, dimuon->get_mass());
	_h_pT_s->Fill(ep, pT_DiMuon);
      }

      // same arm North dimuons Q(++)
      if (type == 2 and arm0 == arm1)
	{
        _h_mass_p->Fill(ep, dimuon->get_mass());
	_h_pT_p->Fill(ep, pT_DiMuon);
	}
      if (type == 2 and arm0==1 and arm1==1){    // north
	_h_mass_p_n->Fill(ep, dimuon->get_mass());
	_h_pT_p_n->Fill(ep, pT_DiMuon);
      }
      if (type == 2 and arm0==0 and arm1==0){    // south
	_h_mass_p_s->Fill(ep, dimuon->get_mass());
   	_h_pT_p_s->Fill(ep, pT_DiMuon);
      }

      // same arm dimuons Q(--)
      if (type == 1 and arm0 == arm1)
	{
        _h_mass_m->Fill(ep, dimuon->get_mass());
	_h_pT_m->Fill(ep, pT_DiMuon);
	}
      if (type ==1 and arm0==1 and arm1==1){    // north
	_h_mass_m_n->Fill(ep, dimuon->get_mass());
	_h_pT_m_n->Fill(ep, pT_DiMuon);
      }
      if (type ==1 and arm0==0 and arm1==0){    // south
	_h_mass_m_s->Fill(ep, dimuon->get_mass());
	_h_pT_m_s->Fill(ep, pT_DiMuon);
      }

      //
      //opposite arm B2B dimuons Q(+-)
      //
      if (type == 0 and arm0 != arm1)
	{
        _h_b2b->Fill(ep, dimuon->get_mass());
	}
      //opposite arm dimuons Q(++)
      if (type == 2 and arm0 != arm1)
	{
        _h_b2b_p->Fill(ep, dimuon->get_mass());
	}
      //opposite arm dimuons Q(--)
      if (type == 1 and arm0 != arm1)
	{
        _h_b2b_m->Fill(ep, dimuon->get_mass());
	}

      fill_ok = true;
    }

  return fill_ok? EVENT_OK : DISCARDEVENT;
}

//! global termination
int
saModuleDimuonMing::end(PHCompositeNode *topNode, sa_hist_mangager_ptr hm)
{

  //calculate asymmetry with relative lumi of BbcNoCutPHENIX
  _h_mass->CalcAsymmetry(hm.get()->getHisto_DefaultLumi(),
      saHist::BbcNoCutPHENIX);
  _h_mass_n->CalcAsymmetry(hm.get()->getHisto_DefaultLumi(),
      saHist::BbcNoCutPHENIX);
  _h_mass_s->CalcAsymmetry(hm.get()->getHisto_DefaultLumi(),
      saHist::BbcNoCutPHENIX);


  _h_mass_p->CalcAsymmetry(hm.get()->getHisto_DefaultLumi(),
      saHist::BbcNoCutPHENIX);
  _h_mass_p_n->CalcAsymmetry(hm.get()->getHisto_DefaultLumi(),
      saHist::BbcNoCutPHENIX);
  _h_mass_p_s->CalcAsymmetry(hm.get()->getHisto_DefaultLumi(),
      saHist::BbcNoCutPHENIX);


  _h_mass_m->CalcAsymmetry(hm.get()->getHisto_DefaultLumi(),
      saHist::BbcNoCutPHENIX);
  _h_mass_m_n->CalcAsymmetry(hm.get()->getHisto_DefaultLumi(),
      saHist::BbcNoCutPHENIX);
  _h_mass_m_s->CalcAsymmetry(hm.get()->getHisto_DefaultLumi(),
      saHist::BbcNoCutPHENIX);

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


  // B2B dimuons
  _h_b2b->CalcAsymmetry(hm.get()->getHisto_DefaultLumi(),
      saHist::BbcNoCutPHENIX);

  _h_b2b_p->CalcAsymmetry(hm.get()->getHisto_DefaultLumi(),
      saHist::BbcNoCutPHENIX);

  _h_b2b_m->CalcAsymmetry(hm.get()->getHisto_DefaultLumi(),
      saHist::BbcNoCutPHENIX);





  return EVENT_OK;
}

