// $Id: saModuleDimuonDYDarshana.C created by jinhuang modified for Drell-Yan Analysis By Darshana Perera 
#include <Fun4AllReturnCodes.h>
#include <PHCompositeNode.h>
#include <getClass.h>
#include <MWGConsts.h>
#include <Tools.h>
#include <Fun4AllServer.h>
#include <Fun4AllHistoManager.h>
#include <TH1F.h>
#include <TriggerHelper.h>
#include <TrigRunLvl1.h>
#include <TrigLvl1.h>
#include <TH3.h>
#include <TH2.h>
#include <TF1.h>
#include <SyncObject.h>
#include <RunHeader.h>
#include <TMath.h>
#include <SingleMuonContainer.h>
#include <DiMuonContainer.h>
#include <DiMuon.h>
#include <MCDiMuonContainer.h>
#include <MCDiMuon.h>
#include <TFvtxCompactTrk.h>
#include <TFvtxCompactTrkMap.h>
#include "saFlag.h"
#include <SingleMuon.h>
#include <Fun4AllReturnCodes.h>
#include <PHNodeIterator.h>
#include <getClass.h>
#include <TClonesArray.h>
#include "saModuleDimuonDYDarshana.h"

using namespace std;

saModuleDimuonDYDarshana::saModuleDimuonDYDarshana(const std::string &name) :
	saModuleBase(name),
        //_h_mass(NULL),
	_h_pT_os_N(NULL),  _h_pT_os_S(NULL),
	_h_pT_sm_N(NULL),  _h_pT_sm_S(NULL)
	//_h_pT_al_os(NULL), _h_pT_al_ss(NULL),
	//_h_MS_os_N(NULL),  _h_MS_os_S(NULL),
	//_h_MS_sm_N(NULL),  _h_MS_sm_S(NULL)
	//_h_MS_al_os(NULL), _h_MS_al_ss(NULL)

{
 	_use_cut_tracklet_chi2=false;
	_cut_tracklet_chi2=4.;
	_cut_tracklet_dcar=1.5;
	_use_2_hit_tracklet=false;

	_pT_bins[0] = 0.0;
	_pT_bins[1] = 2.0;
	_pT_bins[2] = 10.0;
	//_pT_bins[3] = 10.;

	//_MS_bins[0] = 4.0;
	//_MS_bins[1] = 5.0;
	//_MS_bins[2] = 8.0;
	//_MS_bins[3] = 8.;

}

saModuleDimuonDYDarshana::~saModuleDimuonDYDarshana()
{

}

//! global initialization
int saModuleDimuonDYDarshana::init(PHCompositeNode *topNode, sa_hist_mangager_ptr hm)
{
//-------------------------------------------------------------------
	TH1F * h_pT_os_N = new TH1F("pT_os_sig_n","opposite sign, signal region, north", 2, _pT_bins);
	_h_pT_os_N = new saHist(h_pT_os_N,saHist::DEFAULT_FLAG | saHist::RUN_PROPERTY | 
        saHist::FILL_PROPERTY| saHist::CROSSING_CHECKS,Verbosity());
	hm.get()->registerHisto(_h_pT_os_N);

	TH1F * h_pT_os_S = new TH1F("pT_os_sig_s","opposite sign, signal region, south", 2, _pT_bins);
	_h_pT_os_S = new saHist(h_pT_os_S,saHist::DEFAULT_FLAG | saHist::RUN_PROPERTY | 
        saHist::FILL_PROPERTY| saHist::CROSSING_CHECKS,Verbosity());
	hm.get()->registerHisto(_h_pT_os_S);

	TH1F * h_pT_sm_N = new TH1F("pT_ss_n","same sign, signal region, north", 2, _pT_bins);
	_h_pT_sm_N = new saHist(h_pT_sm_N,saHist::DEFAULT_FLAG | saHist::RUN_PROPERTY | 
        saHist::FILL_PROPERTY| saHist::CROSSING_CHECKS,Verbosity());
	hm.get()->registerHisto(_h_pT_sm_N);

	TH1F * h_pT_sm_S = new TH1F("pT_ss_s","same sign, signal region, south", 2, _pT_bins);
	_h_pT_sm_S = new saHist(h_pT_sm_S,saHist::DEFAULT_FLAG | saHist::RUN_PROPERTY | 
        saHist::FILL_PROPERTY| saHist::CROSSING_CHECKS,Verbosity());
	hm.get()->registerHisto(_h_pT_sm_S);

/*	TH1F * h_pT_al_os = new TH1F("pT_os_all","opposite sign, all ", 3, _pT_bins);
	_h_pT_al_os = new saHist(h_pT_al_os,saHist::DEFAULT_FLAG | saHist::RUN_PROPERTY | 
        saHist::FILL_PROPERTY| saHist::CROSSING_CHECKS,Verbosity());
	hm.get()->registerHisto(_h_pT_al_os);

	TH1F * h_pT_al_ss = new TH1F("pT_ss_all","same sign, all ", 3, _pT_bins);
	_h_pT_al_ss = new saHist(h_pT_al_ss,saHist::DEFAULT_FLAG | saHist::RUN_PROPERTY | 
        saHist::FILL_PROPERTY| saHist::CROSSING_CHECKS,Verbosity());
	hm.get()->registerHisto(_h_pT_al_ss);

*//*
	TH1F * h_MS_os_N = new TH1F("Ms_os_sig_n","opposite sign, north arm",2,_MS_bins);
	_h_MS_os_N = new saHist(h_MS_os_N,saHist::DEFAULT_FLAG | saHist::RUN_PROPERTY | 
        saHist::FILL_PROPERTY| saHist::CROSSING_CHECKS,Verbosity());
	hm.get()->registerHisto(_h_MS_os_N);

	TH1F * h_MS_os_S = new TH1F("Ms_os_sig_s","opposite sign, south arm",2,_MS_bins);
	_h_MS_os_S = new saHist(h_MS_os_S,saHist::DEFAULT_FLAG | saHist::RUN_PROPERTY | 
        saHist::FILL_PROPERTY| saHist::CROSSING_CHECKS,Verbosity());
	hm.get()->registerHisto(_h_MS_os_S);

	TH1F * h_MS_sm_N = new TH1F("Ms_ss_n","same sign, north arm",2,_MS_bins);
	_h_MS_sm_N = new saHist(h_MS_sm_N,saHist::DEFAULT_FLAG | saHist::RUN_PROPERTY | 
        saHist::FILL_PROPERTY| saHist::CROSSING_CHECKS,Verbosity());
	hm.get()->registerHisto(_h_MS_sm_N);

	TH1F * h_MS_sm_S = new TH1F("Ms_ss_s","same sign, south arm",2,_MS_bins);
	_h_MS_sm_S = new saHist(h_MS_sm_S,saHist::DEFAULT_FLAG | saHist::RUN_PROPERTY | 
        saHist::FILL_PROPERTY| saHist::CROSSING_CHECKS,Verbosity());
	hm.get()->registerHisto(_h_MS_sm_S);*/
/*
	TH1F * h_MS_al_os = new TH1F("Ms_os_all","opposite sign, all",2,4, 8);
	_h_MS_al_os = new saHist(h_MS_al_os,saHist::DEFAULT_FLAG | saHist::RUN_PROPERTY | 
        saHist::FILL_PROPERTY| saHist::CROSSING_CHECKS,Verbosity());
	hm.get()->registerHisto(_h_MS_al_os);

	TH1F * h_MS_al_ss = new TH1F("Ms_ss_all","same sign, all", 2,4, 8);
	_h_MS_al_ss = new saHist(h_MS_al_ss,saHist::DEFAULT_FLAG | saHist::RUN_PROPERTY | 
        saHist::FILL_PROPERTY| saHist::CROSSING_CHECKS,Verbosity());
	hm.get()->registerHisto(_h_MS_al_ss);
*/

//-------------------------------------------------------------------


 /* // make a simple histogram //the template histogram, any TH1 and derivatives is accepted
  TH1F * h_mass = new TH1F("InvMass", "Invariant Mass;Invariant Mass (GeV)", 2,4, 8);
  _h_mass = new saHist(h_mass, saHist::DEFAULT_FLAG | saHist::RUN_PROPERTY | 
  saHist::FILL_PROPERTY| saHist::CROSSING_CHECKS,Verbosity());
  hm.get()->registerHisto(_h_mass);
*/
  // make the flag node which is syncronized with cut on picodst_object
  PHCompositeNode* dstNode = NULL;
  PHNodeIterator nodeItr(topNode);
  dstNode = static_cast<PHCompositeNode*>(nodeItr.findFirst("PHCompositeNode","DST"));
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
          "SimpleDimuonFlag", "PHObject");

      if (!node)
        {

          cout
              << "saModuleDimuonDYDarshana::Init failed to create saEventProperty Node"
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
      cout << "saModuleDimuonDYDarshana::Init failed to create saEventProperty"
          << endl;
      return ABORTRUN;
    }

  return EVENT_OK;
}

//! Run initialization
int saModuleDimuonDYDarshana::init_run(PHCompositeNode *topNode,sa_hist_mangager_ptr hm)
{
  return EVENT_OK;
}

//! event method
int saModuleDimuonDYDarshana::event(PHCompositeNode *topNode, sa_hist_mangager_ptr hm)
{
//---------------------------------------------------------------------------------------------------------
   dimuons = findNode::getClass<DiMuonContainer>(topNode,"DiMuonContainer");
  if (!dimuons)
    {
      cout<< "saModuleDimuonDYDarshana:: DiMuonContainer - ERROR - not in Node Tree"<< endl;
      return ABORTRUN;
    }
//---------------------------------------------------------------------------------------------------------
  const saEventProperty * ep = findNode::getClass<saEventProperty>(topNode,"saEventProperty");
  if (!ep)
    {
      cout<< "saModuleDimuonDYDarshana::event - ERROR  - Cannot find EventProperty node in Top Node"<< endl;
      return ABORTRUN;
    }
//---------------------------------------------------------------------------------------------------------

	singlemuoncontainer = findNode::getClass<SingleMuonContainer>(topNode,"SingleMuonContainer");
	if (!singlemuoncontainer){
		cout << "dnpfilter:: SingleMuonContainer not in Node Tree" << endl;
		return ABORTRUN;
	}

//-------------------------------------------------------------------------------------------
	fvtx_trk_map = findNode::getClass<DSTReader>(topNode,"DSTReader");
	if(!fvtx_trk_map){
		cout << PHWHERE << "dnpfilter:: DSTReader not in Node Tree" << endl;
		return ABORTRUN;
	}
//------------------------------------------------------------------------------------------------
  if (Verbosity() >= 2)
    {
      cout << "saModuleDimuonDYDarshana::event - INFO  - ";
      ep->identify();
    }

//---------------------------------------------------------------------------------------------------------
  saFlagC *flag = findNode::getClass<saFlagC>(topNode, "SimpleDimuonFlag");
  if (!flag)
    {
      cout << "saModuleDimuonDYDarshana::event - SimpleDimuonFlag not in Node Tree"<< endl;
      return ABORTRUN;
    }

//---------------------------------------------------------------------------------------------------------
	TrigLvl1* triglv1 = findNode::getClass<TrigLvl1>(topNode,"TrigLvl1");
	if (!triglv1)
	{
		cout << "saModuleDimuonDYDarshana:: TrigLvl1 - ERROR - not in Node Tree"<< endl;
		return ABORTRUN;
	}
//---------------------------------------------------------------------------------------------------------
  // syncronizely set the cut flag with dimuon container
  flag->flags.Set(dimuons->get_nDiMuons());
  flag->flags.Reset(0);
//---------------------------------------------------------------------------------------------------------

  bool save_event = false;

  double Evt_bbcZ = dimuons->get_Evt_bbcZ();

  if ( fabs(Evt_bbcZ) >= 10 ) return DISCARDEVENT;
    
  for (unsigned int imuon = 0; imuon < dimuons->get_nDiMuons(); imuon++)
    {

      //DiMuon *dimuon = dynamic_cast<DiMuon*> (dimu_container->get_DiMuon(idimu));

      dimuon = dynamic_cast<DiMuon*> (dimuons->get_DiMuon(imuon));

      enum EVENT_TYPE type = BAD_EVENT;
      enum ARM_TYPE arm;
      arm = (dimuon->get_rapidity() > 0) ? NORTH : SOUTH;


      nTracklets = GetNumberofTracklets(fvtx_trk_map,dimuon);

      if ( !passesCuts() ) continue;

     // if (nTracklets >= 5) continue;
      if (nTracklets <= 4 || nTracklets >=10 ) continue;
     
      cout<< "*****************Number of tracklets = "<<nTracklets<<" **************Mass = "<<dimuon->get_mass()<<" ******evt  bbc = "<<dimuons->get_Evt_bbcZ()<<endl;
      

if (dimuon->get_charge() == 0 && arm==NORTH) {type = OS_SIG_N;}

if (dimuon->get_charge() == 0 && arm==SOUTH) {type = OS_SIG_S;}

if (dimuon->get_charge() != 0 && arm==NORTH) {type = SS_SIG_N;}

if (dimuon->get_charge() != 0 && arm==SOUTH) {type = SS_SIG_S;}

//if (dimuon->get_charge() == 0) {type = OS_ALL_;}

//if (dimuon->get_charge() == 0) {type = SS_ALL_;}

save_event = true;

flag->flags[imuon] = 1;

if(type==OS_SIG_N) {
_h_pT_os_N->Fill(ep,dimuon->get_pT());}
//_h_MS_os_N->Fill(ep,dimuon->get_mass());}

if(type==OS_SIG_S) {
_h_pT_os_S->Fill(ep,dimuon->get_pT());}
//_h_MS_os_S->Fill(ep,dimuon->get_mass());}

if(type==SS_SIG_N) {
_h_pT_sm_N->Fill(ep,dimuon->get_pT());}
//_h_MS_sm_N->Fill(ep,dimuon->get_mass());}

if(type==SS_SIG_S) {
_h_pT_sm_S->Fill(ep,dimuon->get_pT());}
//_h_MS_sm_S->Fill(ep,dimuon->get_mass());}

//if(type==OS_ALL_)  {
//_h_pT_al_os->Fill(ep,dimuon->get_pT());
//_h_MS_al_os->Fill(ep,dimuon->get_mass());}

//if(type==SS_ALL_)  {
//_h_pT_al_ss->Fill(ep,dimuon->get_pT());
//_h_MS_al_ss->Fill(ep,dimuon->get_mass());}

//if(type==OS_SIG_N) _h_mass->Fill(ep, dimuon->get_mass());

    }

  return save_event ? EVENT_OK : DISCARDEVENT;
}


bool saModuleDimuonDYDarshana::passesCuts() const{
  // Cuts

  if ( dimuon->get_mass() >= 8 || dimuon->get_mass() <= 4.0 )
    return false;

  if ( !dimuon->get_same_event() ) 
    return false;

  int arm0 = (dimuon->get_Tr0_pz() > 0) ? 1 : 0;
  int arm1 = (dimuon->get_Tr1_pz() > 0) ? 1 : 0;

  if ( arm0 != arm1 ) 
    return false;

  //if ( arm0 != 0 ) 
  //  return false;
      
  if ( dimuon->get_Evt_vtxchi2() >= 4 ) 
    return false;

  //if ( dimuon->get_charge() != 0 ) 
  // return false;

  if ( fabs(dimuon->get_rapidity()) <= 1.2 || fabs(dimuon->get_rapidity()) >= 2.4 )
    return false;

  if ( fabs(dimuon->get_Tr0_pz()) <= 2 || fabs(dimuon->get_Tr1_pz()) <= 2.0 )
    return false;

  if ( dimuon->get_Tr0_DG0() >= 20 || dimuon->get_Tr1_DG0() >= 20 )
    return false;

  if ( dimuon->get_Tr0_DDG0() >= 8 || dimuon->get_Tr1_DDG0() >= 8 )
    return false;

  if ( dimuon->get_Tr0_ntrhits() <= 10 || dimuon->get_Tr1_ntrhits() <= 10 )
    return false;
  if ( dimuon->get_Tr0_nidhits() <= 6 || dimuon->get_Tr1_nidhits() <= 6 )
    return false;
      
  if ( dimuon->get_Tr0_lastgap() <= 3 || dimuon->get_Tr1_lastgap() <= 3 ) 
    return false;

  if ( dimuon->get_Tr0_trchi2() >= 10 || dimuon->get_Tr1_trchi2() >= 10 )
    return false;


  return true;
}


int saModuleDimuonDYDarshana::GetNumberofTracklets(DSTReader *fvtx_trk_map, const DiMuon *dimuon)
{
	if(!fvtx_trk_map){
	cout<<"EXCEPTION: "<<PHWHERE<<endl;
	return NULL;
	}
	
	int ntrklets = 0;
	TClonesArray *array = fvtx_trk_map->get_FvtxCompactTrk();
	for (int i = 0; i < array->GetSize(); i++) {
	TFvtxCompactTrk *tracklet = dynamic_cast<TFvtxCompactTrk*> (array->At(i));
	
	if(!tracklet){
	//cout<<"No tracklet"<<__LINE__<<"size: "<<array->GetSize()<<endl;
	break;
	}
	
	if(_use_cut_tracklet_chi2 && (tracklet->get_chi2_ndf() > _cut_tracklet_chi2)) continue;
	if(!_use_2_hit_tracklet && tracklet->get_nhits() <= 2) continue;

	SingleMuon *muon0 = singlemuoncontainer->get_SingleMuon(0);
	SingleMuon *muon1 = singlemuoncontainer->get_SingleMuon(1);

	float xx0 = tracklet->get_fvtx_vtx().getX()-(tracklet->get_fvtx_vtx().getZ()- dimuons->get_Evt_fvtxZ())*
	tan(tracklet->get_fvtx_theta())*cos(tracklet->get_fvtx_phi());
	float yy0 = tracklet->get_fvtx_vtx().getY()-(tracklet->get_fvtx_vtx().getZ()- dimuons->get_Evt_fvtxZ())*
	tan(tracklet->get_fvtx_theta())*sin(tracklet->get_fvtx_phi());
	float x0y0=sqrt((xx0 - dimuons->get_Evt_fvtxX())*(xx0 - dimuons->get_Evt_fvtxX()) +
	(yy0 - dimuons->get_Evt_fvtxY())*(yy0 - dimuons->get_Evt_fvtxY()));

	if (fabs(TMath::ATan2(sqrt(muon0->get_px_fvtxmutr()*muon0->get_px_fvtxmutr()+
	muon0->get_py_fvtxmutr()*muon0->get_py_fvtxmutr()),muon0->get_pz_fvtxmutr())-
	tracklet->get_fvtx_theta()) >0.001 && fabs(TMath::ATan2(sqrt(muon1->get_px_fvtxmutr()*
	muon1->get_px_fvtxmutr()+muon1->get_py_fvtxmutr()*muon1->get_py_fvtxmutr()),
	muon1->get_pz_fvtxmutr())-tracklet->get_fvtx_theta()) >0.001 && fabs(TMath::ATan2(muon0->get_py_fvtxmutr(),
	muon0->get_px_fvtxmutr())-tracklet->get_fvtx_phi())>0.001 && fabs(TMath::ATan2(muon1->get_py_fvtxmutr(),
	muon1->get_px_fvtxmutr())-tracklet->get_fvtx_phi())>0.001 && tracklet->get_fvtx_theta()+0!=0 && 
	x0y0 < 1.5){
	ntrklets++;
	}

	}

	return ntrklets;
        
}


//! global termination
int saModuleDimuonDYDarshana::end(PHCompositeNode *topNode, sa_hist_mangager_ptr hm)
{

  //calculate asymmetry with relative lumi of BbcNoCutPHENIX
  //_h_mass->CalcAsymmetry(hm.get()->getHisto_DefaultLumi(),saHist::BbcNoCutPHENIX);
  //calculate asymmetry with relative lumi of BbcNoCutPHENIX
  _h_pT_os_N->CalcAsymmetry(hm.get()->getHisto_DefaultLumi(),saHist::BbcVertexCut);
  _h_pT_os_S->CalcAsymmetry(hm.get()->getHisto_DefaultLumi(),saHist::BbcVertexCut);

  _h_pT_sm_N->CalcAsymmetry(hm.get()->getHisto_DefaultLumi(),saHist::BbcVertexCut);
  _h_pT_sm_S->CalcAsymmetry(hm.get()->getHisto_DefaultLumi(),saHist::BbcVertexCut);

//  _h_pT_al_os->CalcAsymmetry(hm.get()->getHisto_DefaultLumi(),saHist::BbcNoCutPHENIX);
//  _h_pT_al_ss->CalcAsymmetry(hm.get()->getHisto_DefaultLumi(),saHist::BbcNoCutPHENIX);


/*  _h_MS_os_N->CalcAsymmetry(hm.get()->getHisto_DefaultLumi(),saHist::BbcVertexCut);
  _h_MS_os_S->CalcAsymmetry(hm.get()->getHisto_DefaultLumi(),saHist::BbcVertexCut);

  _h_MS_sm_N->CalcAsymmetry(hm.get()->getHisto_DefaultLumi(),saHist::BbcVertexCut);
  _h_MS_sm_S->CalcAsymmetry(hm.get()->getHisto_DefaultLumi(),saHist::BbcVertexCut);*/
/*
  _h_MS_os_N->CalcAsymmetry(hm.get()->getHisto_DefaultLumi(),saHist::BbcNoCutPHENIX);
  _h_MS_os_S->CalcAsymmetry(hm.get()->getHisto_DefaultLumi(),saHist::BbcNoCutPHENIX);

  _h_MS_sm_N->CalcAsymmetry(hm.get()->getHisto_DefaultLumi(),saHist::BbcNoCutPHENIX);
  _h_MS_sm_S->CalcAsymmetry(hm.get()->getHisto_DefaultLumi(),saHist::BbcNoCutPHENIX);
*/
 /* _h_MS_al_os->CalcAsymmetry(hm.get()->getHisto_DefaultLumi(),saHist::BbcNoCutPHENIX);
  _h_MS_al_ss->CalcAsymmetry(hm.get()->getHisto_DefaultLumi(),saHist::BbcNoCutPHENIX);

  hm.get()->registerHisto(_h_pT_os_N->CalcAsymmetry_Chi2Fit(false));
  hm.get()->registerHisto(_h_pT_os_S->CalcAsymmetry_Chi2Fit(false));

  hm.get()->registerHisto(_h_pT_sm_N->CalcAsymmetry_Chi2Fit(false));
  hm.get()->registerHisto(_h_pT_sm_S->CalcAsymmetry_Chi2Fit(false));

  hm.get()->registerHisto(_h_pT_al_os->CalcAsymmetry_Chi2Fit(false));
  hm.get()->registerHisto(_h_pT_al_ss->CalcAsymmetry_Chi2Fit(false));

  hm.get()->registerHisto(_h_MS_os_N->CalcAsymmetry_Chi2Fit(false));
  hm.get()->registerHisto(_h_MS_os_S->CalcAsymmetry_Chi2Fit(false));

  hm.get()->registerHisto(_h_MS_sm_N->CalcAsymmetry_Chi2Fit(false));
  hm.get()->registerHisto(_h_MS_sm_S->CalcAsymmetry_Chi2Fit(false));

  hm.get()->registerHisto(_h_MS_al_os->CalcAsymmetry_Chi2Fit(false));
  hm.get()->registerHisto(_h_MS_al_ss->CalcAsymmetry_Chi2Fit(false));
*/
  return EVENT_OK;
}
