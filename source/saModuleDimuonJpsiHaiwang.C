// $Id: saModuleDimuonJpsiHaiwang.C,v 1.4 2014/12/11 22:33:08 yuhw Exp $                                                                                             

/*!
 * \file saModuleDimuonJpsiHaiwang.C
 * \brief modified from Ming's saModuleDimuonMing. for DiMuon Jpsi Analysis 
 * \author Haiwang Yu <yuhw@rcf.rhic.bnl.gov>
 * \version $Revision: 1.4 $
 * \date $Date: 2014/12/11 22:33:08 $
 */
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
#include <RpcDiMuonContainer.h>
#include <RpcDiMuon.h>
#include <TrigLvl1.h>

#include "saFlag.h"

#include "saModuleDimuonJpsiHaiwang.h"


using namespace std;

//-- DiMuon Mass
saModuleDimuonJpsiHaiwang::saModuleDimuonJpsiHaiwang(const std::string &name, bool doControlData) :
	saModuleBase(name),                                     //
	_h_Xing_Check(NULL),
	_h_pT_os_sig_n(NULL), _h_pT_os_sig_s(NULL),
	_h_pT_os_lsb_n(NULL), _h_pT_os_lsb_s(NULL),
	_h_pT_ss_all_n(NULL), _h_pT_ss_all_s(NULL),
	_h_pT_cb_bkg_n(NULL), _h_pT_cb_bkg_s(NULL)
{
	verbosity = 0;

	_use_fixed_mass_win = false;
	_use_bbc_cut = true;
	_use_rpc_cut = true;

	_doControlData = doControlData;
	_FillCounter = 0;

	JPSI_SIG_MASS_MIN_N = 2.7;
	JPSI_SIG_MASS_MAX_N = 3.5;
	JPSI_SIG_MASS_MIN_S = 2.7;
	JPSI_SIG_MASS_MAX_S = 3.5;
	JPSI_LSB_MASS_MIN = 2.00;
	JPSI_LSB_MASS_MAX = 2.50;
	JPSI_ALL_MASS_MIN = 2.00;
	JPSI_ALL_MASS_MAX = 3.60;

	_pT_bins[0] = 0.0;
	_pT_bins[1] = 2.0;
	_pT_bins[2] = 4.0;
	_pT_bins[3] = 10.;

	for(int i=0;i<2;i++)
	{
		for(int j=0;j<3;j++)
		{
			_pT_MinJpsiMs[i][j] = 2.7;
			_pT_MaxJpsiMs[i][j] = 3.5;
		}
	}

	ifs_fitparam.open("fitparam.dat");

	if(ifs_fitparam.is_open()) set_fitparam(2);

	//! for Run 12
	//MIN_FILL_NUM = 16593;
	//MAX_FILL_NUM = 16735;

	//! for Run 12
	MIN_FILL_NUM = 17201;
	MAX_FILL_NUM = 17601;

	fout.open("ControlData.txt");
}

saModuleDimuonJpsiHaiwang::~saModuleDimuonJpsiHaiwang()
{

}

//! global initialization
int saModuleDimuonJpsiHaiwang::init(PHCompositeNode *topNode, sa_hist_mangager_ptr hm)
{
	set_fitparam(2);

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
		cout << "saModuleDimuonJpsiHaiwang::Init failed to create saEventProperty"
			<< endl;
		return ABORTRUN;
	}

	//
	// --- define user histograms ----
	//

	TH2F * h_Xing_Check = new TH2F("Xing_Check","Xing Check", 121,-0.5,120.5,(MAX_FILL_NUM-MIN_FILL_NUM+21),MIN_FILL_NUM-10.5,MAX_FILL_NUM+10.5);
	_h_Xing_Check = new saHist(
			h_Xing_Check,
			saHist::SPIN_DEPENDENT |
			saHist::RUN_PROPERTY | saHist::FILL_PROPERTY
			| saHist::CROSSING_CHECKS, // flags
			Verbosity() // verbosity
			);
	hm.get()->registerHisto(_h_Xing_Check);

	TH1F * h_pT_os_sig_n = new TH1F("pT_os_sig_n","opposite sign, signal region, north", 3, _pT_bins);
	_h_pT_os_sig_n = new saHist(
			h_pT_os_sig_n,
			saHist::SPIN_DEPENDENT |
			saHist::RUN_PROPERTY | saHist::FILL_PROPERTY
			| saHist::CROSSING_CHECKS, // flags
			Verbosity() // verbosity
			);
	hm.get()->registerHisto(_h_pT_os_sig_n);

	TH1F * h_pT_os_sig_s = new TH1F("pT_os_sig_s","opposite sign, signal region, south", 3, _pT_bins);
	_h_pT_os_sig_s = new saHist(
			h_pT_os_sig_s,
			saHist::SPIN_DEPENDENT |
			saHist::RUN_PROPERTY | saHist::FILL_PROPERTY
			| saHist::CROSSING_CHECKS, // flags
			Verbosity() // verbosity
			);
	hm.get()->registerHisto(_h_pT_os_sig_s);

	TH1F * h_pT_os_lsb_n = new TH1F("pT_os_lsb_n","opposite sign, lower sideband, north", 3, _pT_bins);
	_h_pT_os_lsb_n = new saHist(
			h_pT_os_lsb_n,
			saHist::SPIN_DEPENDENT |
			saHist::RUN_PROPERTY | saHist::FILL_PROPERTY
			| saHist::CROSSING_CHECKS, // flags
			Verbosity() // verbosity
			);
	hm.get()->registerHisto(_h_pT_os_lsb_n);

	TH1F * h_pT_os_lsb_s = new TH1F("pT_os_lsb_s","opposite sign, lower sideband, south", 3, _pT_bins);
	_h_pT_os_lsb_s = new saHist(
			h_pT_os_lsb_s,
			saHist::SPIN_DEPENDENT |
			saHist::RUN_PROPERTY | saHist::FILL_PROPERTY
			| saHist::CROSSING_CHECKS, // flags
			Verbosity() // verbosity
			);
	hm.get()->registerHisto(_h_pT_os_lsb_s);

	TH1F * h_pT_ss_all_n = new TH1F("pT_ss_all_n","same sign, all region, north", 3, _pT_bins);
	_h_pT_ss_all_n = new saHist(
			h_pT_ss_all_n,
			saHist::SPIN_DEPENDENT |
			saHist::RUN_PROPERTY | saHist::FILL_PROPERTY
			| saHist::CROSSING_CHECKS, // flags
			Verbosity() // verbosity
			);
	hm.get()->registerHisto(_h_pT_ss_all_n);

	TH1F * h_pT_ss_all_s = new TH1F("pT_ss_all_s","same sign, all region, south", 3, _pT_bins);
	_h_pT_ss_all_s = new saHist(
			h_pT_ss_all_s,
			saHist::SPIN_DEPENDENT |
			saHist::RUN_PROPERTY | saHist::FILL_PROPERTY
			| saHist::CROSSING_CHECKS, // flags
			Verbosity() // verbosity
			);
	hm.get()->registerHisto(_h_pT_ss_all_s);

	TH1F * h_pT_cb_bkg_n = new TH1F("pT_cb_bkg_n","combined bkg, north", 3, _pT_bins);
	_h_pT_cb_bkg_n = new saHist(
			h_pT_cb_bkg_n,
			saHist::SPIN_DEPENDENT |
			saHist::RUN_PROPERTY | saHist::FILL_PROPERTY
			| saHist::CROSSING_CHECKS, // flags
			Verbosity() // verbosity
			);
	hm.get()->registerHisto(_h_pT_cb_bkg_n);

	TH1F * h_pT_cb_bkg_s = new TH1F("pT_cb_bkg_s","combined bkg, south", 3, _pT_bins);
	_h_pT_cb_bkg_s = new saHist(
			h_pT_cb_bkg_s,
			saHist::SPIN_DEPENDENT |
			saHist::RUN_PROPERTY | saHist::FILL_PROPERTY
			| saHist::CROSSING_CHECKS, // flags
			Verbosity() // verbosity
			);
	hm.get()->registerHisto(_h_pT_cb_bkg_s);

	return EVENT_OK;
}

//! Run initialization
int saModuleDimuonJpsiHaiwang::init_run(PHCompositeNode *topNode,
		sa_hist_mangager_ptr hm)
{

	return EVENT_OK;
}

//! event method
int saModuleDimuonJpsiHaiwang::event(PHCompositeNode *topNode, sa_hist_mangager_ptr hm)
{
	DiMuonContainer* dimuons = findNode::getClass<DiMuonContainer>(topNode,
			"DiMuonContainer");
	if (!dimuons)
	{
		cout << "saModuleDimuonJpsiHaiwang:: DiMuonContainer - ERROR - not in Node Tree"
			<< endl;
		return ABORTRUN;
	}

	RpcDiMuonContainer* rpcdimuons = findNode::getClass<RpcDiMuonContainer>(topNode,
			"RpcDiMuonContainer");
	if (!rpcdimuons)
	{
		cout << "saModuleDimuonJpsiHaiwang:: RpcDiMuonContainer - ERROR - not in Node Tree"
			<< endl;
		return ABORTRUN;
	}

	const saEventProperty * ep = findNode::getClass<saEventProperty>(topNode,
			"saEventProperty");
	if (!ep)
	{
		cout << "saModuleDimuonJpsiHaiwang::event - ERROR  - Cannot find EventProperty node in Top Node"
			<< endl;
		return ABORTRUN;
	}

	if (Verbosity() >= 2)
	{
		cout << "saModuleDimuonJpsiHaiwang::event - INFO  - ";
		ep->identify();
	}

	saFlagC *flag = findNode::getClass<saFlagC>(topNode, "DimuonFlag");
	if (!flag)
	{
		cout << "saModuleDimuonJpsiHaiwang::event - DimuonFlag not in Node Tree"
			<< endl;
		return ABORTRUN;
	}

	TrigLvl1* triglv1 = findNode::getClass<TrigLvl1>(topNode,
			"TrigLvl1");
	if (!triglv1)
	{
		cout << "saModuleDimuonJpsiHaiwang:: TrigLvl1 - ERROR - not in Node Tree"
			<< endl;
		return ABORTRUN;
	}

	// syncronizely set the cut flag with dimuon container
	flag->flags.Set(dimuons->get_nDiMuons());
	flag->flags.Reset(0);

	// save the event for DST files?
	// int return_code = dimuons->get_nDiMuons() >0 ? EVENT_OK : DISCARDEVENT;

	bool fill_ok = false;   // save the event for DST output 

	if(_use_bbc_cut) if(!(TMath::Abs(dimuons->get_Evt_bbcZ())<30.)) return DISCARDEVENT;

	for (unsigned int imuon = 0; imuon < dimuons->get_nDiMuons(); imuon++)
	{
		DiMuon* dimuon = dimuons->get_DiMuon(imuon);

		RpcDiMuon* rpcdimuon = rpcdimuons->get_RpcDiMuon(imuon);

		if(!pass_dimuon_cut(dimuon,rpcdimuon)) continue;

		enum EVENT_TYPE type = BAD_EVENT;

		float mass = dimuon->get_mass();

		enum ARM_TYPE arm;
		arm = (dimuon->get_rapidity() > 0) ? NORTH : SOUTH;

		if(!_use_fixed_mass_win)
			for(int i=0;i<3;i++)
				if(dimuon->get_pT()>_pT_bins[i] 
						&& dimuon->get_pT()<=_pT_bins[i+1])
				{
					JPSI_SIG_MASS_MAX_N = _pT_MaxJpsiMs[NORTH][i];
					JPSI_SIG_MASS_MIN_N = _pT_MinJpsiMs[NORTH][i];
					JPSI_SIG_MASS_MAX_S = _pT_MaxJpsiMs[SOUTH][i];
					JPSI_SIG_MASS_MIN_S = _pT_MinJpsiMs[SOUTH][i];
				}

		if (dimuon->get_charge() == 0
				&& arm==NORTH 
				&& mass<JPSI_SIG_MASS_MAX_N && mass>JPSI_SIG_MASS_MIN_N
			 )
			type = OS_SIG_N;

		if (dimuon->get_charge() == 0
				&& arm==SOUTH 
				&& mass<JPSI_SIG_MASS_MAX_S && mass>JPSI_SIG_MASS_MIN_S
			 )
			type = OS_SIG_S;

		if (dimuon->get_charge() == 0
				&& arm==NORTH 
				&& mass<JPSI_LSB_MASS_MAX && mass>JPSI_LSB_MASS_MIN
			 )
			type = OS_LSB_N;

		if (dimuon->get_charge() == 0
				&& arm==SOUTH 
				&& mass<JPSI_LSB_MASS_MAX && mass>JPSI_LSB_MASS_MIN
			 )
			type = OS_LSB_S;

		if (dimuon->get_charge() != 0
				&& arm==NORTH
				//&& mass<JPSI_ALL_MASS_MAX && mass>JPSI_ALL_MASS_MIN
				&& mass<JPSI_SIG_MASS_MAX_N && mass>JPSI_SIG_MASS_MIN_N
			 )
			type = SS_ALL_N;

		if (dimuon->get_charge() != 0
				&& arm==SOUTH 
				//&& mass<JPSI_ALL_MASS_MAX && mass>JPSI_ALL_MASS_MIN
				&& mass<JPSI_SIG_MASS_MAX_S && mass>JPSI_SIG_MASS_MIN_S
			 )
			type = SS_ALL_S;

		if(verbosity>1)
		{
			cout
				<<"DEBUG: "<<__FILE__<<" : "<<__LINE__
				<<" pT: "<<dimuon->get_pT()
				<<" mass: "<<dimuon->get_mass()
				<<" rpcdimuon->get_Tr0_Rpc1St1Time: "<<rpcdimuon->get_Tr0_Rpc1St1Time()
				<<" type: "<<type
				<<endl;
		}

		// this is an event used to fill the histogram 
		flag->flags[imuon] = 1;   

		_h_Xing_Check->Fill(ep,ep->get_crossing_id_RHIC(),ep->get_fill_number());

		if(type==OS_SIG_N) _h_pT_os_sig_n->Fill(ep,dimuon->get_pT());
		if(type==OS_SIG_S) _h_pT_os_sig_s->Fill(ep,dimuon->get_pT());

		if(type==OS_LSB_N) _h_pT_os_lsb_n->Fill(ep,dimuon->get_pT());
		if(type==OS_LSB_S) _h_pT_os_lsb_s->Fill(ep,dimuon->get_pT());

		if(type==SS_ALL_N) _h_pT_ss_all_n->Fill(ep,dimuon->get_pT());
		if(type==SS_ALL_S) _h_pT_ss_all_s->Fill(ep,dimuon->get_pT());

		if(type==OS_LSB_N||type==SS_ALL_N) _h_pT_cb_bkg_n->Fill(ep,dimuon->get_pT());
		if(type==OS_LSB_S||type==SS_ALL_S) _h_pT_cb_bkg_s->Fill(ep,dimuon->get_pT());

		if(_doControlData
				&&dimuon->get_pT()>0.&&dimuon->get_pT()<2.
				//&&ep->get_spin_state()==0
				&&(
					type==OS_SIG_N
					//|| type==OS_SIG_S
					)
				)
		{
			fout
				//<<setw(15)<<" run num: "<<setw(10)<<ep->get_run_number()
				//<<setw(15)<<" event num: "<<setw(10)<<ep->get_event_number()
				//<<setw(15)<<" FillCounter: "<<setw(10)<<_FillCounter++
				//<<setw(15)<<" XingPHENIX: "<<setw(4)<<ep->get_crossing_id()
				//<<setw(15)<<" spin stats: "<<setw(3)<<ep->get_spin_state()
				<<setw(15)<<" beamclk0: "<<setw(10)<<triglv1->get_lvl1_beam_clk(0)
				//<<setw(15)<<" beamclk1: "<<setw(10)<<triglv1->get_lvl1_beam_clk(1)
				//<<setw(15)<<" Tr0_idhits: "<<setw(10)<<dimuon->get_Tr0_idhits()
				//<<setw(15)<<" Tr1_idhits: "<<setw(10)<<dimuon->get_Tr1_idhits()
				<<setw(15)<<" dimuon->get_mass(): "<<setw(10)<<dimuon->get_mass()
				<<setw(15)<<" dimuon->get_pT(): "<<setw(10)<<dimuon->get_pT()
				<<endl;
		}

		fill_ok = true;
	}

	return fill_ok? EVENT_OK : DISCARDEVENT;
}

//! global termination
int saModuleDimuonJpsiHaiwang::end(PHCompositeNode *topNode, sa_hist_mangager_ptr hm)
{

	//calculate asymmetry with relative lumi of BbcVertexCut
	//_h_pT_os_sig_n->CalcAsymmetry(hm.get()->getHisto_DefaultLumi(),saHist::BbcVertexCut);
	//_h_pT_os_sig_s->CalcAsymmetry(hm.get()->getHisto_DefaultLumi(),saHist::BbcVertexCut);

	//_h_pT_os_lsb_n->CalcAsymmetry(hm.get()->getHisto_DefaultLumi(),saHist::BbcVertexCut);
	//_h_pT_os_lsb_s->CalcAsymmetry(hm.get()->getHisto_DefaultLumi(),saHist::BbcVertexCut);

	//_h_pT_ss_all_n->CalcAsymmetry(hm.get()->getHisto_DefaultLumi(),saHist::BbcVertexCut);
	//_h_pT_ss_all_s->CalcAsymmetry(hm.get()->getHisto_DefaultLumi(),saHist::BbcVertexCut);

	//hm.get()->registerHisto(_h_pT_os_sig_n->CalcAsymmetry_Chi2Fit(false));
	//hm.get()->registerHisto(_h_pT_os_sig_s->CalcAsymmetry_Chi2Fit(false));

	//hm.get()->registerHisto(_h_pT_os_lsb_n->CalcAsymmetry_Chi2Fit(false));
	//hm.get()->registerHisto(_h_pT_os_lsb_s->CalcAsymmetry_Chi2Fit(false));

	//hm.get()->registerHisto(_h_pT_ss_all_n->CalcAsymmetry_Chi2Fit(false));
	//hm.get()->registerHisto(_h_pT_ss_all_s->CalcAsymmetry_Chi2Fit(false));

	fout.close();

	return EVENT_OK;
}

//! 
bool saModuleDimuonJpsiHaiwang::pass_dimuon_cut(const DiMuon *dimuon, const RpcDiMuon *rpcdimuon)
{
	if(!(dimuon->get_same_event()==true)) return false;
	if(!(dimuon->get_Evt_vtxchi2()<5.)) return false;
	if(!(dimuon->get_Evt_vtxoor()<1.)) return false;
	if(!(sqrt(dimuon->get_X0()*dimuon->get_X0()+dimuon->get_Y0()*dimuon->get_Y0())<2.)) return false;

	if(!(dimuon->get_Tr0_DG0()<15. && dimuon->get_Tr1_DG0()<15.)) return false;
	if(!(dimuon->get_Tr0_DDG0()<10. && dimuon->get_Tr1_DDG0()<10.)) return false;
	if(!(dimuon->get_Tr0_ntrhits()>9 && dimuon->get_Tr1_ntrhits()>9)) return false;
	if(!(dimuon->get_Tr0_nidhits()>5 && dimuon->get_Tr1_nidhits()>5)) return false;
	if(!(dimuon->get_Tr0_lastgap()>2 && dimuon->get_Tr1_lastgap()>2)) return false;
	if(!(dimuon->get_Tr0_dca_r()<5. && dimuon->get_Tr1_dca_r()<5.)) return false;

	if(!(dimuon->get_Tr0_pz()*dimuon->get_Tr1_pz()>0.)) return false; // both tracks in same arm

	if(!(dimuon->get_pT()<10.)) return false;
	if(!(dimuon->get_pz()<100.)) return false;

	if(!(dimuon->get_mass()>0.5 && dimuon->get_mass()<9.5)) return false;
	if(!(TMath::Abs(dimuon->get_rapidity()) > 1.2 && TMath::Abs(dimuon->get_rapidity()) < 2.2)) return false;

	if(_use_rpc_cut)
	{
		if(!(
					(rpcdimuon->get_Tr0_Rpc1St1Time()>-1 && rpcdimuon->get_Tr0_Rpc1St1Time()<44) || 
					(rpcdimuon->get_Tr0_Rpc3St3Time()>-1 && rpcdimuon->get_Tr0_Rpc3St3Time()<44) ||
					(rpcdimuon->get_Tr1_Rpc1St1Time()>-1 && rpcdimuon->get_Tr1_Rpc1St1Time()<44) || 
					(rpcdimuon->get_Tr1_Rpc3St3Time()>-1 && rpcdimuon->get_Tr1_Rpc3St3Time()<44)
				)) return false;
	}

	return true;
}

bool saModuleDimuonJpsiHaiwang::set_fitparam(int nSigma)
{
	if(!ifs_fitparam.is_open())
	{
		cout<<"ERROR: "<<__FILE__<<" : "<<__LINE__<<endl;
		return false;
	}

	int arm = 0;
	float pTmin = 0;
	float pTmax = 0;
	float r_2sigma = 0;
	float r_3sigma = 0;
	float min_2sigma = 0;
	float max_2sigma = 0;
	float min_3sigma = 0;
	float max_3sigma = 0;

	while(ifs_fitparam
			>>arm>>pTmin>>pTmax
			>>r_2sigma>>r_3sigma
			>>min_2sigma>>max_2sigma
			>>min_3sigma>>max_3sigma)
	{
		if(nSigma == 2)
		{
			for(int i=0;i<3;i++)
				if(TMath::Abs(pTmax-_pT_bins[i+1]) < 0.01)
				{
					_pT_MinJpsiMs[arm][i] = min_2sigma;
					_pT_MaxJpsiMs[arm][i] = max_2sigma;

					cout
						<<arm<<" "<<i<<" "
						<<_pT_MinJpsiMs[arm][i]<<" "<<_pT_MaxJpsiMs[arm][i]<<endl;
				}
		}
		if(nSigma == 3)
		{
			for(int i=0;i<3;i++)
				if(TMath::Abs(pTmax-_pT_bins[i+1]) < 0.01)
				{
					_pT_MinJpsiMs[arm][i] = min_3sigma;
					_pT_MaxJpsiMs[arm][i] = max_3sigma;
				}
		}
	}

	return true;
}

