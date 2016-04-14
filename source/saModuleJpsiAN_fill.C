// $Id: saModuleJpsiAN.C,v 1.4 2015/09/03 20:45:38 yuhw Exp $                                                                                             

/*!
 * \file saModuleJpsiAN.C
 * \brief modified from Ming's saModuleDimuonMing. for DiMuon Jpsi Analysis 
 * \author Haiwang Yu <yuhw@rcf.rhic.bnl.gov>
 * \version $Revision: 1.4 $
 * \date $Date: 2015/09/03 20:45:38 $
 */

//STL

//BOOST
#include<boost/make_shared.hpp>
#include<boost/foreach.hpp>

//ROOT
#include <TH1F.h>
#include <TH3.h>
#include <TH2.h>
#include <TF1.h>
#include <TMath.h>

//PHENIX
#include <Fun4AllReturnCodes.h>
#include <PHCompositeNode.h>
#include <PHNodeIterator.h>
#include <getClass.h>
#include <MWGConsts.h>
#include <Tools.h>
#include <Fun4AllServer.h>
#include <Fun4AllHistoManager.h>
#include <TrigRunLvl1.h>
#include <SyncObject.h>
#include <RunHeader.h>
#include <DiMuonContainer.h>
#include <DiMuon.h>
#include <MCDiMuonContainer.h>
#include <MCDiMuon.h>
#include <TrigLvl1.h>

#include "saFlag.h"
#include "saModuleJpsiAN.h"

//Black Magic
#define __ERROR__ "ERROR: "<<__FILE__<<": "<<__LINE__<<endl

#define SMART(expr) boost::shared_ptr<expr>
#define NEW(expr) boost::make_shared<expr>

using namespace std;

//-- DiMuon Mass
saModuleJpsiAN::saModuleJpsiAN(const std::string &name, bool doControlData) :
	saModuleBase(name),                                     //
	_sah_pT_os_sig_n(NULL), _sah_pT_os_sig_s(NULL),
	_sah_pT_os_lsb_n(NULL), _sah_pT_os_lsb_s(NULL)
{
	verbosity = 0;

	_use_fixed_mass_win = true;
	_use_bbc_cut = true;

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

	//ifs_fitparam.open("fitparam.dat");

	//if(ifs_fitparam.is_open()) set_fitparam(2);

	//! for Run 12
	//MIN_FILL_NUM = 16593;
	//MAX_FILL_NUM = 16735;

	//! for Run 13
	//MIN_FILL_NUM = 17201;
	//MAX_FILL_NUM = 17601;

	fout.open("ControlData.txt");
}

saModuleJpsiAN::~saModuleJpsiAN()
{

}

//! global initialization
int saModuleJpsiAN::init(PHCompositeNode *topNode, sa_hist_mangager_ptr hm)
{
	//set_fitparam(2);

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
		cout << "saModuleJpsiAN::Init failed to create saEventProperty"
			<< endl;
		return ABORTRUN;
	}

	//
	// --- define user histograms ----
	//
	_v_sahist.clear();

	_sah_pT_os_sig_n = makesaHist_pT_phi(hm,makeTHist_pT_phi("_sah_pT_os_sig_n","opposite sign, signal region, north"));
	_sah_pT_os_sig_s = makesaHist_pT_phi(hm,makeTHist_pT_phi("_sah_pT_os_sig_s","opposite sign, signal region, south"));
	_sah_pT_os_lsb_n = makesaHist_pT_phi(hm,makeTHist_pT_phi("_sah_pT_os_lsb_n","same sign, lower side band, north"));
	_sah_pT_os_lsb_s = makesaHist_pT_phi(hm,makeTHist_pT_phi("_sah_pT_os_lsb_s","same sign, lower side band, south"));

	return EVENT_OK;
}

//! Run initialization
int saModuleJpsiAN::init_run(PHCompositeNode *topNode,
		sa_hist_mangager_ptr hm)
{
	return EVENT_OK;
}

//! event method
int saModuleJpsiAN::event(PHCompositeNode *topNode, sa_hist_mangager_ptr hm)
{
	DiMuonContainer* dimuons = findNode::getClass<DiMuonContainer>(topNode,
			"DiMuonContainer");
	if (!dimuons)
	{
		cout << "saModuleJpsiAN:: DiMuonContainer - ERROR - not in Node Tree"
			<< endl;
		return ABORTRUN;
	}

	const saEventProperty * ep = findNode::getClass<saEventProperty>(topNode,
			"saEventProperty");
	if (!ep)
	{
		cout << "saModuleJpsiAN::event - ERROR  - Cannot find EventProperty node in Top Node"
			<< endl;
		return ABORTRUN;
	}

	if (Verbosity() >= 2)
	{
		cout << "saModuleJpsiAN::event - INFO  - ";
		ep->identify();
	}

	saFlagC *flag = findNode::getClass<saFlagC>(topNode, "DimuonFlag");
	if (!flag)
	{
		cout << "saModuleJpsiAN::event - DimuonFlag not in Node Tree"
			<< endl;
		return ABORTRUN;
	}

	TrigLvl1* triglv1 = findNode::getClass<TrigLvl1>(topNode,
			"TrigLvl1");
	if (!triglv1)
	{
		cout << "saModuleJpsiAN:: TrigLvl1 - ERROR - not in Node Tree"
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

		if(!pass_dimuon_cut(dimuon)) continue;

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
				<<" pT: "<<dimuon->get_pT()
				<<" mass: "<<dimuon->get_mass()
				<<" type: "<<type
				<<endl;
		}

		// this is an event used to fill the histogram 
		flag->flags[imuon] = 1;   

		float py = dimuon->get_py();
		float px = dimuon->get_px();
		double phi = TMath::ATan(py/px)+(px<0&&py>0)*TMath::Pi()-(px<0&&py<0)*TMath::Pi();


		if(type==OS_SIG_N) _sah_pT_os_sig_n->Fill(ep,dimuon->get_pT(),phi);
		if(type==OS_SIG_S) _sah_pT_os_sig_s->Fill(ep,dimuon->get_pT(),phi);

		if(type==OS_LSB_N) _sah_pT_os_lsb_n->Fill(ep,dimuon->get_pT(),phi);
		if(type==OS_LSB_S) _sah_pT_os_lsb_s->Fill(ep,dimuon->get_pT(),phi);

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
int saModuleJpsiAN::end(PHCompositeNode *topNode, sa_hist_mangager_ptr hm)
{

	fout.close();

	BOOST_FOREACH(saHist *h, _v_sahist)
	{
		h->CalcAsymmetry(hm.get()->getHisto_DefaultLumi(),saHist::BbcVertexCut);
		hm.get()->registerHisto(h->CalcAsymmetry_Chi2Fit(false));
	}

	_v_sahist.clear();

	return EVENT_OK;
}

//! 
bool saModuleJpsiAN::pass_dimuon_cut(const DiMuon *dimuon)
{
	if(!(dimuon->get_same_event()==true)) return false;
	if(!(dimuon->get_Evt_vtxchi2()<5.)) return false;
	//if(!(dimuon->get_Evt_vtxoor()<1.)) return false;
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

	return true;
}

bool saModuleJpsiAN::set_fitparam(int nSigma)
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

TH1 *saModuleJpsiAN::makeTHist_pT_phi(
		std::string name,
		std::string title
		)
{
	TH1 * h = new TH2F(name.data(),title.data(), 3, _pT_bins,100, -1.*TMath::Pi(), TMath::Pi());
	return h;
}

saHist *saModuleJpsiAN::makesaHist_pT_phi(
		sa_hist_mangager_ptr hm,
		TH1 *h
		)
{
	//TH2F * h = new TH2F(name.data(),title.data(), 3, _pT_bins,100, -1.*TMath::Pi(), TMath::Pi());
	//SMART(TH2F) h = NEW(TH2F)(name.data(),title.data(), 3, _pT_bins,100, -1.*TMath::Pi(), TMath::Pi());
	saHist *sahist = new saHist(
			h,
			saHist::SPIN_DEPENDENT |
			saHist::RUN_PROPERTY | saHist::FILL_PROPERTY
			| saHist::CROSSING_CHECKS, // flags
			Verbosity() // verbosity
			);
	hm.get()->registerHisto(sahist);

	_v_sahist.push_back(sahist);

	return sahist;
}

