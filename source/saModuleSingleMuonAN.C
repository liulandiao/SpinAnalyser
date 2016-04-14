/*!
 * \file saModuleSingleMuonAN.C
 * \brief modified from Haiwang for SingleMuon Analysis
 * \author Landiao <liulandiao@gmail.com>
 * \version $Revision: 1.0 $
 * \date $Date: 2016/4/13  $
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
#include <SingleMuonContainer.h>
#include <SingleMuon.h>
#include <TrigLvl1.h>

#include "saFlag.h"
#include "saModuleSingleMuonAN.h"

//Black Magic
#define __ERROR__ "ERROR: "<<__FILE__<<": "<<__LINE__<<endl

#define SMART(expr) boost::shared_ptr<expr>
#define NEW(expr) boost::make_shared<expr>

using namespace std;

//-- SingleMuon Mass
saModuleSingleMuonAN::saModuleSingleMuonAN(const std::string &name, bool doControlData) :
	saModuleBase(name),                                     //
	_sah_pT_phi_MuM_phi_N_LG4(NULL), _sah_pT_phi_MuM_phi_N_LG23(NULL),
	_sah_pT_phi_MuM_phi_S_LG4(NULL), _sah_pT_phi_MuM_phi_S_LG23(NULL)
{
	verbosity = 0;

	_use_bbc_cut = true;

	_doControlData = doControlData;
	_FillCounter = 0;

	_pT_bins[0] = 0.0;
	_pT_bins[1] = 2.0;
	_pT_bins[2] = 4.0;
	_pT_bins[3] = 10.;

	fout.open("ControlData.txt");
}

saModuleSingleMuonAN::~saModuleSingleMuonAN()
{

}

//! global initialization
int saModuleSingleMuonAN::init(PHCompositeNode *topNode, sa_hist_mangager_ptr hm)
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
				"SinglemuonFlag", "PHObject");

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
		cout << "saModuleSingleMuonAN::Init failed to create saEventProperty"
			<< endl;
		return ABORTRUN;
	}

	//
	// --- define user histograms ----
	//
	_v_sahist.clear();

	_sah_pT_phi_MuM_phi_N_LG4  = makesaHist_pT_phi(hm,makeTHist_pT_phi("_sah_pT_phi_MuM_phi_N_LG4", "North, Last Gap = 4"));
        _sah_pT_phi_MuM_phi_N_LG23 = makesaHist_pT_phi(hm,makeTHist_pT_phi("_sah_pT_phi_MuM_phi_N_LG23","North, Last Gap = 23"));
        _sah_pT_phi_MuM_phi_S_LG4  = makesaHist_pT_phi(hm,makeTHist_pT_phi("_sah_pT_phi_MuM_phi_S_LG4", "South, Last Gap = 4"));
        _sah_pT_phi_MuM_phi_S_LG23 = makesaHist_pT_phi(hm,makeTHist_pT_phi("_sah_pT_phi_MuM_phi_S_LG23","South, Last Gap = 23"));
	return EVENT_OK;
}

//! Run initialization
int saModuleSingleMuonAN::init_run(PHCompositeNode *topNode,
		sa_hist_mangager_ptr hm)
{
	return EVENT_OK;
}

//! event method
int saModuleSingleMuonAN::event(PHCompositeNode *topNode, sa_hist_mangager_ptr hm)
{
	SingleMuonContainer* sngmuons = findNode::getClass<SingleMuonContainer>(topNode,
			"SingleMuonContainer");
	if (!sngmuons)
	{
		cout << "saModuleSingleMuonAN:: SingleMuonContainer - ERROR - not in Node Tree"
			<< endl;
		return ABORTRUN;
	}

	const saEventProperty * ep = findNode::getClass<saEventProperty>(topNode,
			"saEventProperty");
	if (!ep)
	{
		cout << "saModuleSingleMuonAN::event - ERROR  - Cannot find EventProperty node in Top Node"
			<< endl;
		return ABORTRUN;
	}

	if (Verbosity() >= 2)
	{
		cout << "saModuleSingleMuonAN::event - INFO  - ";
		ep->identify();
	}

	saFlagC *flag = findNode::getClass<saFlagC>(topNode, "SinglemuonFlag");
	if (!flag)
	{
		cout << "saModuleSingleMuonAN::event - SinglemuonFlag not in Node Tree"
			<< endl;
		return ABORTRUN;
	}

	TrigLvl1* triglv1 = findNode::getClass<TrigLvl1>(topNode,
			"TrigLvl1");
	if (!triglv1)
	{
		cout << "saModuleSingleMuonAN:: TrigLvl1 - ERROR - not in Node Tree"
			<< endl;
		return ABORTRUN;
	}

	// syncronizely set the cut flag with singlemuon container
	flag->flags.Set(sngmuons->get_nSingleMuons());
	flag->flags.Reset(0);

	// save the event for DST files?
	// int return_code = sngmuons->get_nSingleMuons() >0 ? EVENT_OK : DISCARDEVENT;

	bool fill_ok = false;   // save the event for DST output 

	if(_use_bbc_cut) if(!(TMath::Abs(sngmuons->get_Evt_bbcZ())<30.)) return DISCARDEVENT;

	for (unsigned int imuon = 0; imuon < sngmuons->get_nSingleMuons(); imuon++)
	{
		SingleMuon* singlemuon = sngmuons->get_SingleMuon(imuon);

		if(!pass_sngmuon_cut(singlemuon)) continue;

		enum EVENT_TYPE type = BAD_EVENT;

		int last_gap = singlemuon->get_lastgap();

		enum ARM_TYPE arm;
		arm = (singlemuon->get_rapidity() > 0) ? NORTH : SOUTH;

		if (singlemuon->get_charge() == 0
                                && arm==NORTH
                                && last_gap==4
                         )
                        type = MuM_N_LG4;

                if (singlemuon->get_charge() == 0
                                && arm==NORTH
                                && (last_gap==2 || last_gap==3)
                         )
                        type = MuM_N_LG23;

		if (singlemuon->get_charge() == 0
				&& arm==SOUTH
				&& last_gap==4
			 )
			type = MuM_S_LG4;

                if (singlemuon->get_charge() == 0
                                && arm==SOUTH
                                && (last_gap==2 || last_gap==3)
                         )
                        type = MuM_S_LG23;

		if(verbosity>1)
		{
			cout
				<<" pT: "<<singlemuon->get_pT()
				<<" type: "<<type
				<<endl;
		}

		// this is an event used to fill the histogram 
		flag->flags[imuon] = 1;   

		float py = singlemuon->get_py();
		float px = singlemuon->get_px();
		double phi = TMath::ATan(py/px)+(px<0&&py>0)*TMath::Pi()-(px<0&&py<0)*TMath::Pi();


		if(type==MuM_N_LG4)  _sah_pT_phi_MuM_phi_N_LG4->Fill( ep,singlemuon->get_pT(),phi);
		if(type==MuM_N_LG23) _sah_pT_phi_MuM_phi_N_LG23->Fill(ep,singlemuon->get_pT(),phi);
		if(type==MuM_S_LG4)  _sah_pT_phi_MuM_phi_S_LG4->Fill( ep,singlemuon->get_pT(),phi);
		if(type==MuM_S_LG23) _sah_pT_phi_MuM_phi_S_LG23->Fill(ep,singlemuon->get_pT(),phi);


		if(_doControlData
				&&singlemuon->get_pT()>0.&&singlemuon->get_pT()<2.
				//&&ep->get_spin_state()==0
				&&(
					type==MuM_N_LG4
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
				//<<setw(15)<<" Tr0_idhits: "<<setw(10)<<singlemuon->get_Tr0_idhits()
				//<<setw(15)<<" Tr1_idhits: "<<setw(10)<<singlemuon->get_Tr1_idhits()
				<<endl;

		}

		fill_ok = true;
	}

	return fill_ok? EVENT_OK : DISCARDEVENT;
}

//! global termination
int saModuleSingleMuonAN::end(PHCompositeNode *topNode, sa_hist_mangager_ptr hm)
{

	fout.close();

	BOOST_FOREACH( saHist* h, _v_sahist )
	{
		h->CalcAsymmetry(hm.get()->getHisto_DefaultLumi(),saHist::BbcVertexCut);
		hm.get()->registerHisto(h->CalcAsymmetry_Chi2Fit(false));
	}

	_v_sahist.clear();

	return EVENT_OK;
}

//! 
bool saModuleSingleMuonAN::pass_sngmuon_cut(const SingleMuon *singlemuon)
{
	return true;
}

TH1 *saModuleSingleMuonAN::makeTHist_pT_phi(
		std::string name,
		std::string title
		)
{
	TH1 * h = new TH2F(name.data(),title.data(), 3, _pT_bins,100, -1.*TMath::Pi(), TMath::Pi());
	return h;
}

saHist *saModuleSingleMuonAN::makesaHist_pT_phi(
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

