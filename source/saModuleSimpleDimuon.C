// $Id: saModuleSimpleDimuon.C,v 1.5 2013/08/14 07:07:59 jinhuang Exp $                                                                                             

/*!
 * \file saModuleSimpleDimuon.C
 * \brief 
 * \author Jin Huang <jhuang@bnl.gov>
 * \version $Revision: 1.5 $
 * \date $Date: 2013/08/14 07:07:59 $
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

#include "saFlag.h"

#include <Fun4AllReturnCodes.h>
#include <PHNodeIterator.h>
#include <getClass.h>

#include "saModuleSimpleDimuon.h"

using namespace std;

saModuleSimpleDimuon::saModuleSimpleDimuon(const std::string &name) :
    saModuleBase(name), _h_mass(NULL)
{

}

saModuleSimpleDimuon::~saModuleSimpleDimuon()
{

}

//! global initialization
int
saModuleSimpleDimuon::init(PHCompositeNode *topNode, sa_hist_mangager_ptr hm)
{

  // make a simple histogram
  TH1F * h_mass = new TH1F("InvMass", "Invariant Mass;Invariant Mass (GeV)", 20,
      1, 10);
  _h_mass = new saHist(
      //
      h_mass, //the template histogram, any TH1 and derivatives is accepted
      saHist::DEFAULT_FLAG | saHist::RUN_PROPERTY | saHist::FILL_PROPERTY
          | saHist::CROSSING_CHECKS, // flags
      Verbosity() // verbosity
      );
  hm.get()->registerHisto(_h_mass);

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
          "SimpleDimuonFlag", "PHObject");

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
      cout << "saModuleSimpleDimuon::Init failed to create saEventProperty"
          << endl;
      return ABORTRUN;
    }

  return EVENT_OK;
}

//! Run initialization
int
saModuleSimpleDimuon::init_run(PHCompositeNode *topNode,
    sa_hist_mangager_ptr hm)
{

  return EVENT_OK;
}

//! event method
int
saModuleSimpleDimuon::event(PHCompositeNode *topNode, sa_hist_mangager_ptr hm)
{

  DiMuonContainer* dimuons = findNode::getClass<DiMuonContainer>(topNode,
      "DiMuonContainer");
  if (!dimuons)
    {
      cout
          << "saModuleSimpleDimuon:: DiMuonContainer - ERROR - not in Node Tree"
          << endl;
      return ABORTRUN;
    }

  const saEventProperty * ep = findNode::getClass<saEventProperty>(topNode,
      "saEventProperty");
  if (!ep)
    {
      cout
          << "saModuleSimpleDimuon::event - ERROR  - Cannot find EventProperty node in Top Node"
          << endl;
      return ABORTRUN;
    }

  if (Verbosity() >= 2)
    {
      cout << "saModuleSimpleDimuon::event - INFO  - ";
      ep->identify();
    }

  saFlagC *flag = findNode::getClass<saFlagC>(topNode, "SimpleDimuonFlag");
  if (!flag)
    {
      cout << "saModuleSimpleDimuon::event - SimpleDimuonFlag not in Node Tree"
          << endl;
      return ABORTRUN;
    }

  // syncronizely set the cut flag with dimuon container
  flag->flags.Set(dimuons->get_nDiMuons());
  flag->flags.Reset(0);

  // save the event for DST files?
  int return_code = dimuons->get_nDiMuons() >0 ? EVENT_OK : DISCARDEVENT;

  for (unsigned int imuon = 0; imuon < dimuons->get_nDiMuons(); imuon++)
    {
      DiMuon* dimuon = dimuons->get_DiMuon(imuon);

      int type = 0;
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

      if (type == 0 and arm0 == arm1)
        {
          flag->flags[imuon] = 1;
          _h_mass->Fill(ep, dimuon->get_mass());
        }
    }

  return return_code;
}

//! global termination
int
saModuleSimpleDimuon::end(PHCompositeNode *topNode, sa_hist_mangager_ptr hm)
{

  //calculate asymmetry with relative lumi of BbcNoCutPHENIX
  _h_mass->CalcAsymmetry(hm.get()->getHisto_DefaultLumi(),
      saHist::BbcNoCutPHENIX);

  return EVENT_OK;
}

