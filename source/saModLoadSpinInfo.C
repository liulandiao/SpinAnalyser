// $Id: saModLoadSpinInfo.C,v 1.10 2013/08/08 05:34:27 jinhuang Exp $                                                                                             

/*!
 * \file saModLoadSpinInfo.C
 * \brief 
 * \author Jin Huang <jhuang@bnl.gov>
 * \version $Revision: 1.10 $
 * \date $Date: 2013/08/08 05:34:27 $
 */

#include <fstream>
#include <iostream>     // std::cout
#include <sstream>      // std::istringstream
#include <string>       // std::string
#include <TH1D.h>
#include <TH2D.h>

#include <SpinDBContent.hh>
#include <SpinDBOutput.hh>
#include <EventHeader.h>
#include <RunHeader.h>
#include <SyncObject.h>
#include <TMutNode.h>
#include <recoConsts.h>
#include <TrigLvl1.h>

#include <Fun4AllReturnCodes.h>
#include <PHNodeIterator.h>
#include <getClass.h>

#include "saEventProperty_v1.h"
typedef saEventProperty_v1 saEventProperty_t;

#include "saModLoadSpinInfo.h"

using namespace std;

saModLoadSpinInfo::saModLoadSpinInfo(const std::string &name) :
    saModuleBase(name), _run_header(NULL), _event_header(NULL), _sync(NULL), _trig_lv1(
        NULL), _run_num(-1), //
    _spin_cont(static_cast<SpinDBContent *>(NULL)), //
    _h_lumi(NULL), _use_spin_QA(false)
{

}

saModLoadSpinInfo::~saModLoadSpinInfo()
{
  // TODO Auto-generated destructor stub
}

//! global initialization
int
saModLoadSpinInfo::init(PHCompositeNode *topNode, sa_hist_mangager_ptr hm)
{

  PHCompositeNode* dstNode = NULL;
  PHNodeIterator nodeItr(topNode);
  dstNode = static_cast<PHCompositeNode*>(nodeItr.findFirst("PHCompositeNode",
      "DST"));
  if (!dstNode)
    {
      dstNode = new PHCompositeNode("DST");
      topNode->addNode(dstNode);
    }

  saEventProperty *ep = new saEventProperty_t();
  if (ep)
    {
      PHIODataNode<PHObject> * node = new PHIODataNode<PHObject>(ep,
          "saEventProperty", "PHObject");

      if (!node)
        {

          cout << "saSpinAnalyzer::Init failed to create saEventProperty Node"
              << endl;
          return ABORTRUN;

        }
      else
        {

          dstNode->addNode(node);
          cout << "saEventProperty Node is added with version "
              << ep->ClassName() << endl;
        }
    }
  else
    {
      cout << "saSpinAnalyzer::Init failed to create saEventProperty" << endl;
      return ABORTRUN;
    }

  TH1D * h_lumi = new TH1D("RelLumi",
      "Relative Luminorsity in array of enu_lumi_id from saModLoadSpinInfo", //
      saHist::LUMI_COUNT, .5, saHist::LUMI_COUNT + .5);
  _h_lumi = new saHist(h_lumi, saHist::DEFAULT_FLAG_LUMI, Verbosity() - 1);
  hm.get()->registerHisto_DefaultLumi(_h_lumi);
//  _h_lumi->Verbosity(Verbosity() + 1);

  return EVENT_OK;
}

//! Run initialization
int
saModLoadSpinInfo::init_run(PHCompositeNode *top_node, sa_hist_mangager_ptr hm)
{

  try
    {
      _event_header = TMutNode<EventHeader>::find_io_node(top_node,
          "EventHeader");
    }
  catch (...)
    {
      _event_header = NULL;
    }

  try
    {
      _sync = TMutNode<SyncObject>::find_io_node(top_node, "Sync");
    }
  catch (...)
    {
      _sync = NULL;
    }

  try
    {
      _trig_lv1 = TMutNode<TrigLvl1>::find_io_node(top_node, "TrigLvl1");
    }
  catch (...)
    {
      _trig_lv1 = NULL;
    }

  try
    {
      _run_header = TMutNode<RunHeader>::find_io_node(top_node, "RunHeader");
    }
  catch (...)
    {
      _run_header = NULL;
    }

  if (_run_header)
    {
      _run_num = _run_header->get_RunNumber();

      if (Verbosity() >= 1)
        {
          cout << "saModLoadSpinInfo::init_run - from RunHeader : run number = "
              << _run_num << endl;
        }

    }
  else
    {
      recoConsts *rc = recoConsts::instance();

      _run_num = rc->get_IntFlag("RUNNUMBER");

      if (Verbosity() >= 1)
        {
          cout
              << "saModLoadSpinInfo::init_run - from recoConsts : run number = "
              << _run_num << endl;
        }

    }

  if (Verbosity() >= 1)
    {
      cout << "saModLoadSpinInfo::init_run - loading spin info from DB for run "
          << _run_num << endl;
    }

  _spin_cont = boost::make_shared<SpinDBContent>();

  SpinDBOutput spin_out("phnxrc");

  int ret = 0;

  bool bad_run = false;

  ret = spin_out.StoreDBContent(_run_num, _run_num);
  if (ret <= 0)
    {
      cout
          << "saModLoadSpinInfo::init_run - Error loading SpinDBOutput::StoreDBContent for run "
          << _run_num << " with return = " << ret << endl;

      bad_run = true;
    }

  ret = spin_out.GetDBContentStore(*(_spin_cont.get()), _run_num);
  if (ret <= 0)
    {
      cout
          << "saModLoadSpinInfo::init_run - Error loading SpinDBOutput::GetDBContentStore for run "
          << _run_num << " with return = " << ret << endl;

      bad_run = true;
    }

  if (bad_run)
    {
      cout
                << "saModLoadSpinInfo::init_run - Error - this is an invalid run "
                << _run_num << " in database. Set its spin info to be invalid"<< endl;

      _spin_cont->SetFillNumber(0);
      _spin_cont->SetBadRunFlag(1);
      _spin_cont->SetCrossingShift(0);

      for (int crossing = 0; crossing < saEventProperty::N_CROSSING; crossing++)
        {
          _spin_cont->SetBadBunchFlag(crossing, 1);
        }

    }

  if (Verbosity() >= 1)
    {
      _spin_cont->Print();
    }

  // load per-run property
  _pol_blue.assign(-1);
  _pol_stat_err_blue.assign(-1);
  _pol_sys_err_blue.assign(-1);
  _pol_yellow.assign(-1);
  _pol_stat_err_yellow.assign(-1);
  _pol_sys_err_yellow.assign(-1);

  saEventProperty_t ep_tmp;

  ep_tmp.set_fill_number(_spin_cont->GetFillNumber());

  ep_tmp.set_run_number(_run_num);

  vector<double> lumi_sum;
  lumi_sum.assign(saEventProperty::SUM, 0);

  if (Verbosity() >= 2)
    {
      cout << "saModLoadSpinInfo::init_run - _h_lumi->Fill for run " << _run_num
          << " - ";
      ep_tmp.identify();
    }

  if (_use_spin_QA)
    {
      cout << "saModLoadSpinInfo::init_run - Use custom spin QA for run "
          << _run_num << endl;

      if (!_spin_QA.have_run_crossing_flag(_run_num))
        {

          cout << "saModLoadSpinInfo::init_run - ERROR - Cannot find run "
              << _run_num
              << " in the custom QA list. Will assign unknown spin state for all events!!"
              << endl;

        }

    }
  else
    cout << "saModLoadSpinInfo::init_run - DO NOT use custom spin QA for run "
        << _run_num << ".";

  bool fill_lumi = false;
  saHist * _h_lumi_run = _h_lumi->getsaHist_run(
      static_cast<unsigned int>(_run_num), false);
  if (!_h_lumi_run)
    {
      fill_lumi = true;
      if (Verbosity() >= 1)
        {
          cout << "saModLoadSpinInfo::init_run - lumi information for run "
              << _run_num
              << " has not been saved yet. Will fill it in from database."
              << endl;
        }
    }
  else
    {
      fill_lumi = false;
      if (Verbosity() >= 1)
        {
          cout << "saModLoadSpinInfo::init_run - lumi information for run "
              << _run_num << " has already been filled. Will not fill again. "
              << "Dump of its existing data:" << endl;
          _h_lumi_run->Print();
        }
    }

  for (int crossing = 0; crossing < saEventProperty::N_CROSSING; crossing++)
    {
      _spin_cont->GetPolarizationBluePHENIX(crossing, _pol_blue[crossing],
          _pol_stat_err_blue[crossing], _pol_sys_err_blue[crossing]);
      _spin_cont->GetPolarizationYellowPHENIX(crossing, _pol_yellow[crossing],
          _pol_stat_err_yellow[crossing], _pol_sys_err_yellow[crossing]);

      saEventProperty::sa_spin_state state = get_spin_state(crossing);

      ep_tmp.set_crossing_id(crossing);

      ep_tmp.set_crossing_id_RHIC(
          (crossing + _spin_cont->GetCrossingShift())
              % saEventProperty::N_CROSSING);

      ep_tmp.set_spin_state(state);

      ep_tmp.set_polarization_blue(_pol_blue[crossing]);

      ep_tmp.set_polarization_yellow(_pol_yellow[crossing]);

      ep_tmp.set_polarization_stat_err_blue(_pol_stat_err_blue[crossing]);

      ep_tmp.set_polarization_stat_err_yellow(_pol_stat_err_yellow[crossing]);

      ep_tmp.set_polarization_sys_err_blue(_pol_sys_err_blue[crossing]);

      ep_tmp.set_polarization_sys_err_yellow(_pol_sys_err_yellow[crossing]);

      for (int lumi_id = (int) saHist::BbcVertexCut;
          lumi_id < (int) saHist::LUMI_COUNT; lumi_id++)
        {
          double lumi = 0;

          switch (lumi_id)
            {
          case saHist::BbcVertexCut:
            lumi = _spin_cont->GetScalerBbcVertexCutPHENIX(crossing);
            break;

          case saHist::BbcNoCutPHENIX:
            lumi = _spin_cont->GetScalerBbcNoCutPHENIX(crossing);
            break;

          case saHist::ZdcWidePHENIX:
            lumi = _spin_cont->GetScalerZdcWidePHENIX(crossing);
            break;

          case saHist::ScalerZdcNarrow:
            lumi = _spin_cont->GetScalerZdcNarrowPHENIX(crossing);
            break;

          case saHist::ScalerFVTX:
            // not impelmented right now
            break;
          default:
            break;
            }

          if (Verbosity() >= 2)
            {
              cout << "saModLoadSpinInfo::init_run - _h_lumi->Fill with ";
              ep_tmp.identify();
            }

          if (fill_lumi)
            {
              _h_lumi->Fill(&ep_tmp, lumi_id, lumi);
            }

          if (Verbosity() >= 1)
            {
              if (lumi_id == (int) saHist::BbcVertexCut)
                {
                  lumi_sum[state] += lumi;
                }

            }

        }

    } //   for (int crossing = 0; crossing < saEventProperty::N_CROSSING; crossing++)

  if (Verbosity() >= 1)
    {
      cout << "saModLoadSpinInfo::init_run - loaded lumi info for run "
          << _run_num << ", BbcVertexCut :" << endl << " PP = "
          << lumi_sum[saEventProperty::PP] << endl << " PM = "
          << lumi_sum[saEventProperty::PM] << endl << " MP = "
          << lumi_sum[saEventProperty::MP] << endl << " MM = "
          << lumi_sum[saEventProperty::MM] << endl << " UNKNOWN = "
          << lumi_sum[saEventProperty::UNKNOWN] << endl;

    }

  return EVENT_OK;
}

//! event method
int
saModLoadSpinInfo::event(PHCompositeNode *top_node, sa_hist_mangager_ptr hm)
{

  int event_num = -1;
  if (_event_header)
    event_num = _event_header->get_EvtSequence();
//  else if (_sync)
//    event_num = _sync -> EventNumber();

  int crossing = -1;
  if (_trig_lv1)
    crossing = _trig_lv1->get_lvl1_clock_cross();

  if (crossing < 0 || crossing >= saEventProperty::N_CROSSING)
    {
      cout << "saModLoadSpinInfo::event - Error - Invalid crossing id "
          << crossing << ", abort event" << endl;
      return ABORTEVENT;
    }

  saEventProperty *ep = findNode::getClass<saEventProperty>(top_node,
      "saEventProperty");
  if (!ep)
    {
      cout << "saModLoadSpinInfo::event - saEventProperty not in Node Tree"
          << endl;
      return ABORTRUN;
    }
  ep->Reset();

  ep->set_fill_number(_spin_cont->GetFillNumber());

  ep->set_run_number(_run_num);

  ep->set_event_number(event_num);

  ep->set_crossing_id(crossing);

  ep->set_crossing_id_RHIC(
      (crossing + _spin_cont->GetCrossingShift())
          % saEventProperty::N_CROSSING);

  ep->set_spin_state(get_spin_state(crossing));

  ep->set_polarization_blue(_pol_blue[crossing]);

  ep->set_polarization_yellow(_pol_yellow[crossing]);

  ep->set_polarization_stat_err_blue(_pol_stat_err_blue[crossing]);

  ep->set_polarization_stat_err_yellow(_pol_stat_err_yellow[crossing]);

  ep->set_polarization_sys_err_blue(_pol_sys_err_blue[crossing]);

  ep->set_polarization_sys_err_yellow(_pol_sys_err_yellow[crossing]);

  if (Verbosity() >= 2)
    {
      cout << "saModLoadSpinInfo::event - ";
      ep->identify();
    }

  return EVENT_OK;
}

//! get spin state of a crossing. For bad crossings return saEventProperty::UNKOWN
saEventProperty::sa_spin_state
saModLoadSpinInfo::get_spin_state(int crossing)
{

  //! PHENIX bad bunches
  if (_spin_cont->GetBadBunchFlagPHENIX(crossing))
    {
      if (Verbosity() >= 2)
        cout << "saModLoadSpinInfo::get_spin_state - INFO - spin_state["
            << crossing << "] = UNKNOWN, rejected by spin database"
            <<" | PHENIX Crossing = "<< crossing << endl;

      return saEventProperty::UNKNOWN;
    }

  //! custom bad bunches
  int crossing_RHIC = (crossing + _spin_cont->GetCrossingShift())
              % saEventProperty::N_CROSSING;

  if (_use_spin_QA)
    if (_spin_QA.get_crossing_flag(_run_num, crossing_RHIC) <= 0)
      {
        if (Verbosity() >= 2)
          cout << "saModLoadSpinInfo::get_spin_state - INFO - spin_state_RHIC["
              << crossing_RHIC << "] = UNKNOWN for RHIC, rejected by custom spin QA"
              <<" | PHENIX Crossing = "<< crossing << endl;
        return saEventProperty::UNKNOWN;
      }

  //! Get spin states
  saEventProperty::sa_spin_state state = saEventProperty::UNKNOWN;

  const int sb = _spin_cont->GetSpinPatternBluePHENIX(crossing);
  const int sy = _spin_cont->GetSpinPatternYellowPHENIX(crossing);

  if (sb == 1)
    {
      if (sy == 1)
        state = saEventProperty::PP;
      else if (sy == -1)
        state = saEventProperty::PM;
    }
  else if (sb == -1)
    {
      if (sy == 1)
        state = saEventProperty::MP;
      else if (sy == -1)
        state = saEventProperty::MM;
    }

  if (Verbosity() >= 2)
    cout << "saModLoadSpinInfo::get_spin_state - INFO - spin_state["
        << crossing << "] = "<< state
        <<" | PHENIX Crossing = "<< crossing << endl;
  return state;
}

//! global termination
int
saModLoadSpinInfo::end(PHCompositeNode *topNode, sa_hist_mangager_ptr hm)
{
  return EVENT_OK;
}

//! load data files
void
saModLoadSpinInfo::SpinQA::load_crossing_flag(std::string data_file,
    int verbose)
{

  cout << "saModLoadSpinInfo::SpinQA::load_crossing_flag - INFO - loading "
      << data_file << endl;

  fstream f(data_file.c_str(), ios_base::in);

  if (!f.is_open())
    {

      cout
          << "saModLoadSpinInfo::SpinQA::load_crossing_flag - ERROR - cannot open "
          << data_file << endl;

    }

  int cnt = 0;
  string line;
  while (!f.eof())
    {
      getline(f, line);
      istringstream iss(line);

      int run = -1;
      iss >> run;

      if (run <= 0)
        {
          if (verbose)
            cout
                << "saModLoadSpinInfo::SpinQA::load_crossing_flag - WARNING - read invalid run "
                << run << endl;

          continue;

        }

      SpinQA::vec_crossing_flag flags;
      while (!iss.eof())
        {
          int flag = -9999;
          iss >> flag;
          if (verbose >= 2)
            cout << "     [" << flags.size() << "] = " << flag << " - "
                << (iss.rdstate() & ios_base::failbit) << " - " << iss.rdstate()
                << endl;

          if (iss.rdstate() & ios_base::failbit)
            break;

          flags.push_back(flag);
        }

      if (verbose)
        cout << "saModLoadSpinInfo::SpinQA::load_crossing_flag - INFO - read "
            << flags.size() << " flags for " << run << endl;

      if (flags.size() != (unsigned int) saEventProperty::N_CROSSING)
        cout
            << "saModLoadSpinInfo::SpinQA::load_crossing_flag - WARNING - read "
            << flags.size() << " flags for " << run
            << ", which is not equal to the RHIC crossing numbers" << endl;

      if (have_run_crossing_flag(run))
        {
          cout
              << "saModLoadSpinInfo::SpinQA::load_crossing_flag - WARNING - replacing flags for "
              << run << "." << endl;

          crossing_flag[run].clear();
        }

      crossing_flag[run] = flags;

      cnt++;
    }

  cout << "saModLoadSpinInfo::SpinQA::load_crossing_flag - INFO - loaded "
      << cnt << " records from " << data_file << endl;

}

//! have this run in data file?
bool
saModLoadSpinInfo::SpinQA::have_run_crossing_flag(int run)
{
  return crossing_flag.count(run) > 0;
}

//! get crossing status. 1 = good crossing. assuming RICH crossing numbers
int
saModLoadSpinInfo::SpinQA::get_crossing_flag(int run, int crossing)
{
  if (!have_run_crossing_flag(run))
    return 0;

  const vec_crossing_flag & flag = crossing_flag[run];

  if (crossing < 0 || crossing >= (int) flag.size())
    {
      cout
          << "saModLoadSpinInfo::SpinQA::get_crossing_flag - WARNING - invalid crossing "
          << crossing << " asked with flag array size of " << flag.size()
          << " for run " << run << endl;
      return 0;
    }

  return flag[crossing];

}
