// $Id: saModLoadSpinInfo.h,v 1.4 2013/07/25 08:47:26 jinhuang Exp $                                                                                             

/*!
 * \file saModLoadSpinInfo.h
 * \brief 
 * \author Jin Huang <jhuang@bnl.gov>
 * \version $Revision: 1.4 $
 * \date $Date: 2013/07/25 08:47:26 $
 */

#ifndef SAMODLOADSPININFO_H_
#define SAMODLOADSPININFO_H_

#ifndef __CINT__
#include <boost/smart_ptr.hpp>
#include <boost/array.hpp>
#endif
#include <vector>

#include "saEventProperty.h"
#include "saModuleBase.h"

class SpinDBContent;
class SpinDBOutput;
class EventHeader;
class SyncObject;
class RunHeader;
class TrigLvl1;

/*!
 * \brief saModLoadSpinInfo
 *
 * Load spin info from database following this example:
 *
 *
 gSystem->Load("libuspin");
 int run = 332252;
 SpinDBContent spin_cont;
 SpinDBOutput spin_out("phnxrc");
 spin_out.StoreDBContent(run,run);
 spin_out.GetDBContentStore(spin_cont,run);
 int bspin[120], yspin[120], bbc[120];
 for(int ibunch=0; ibunch<120; ibunch++){
 bspin[ibunch] = 0;
 yspin[ibunch] = 0;
 //      bspin[ibunch] = spin_cont.GetSpinPatternBlue(ibunch);
 //yspin[ibunch] = spin_cont.GetSpinPatternYellow(ibunch);
 //      cout << bspin[ibunch] << endl;
 //cout << yspin[ibunch] << endl;
 bbc[ibunch] = spin_cont.GetScalerBbcVertexCut(ibunch);
 cout << bbc[ibunch] << endl;

 }

 */
class saModLoadSpinInfo : public saModuleBase
{

  // --------------------------------------------------------------------------
  //!@name Reco interface
  //@{

public:
  saModLoadSpinInfo(const std::string &name = "SpinInfo");
  virtual
  ~saModLoadSpinInfo();

#ifndef __CINT__
  //! global initialization
  virtual int
  init(PHCompositeNode* topNode, sa_hist_mangager_ptr hm);
  //! Run initialization
  virtual int
  init_run(PHCompositeNode* topNode, sa_hist_mangager_ptr hm);
  //! event method
  virtual int
  event(PHCompositeNode* topNode, sa_hist_mangager_ptr hm);
  //! global termination
  virtual int
  end(PHCompositeNode* topNode, sa_hist_mangager_ptr hm);

#endif
  //@}

  // --------------------------------------------------------------------------
  //!@name Read spin data base
  //@{

public:

protected:

  //! get spin state of a crossing. For bad crossings return saEventProperty::UNKOWN
  virtual saEventProperty::sa_spin_state
  get_spin_state(int crossing);

  //! header pointers
  RunHeader* _run_header;
  EventHeader* _event_header;
  SyncObject* _sync;
  TrigLvl1* _trig_lv1;

  int _run_num;

#ifndef __CINT__
  typedef boost::shared_ptr<SpinDBContent> _spin_cont_ptr;
  _spin_cont_ptr _spin_cont;

  typedef boost::array<float, saEventProperty::N_CROSSING> crossing_arry;
  crossing_arry _pol_blue;
  crossing_arry _pol_stat_err_blue;
  crossing_arry _pol_sys_err_blue;
  crossing_arry _pol_yellow;
  crossing_arry _pol_stat_err_yellow;
  crossing_arry _pol_sys_err_yellow;
#endif

  //@}

  // --------------------------------------------------------------------------
  //!@name luminorsity histogram
  //@{

public:

  saHist*
  get_h_lumi() const
  {
    return _h_lumi;
  }

  void
  set_h_lumi(saHist* lumi)
  {
    _h_lumi = lumi;
  }

private:

  //! luminorsity in array of saHist::enu_lumi_id
  saHist * _h_lumi;

  //@}

//  //! luminorsity * polarization in 2D array of saHist::enu_pol_id vs saHist::enu_lumi_id
//  saHist * _h_lumi_pol;

  // --------------------------------------------------------------------------
  //!@name custom spin QA
  //@{
public:

  //! call this function to load beam position data
  //! in the order of file >> run >> x >> y >> ex >> ey;
  void
  load_crossing_flag(std::string data_file = "SpinQA.txt")
  {
    _spin_QA.load_crossing_flag(data_file, Verbosity());
    _use_spin_QA = true;
  }

  void
  set_use_spin_QA(bool b)
  {
    _use_spin_QA = b;
  }

  bool
  get_use_spin_QA() const
  {
    return _use_spin_QA;
  }

  //! Managing additional information on spin QA
  class SpinQA
  {
  public:
    //! load data files
    void
    load_crossing_flag(std::string data_file = "SpinQA.txt", int verbose = 0);

    //! have this run in data file?
    bool
    have_run_crossing_flag(int run);

    //! get crossing status. 1 = good crossing. assuming RICH crossing numbers
    int
    get_crossing_flag(int run, int crossing);

    //! position in x y
    typedef std::vector<int> vec_crossing_flag;
    typedef std::map<int, vec_crossing_flag> record_crossing_flag;

    //! map to run -> good crossing list, assuming RICH crossing numbers
    record_crossing_flag crossing_flag;

  };

  SpinQA _spin_QA;

private:

  //! whether to use user defined SpinQA _spin_QA;
  bool _use_spin_QA;

  //@}
};

#endif /* SAMODLOADSPININFO_H_ */
