// $Id: saHist.h,v 1.15 2015/08/17 20:09:10 jinhuang Exp $                                                                                             

/*!
 * \file saHist.h
 * \brief spin analysis histogram, self-managned collection of histograms for spin analysis
 * \author Jin Huang <jhuang@bnl.gov>
 * \version $Revision: 1.15 $
 * \date $Date: 2015/08/17 20:09:10 $
 */

#ifndef SAHIST_H_
#define SAHIST_H_

#include <string>
#include <iomanip>
#include <sstream>
#include <cmath>
#include <map>

#include <TNamed.h>
#include <TRef.h>

#include "saEventProperty.h"

class TH1;

/*!
 * \brief spin analysis histogram, self-managned collection of histograms for spin analysis
 */
class saHist : public TNamed
{

  // --------------------------------------------------------------------------
  //!@name Flags
  //@{

public:


  enum en_flag
  {

    //! this is a sub-histogram of another saHist object
    SUB_HIST = 1 << 0,

    //! depends on spin and for asymmetry analysis?
    SPIN_DEPENDENT = 1 << 1,

    //! filled event by event. e.g. J/Psi count
    EVENT_PROPERTY = 1 << 2,

    //! property for full run and will fill sub histograms of runs
    RUN_PROPERTY = 1 << 3,

    //! property for full fill and will fill sub histograms of fills
    FILL_PROPERTY = 1 << 4,

    //! also produce histograms of product of polarizations
    POL_PRODUCTS = 1 << 5,

    //! produce histograms to check RHIC crossing alignment
    CROSSING_CHECKS = 1 << 6,

    //! if CROSSING_CHECKS is set, produce histograms to check RHIC crossing alignment for each bins
    CROSSING_CHECKS_EACH_BIN = 1 << 7,

    //! This is relative lumi histograms. DO not merge for same run
    LUMI_HIST = 1 << 8,

    DUMMY
  };

  typedef unsigned int sa_flags_t;

  //! default flag for asymmetry analysis
  static const sa_flags_t DEFAULT_FLAG = //
      SPIN_DEPENDENT | //
          EVENT_PROPERTY;
  //! default flag for relative luminorsity
  static const sa_flags_t DEFAULT_FLAG_LUMI = //
      SPIN_DEPENDENT | //
          RUN_PROPERTY | //
          FILL_PROPERTY | //
          POL_PRODUCTS | //
          CROSSING_CHECKS | //
          CROSSING_CHECKS_EACH_BIN | //
          LUMI_HIST;

  virtual void
  set_flag(const en_flag& flag, const bool& value)
  {
    if (value)
      {
        _flags |= flag;
      }
    else
      {
        _flags &= (~flag);
      }
  }

  virtual bool
  get_flag(const en_flag& flag) const
  {
    return (_flags & flag) > 0;
  }

  sa_flags_t
  get_flags() const
  {
    return _flags;
  }

  void
  set_flags(sa_flags_t flags)
  {
    _flags = flags;
  }

  //! Sets the verbosity of this module (0 by default=quiet).
  virtual void
  Verbosity(const int ival)
  {
    verbosity = ival;
  }

  //! set verbosity to sub saHists too
  void Verbosity(int verb, bool recursive);

  //! Gets the verbosity of this module.
  virtual int
  Verbosity() const
  {
    return verbosity;
  }

protected:
  //! The verbosity level. 0 means not verbose at all.
  int verbosity;
  //! flags for this object
  sa_flags_t _flags;

  //@}

  // --------------------------------------------------------------------------
  //!@name Interface to ROOT and C++
  //@{
public:
  saHist(TH1* h, const sa_flags_t flag = DEFAULT_FLAG, int verbosity = 0);
  saHist();

  virtual
  ~saHist();

  //! Print content
  //@param option "EXTEND" means print details of sub spin analysis histograms too
  virtual void
  Print(Option_t* option = "") const;

  virtual Int_t
  Write(const char* name = 0, Int_t option = 0, Int_t bufsize = 0) const;

  //! load all histograms from a root files
  //! @param[in] verbosity negative verbosity means use the verbosity of this saHist
  //! @return true if all operation is successful
  virtual
  bool
  AutoLoad(int verbosity = -1) ;

protected:

  //@}
  // --------------------------------------------------------------------------
  //!@name Interface to saModuleBase
  //@{

public:

  //! base name of this object and all sub histograms
  std::string
  get_name() const
  {
    return _name;
  }

  //! Initalize to a histogram
  virtual
  void
  Init(TH1* h);

  //! fill the histogram and automatically retrive saEventProperty from the TOP node
  virtual Int_t
  Fill(double v1, double v2 = 1, double v3 = 1, double v4 = 1,
      double extra_weight = 1);

  //! fill the histogram
  virtual Int_t
  Fill(const saEventProperty* ep, double v1, double v2 = 1, double v3 = 1,
      double v4 = 1, double extra_weight = 1);

  //! directly fill the internal histogram with a name
  //! @return: the weight for this filling
  virtual
  double
  Fill(const std::string& suffix, double v1, double v2 = 1, double v3 = 1,
      double v4 = 1, double extra_weight = 1) const;

  //! full name from a suffix
  virtual std::string
  make_hist_name(const std::string & suffix) const
  {
    return _name + "_" + suffix;
  }

protected:

  void
  set_name(std::string name)
  {
    _name = name;
  }

  //! base name of this object
  std::string _name;

  //! internal filling function
  //! @return: the weight for this filling
  static double
  Fill(TH1* h, double v1, double v2, double v3, double v4, double extra_weight =
      1);

  //@}

  // --------------------------------------------------------------------------
  //!@name Collection of Histograms
  //@{
public:

  //! whether a histogram (TH1) with a suffix is registered
  int
  isHistoRegistered(const std::string& suffix) const;

  //! get the histogram (TH1) with a suffix
  TH1*
  getHisto(const std::string& suffix) const;

  //! get the histogram (TH1) with index of ihisto
  TH1*
  getHisto(const unsigned int ihisto) const;

  //! get the histogram (TH1) corresponding to a spin state
  TH1*
  getHisto(const saEventProperty::sa_spin_state state) const;

  //! get the name of histogram (TH1) with index of ihisto
  const char*
  getHistoName(const unsigned int ihisto) const;

  //! number of histograms (TH1)
  unsigned int
  nHistos() const;

  //! remove a histogram from the list
  bool
  removeHisto(const std::string& suffix);

protected:

  //! get a clean clone of the base object
  TObject*
  get_clone_hist(const std::string& full_name) const;

  //! clone without dirty association mess with TRef and TDirectory
  static TH1 * CleanClone(TH1 * obj_src, const std::string & name);

  //! insert a clone histogram with suffix into internal hist lists
  // @return true if a new insert was successful, false if the histograme is already there
  bool
  insert_clone_hist(const std::string& suffix);

  //! template for the histogram set in constructor
  typedef std::map<const std::string, TRef> sa_vec_hist;

  //! collection of histograms (TH1) which is copies of the base histogram (getHisto("BASE"))
  sa_vec_hist _hists;

  //@}

  // --------------------------------------------------------------------------
  //!@name Collection of Spin Analysis Histograms
  //@{
public:

  //! whether a spin analysis histogram (saHist) with a suffix is registered
  int
  issaHistRegistered(const std::string& suffix) const;

  //! get the spin analysis histogram (saHist) with a suffix
  saHist *
  getsaHist(const std::string& suffix, bool add = false);

  //! get the spin analysis histogram (saHist) with index of ihisto
  saHist *
  getsaHist(const unsigned int ihisto) const;

  //! get the spin analysis histogram (saHist) for a run
  saHist *
  getsaHist_run(const unsigned int run, bool add = false)
  {
    return getsaHist((get_suffix_run(run)), add);
  }

  //! get the spin analysis histogram (saHist) for a fill
  saHist *
  getsaHist_fill(const unsigned int fill, bool add = false)
  {
    return getsaHist((get_suffix_fill(fill)), add);
  }

  //! get the name of spin analysis histogram (saHist) with index of ihisto
  const char*
  getsaHistName(const unsigned int ihisto) const;

  //! number of spin analysis histograms (saHist)
  unsigned int
  nsaHists() const;

  //! suffix for saHist for one run
  static std::string
  get_suffix_run(const unsigned int run)
  {
    std::ostringstream str;

    str << "RUN_";
    str << std::setfill('0') << std::setw(10);
    str << run;

    return str.str();
  }

  //! suffix for saHist for one fill
  static std::string
  get_suffix_fill(const unsigned int fill)
  {
    std::ostringstream str;

    str << "FILL_";
    str << std::setfill('0') << std::setw(6);
    str << fill;

    return str.str();
  }

  typedef std::vector<int> n_list;

  n_list
  get_run_list(void);

  n_list
  get_fill_list(void);

  //-----------------------------------------------------------------------------------------------
  //! remove a sahist from the list
  bool
  removeFillsaHist(const std::string& suffix);
  
  bool
  removeFillsaHist(const unsigned int fill)
  {
    return removeFillsaHist(get_suffix_fill(fill));
  }

  bool
  removeRunsaHist(const std::string& suffix);
  
  bool
  removeRunsaHist(const unsigned int run)
  {
    return removeRunsaHist(get_suffix_run(run));
  }

protected:

  //! insert a clone spin analysis histogram with suffix into internal hist lists
  // @return true if a new insert was successful, false if the histograme is already there
  bool
  insert_clone_sahist(const std::string& suffix);

  //! collection of spin analysis histogram (saHist)
  sa_vec_hist _sahists;

  //@}

  // --------------------------------------------------------------------------
  //!@name Calculate Asymmetries
  //@{

public:

  enum enu_lumi_id
  {
    BbcVertexCut = 1,
    BbcNoCutPHENIX,
    ZdcWidePHENIX,
    ScalerZdcNarrow,
    ScalerFVTX, //

    //! total count of IDs
    LUMI_COUNT
  };

  enum enu_pol_id
  { //
    //! blue polarization
    Pol_Blue, //

    //! yellow polarization
    Pol_Yellow, //

    //! blue polarization *  yellow polarization
    Pol2_Blue_Yellow, //

    //! total count of IDs
    Pol_COUNT
  };

  //! Calculate asymmetry according to luminorsity histogram h_lumi
  //@return true if successful
  virtual
  bool
  CalcAsymmetry(saHist* h_lumi, enu_lumi_id lumi_type);

  //! Calculate final asymmetry by mini-Chi2 fit on run-by-run or fill-by-fill asymmetry results from the sub saHits
  //! \param[in] run_or_fill true : use the run list, false : use the fill list
  //! \return return a new sub saHist, storing the new asymmetry fit results
  virtual
  saHist *
  CalcAsymmetry_Chi2Fit(bool run_or_fill = true);

  virtual
  bool
  CalcAsymmetry_Chi2Fit_Work(bool run_or_fill, saHist *sahAsym, TString hist_type);

  //! remove previously saved asymmetry histograms
  void
  RemoveAsymmetryHisto();

protected:

  //! Calculate asymmetry according to luminorsity histogram h_lumi, DO NOT iterate through sub-histograms
  //@return true if successful
  virtual
  bool
  CalcAsymmetry_Work(saHist* h_lumi, enu_lumi_id lumi_type);

  //! rename of std::pow()
  static inline double
  Power(const double a, const int n)
  {
    return std::pow(a, n);
  }

  //@}

  // --------------------------------------------------------------------------
  //!@name Merge utilities
  //@{

public:
  //! use this saHist as template to make an unfilled new saHist
  //! @param[in] new_name empty means use the name of the source saHist
  virtual saHist *
  MakeTemplate(std::string new_name = "");

  //! Merge count in hist to this object
  virtual
  void
  MergeHistCount(saHist * hist);

  //! Add a saHist as sub histogram, merge its counts to the master histograms
  //!  and possess the ownership to this object
  virtual
  void
  MergeSubHist(const std::string  & suffix, saHist * sub_hist);

  //! adopt sub saHist of a source master saHist to this object.
  //! Those sub saHist will be removed from the source
  //! @param[in] run_or_fill true : adopt the run list, false : adopt the fill list
  virtual
  void
  AdoptSubHist(saHist * master_hist, bool run_or_fill = true )
  {
    n_list l;
    AdoptSubHist(master_hist, run_or_fill, false, l);
  }

  //! adopt sub saHist of a source master saHist to this object.
  //! Support use of a good run list. Those sub saHist will be removed from the source.
  //! @param[in] run_or_fill true : adopt the run list, false : adopt the fill list
  virtual
  void
  AdoptSubHist(saHist * master_hist, bool run_or_fill , bool use_good_list,
      const n_list & good_list);

  //---------------------------------------------------------------------------------------------------------------

  //! Add a saHist as sub histogram, merge its counts to the master histograms
  //!  and possess the ownership to this object
  virtual
  void
  MergeRunSubHist(const std::string  & suffix, saHist * sub_hist);


  //! adopt run sub saHist of a source master saHist to this object.
  //! Those run sub saHist will be removed from the source
  virtual
  void
  AdoptRunSubHist(saHist * master_hist)
  {
    n_list l;
    AdoptRunSubHist(master_hist, false, l);
  }

  //! adopt run sub saHist of a source master saHist to this object.
  //! Support use of a good run list. Those sub saHist will be removed from the source.
  virtual
  void
  AdoptRunSubHist(saHist * master_hist, bool use_good_list,
      const n_list & good_list);

  //! Add a saHist as sub histogram, merge its counts to the master histograms
  //!  and possess the ownership to this object
  virtual
  void
  MergeFillSubHist(const std::string  & suffix, saHist * sub_hist);

  //! adopt run sub saHist of a source master saHist to this object.
  //! Those run sub saHist will be removed from the source
  virtual
  void
  AdoptFillSubHist(saHist * master_hist)
  {
    n_list l;
    AdoptFillSubHist(master_hist, false, l);
  }

  //! adopt run sub saHist of a source master saHist to this object.
  //! Support use of a good run list. Those sub saHist will be removed from the source.
  virtual
  void
  AdoptFillSubHist(saHist * master_hist, bool use_good_list,
      const n_list & good_list);


  //@}
ClassDef(saHist, 2) // spin analyzer's histogram manager
};

#endif /* SAHIST_H_ */
