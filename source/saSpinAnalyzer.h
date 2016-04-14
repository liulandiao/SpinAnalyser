// $Id: saSpinAnalyzer.h,v 1.4 2013/07/24 07:37:04 jinhuang Exp $                                                                                             

/*!
 * \file saSpinAnalyzer.h
 * \brief 
 * \author Jin Huang <jhuang@bnl.gov>
 * \version $Revision: 1.4 $
 * \date $Date: 2013/07/24 07:37:04 $
 */

#ifndef SASPINANALYZER_H_
#define SASPINANALYZER_H_
#include <SubsysReco.h>
#include <string>
#include <utility>
#include <map>

#ifndef __CINT__
#include <boost/smart_ptr.hpp>
#endif

class saModuleBase;

/*!
 * \brief saSpinAnalyzer
 */
class saSpinAnalyzer : public SubsysReco
{

  // --------------------------------------------------------------------------
  //!@name Fun4All interface
  //@{

public:
  //! constructor
  // @param[in] name          name of this SubsysReco
  // @param[in] mod_loadspin  pointer to the default spin info module, which load the pol. and lumi, etc.
  saSpinAnalyzer(const std::string &name = "SpinAnalyzer",
      saModuleBase * mod_loadspin = NULL);

  //! destructor
  virtual
  ~saSpinAnalyzer();

  //! global initialization
  virtual int
  Init(PHCompositeNode *topNode);

  //! Run initialization
  virtual int
  InitRun(PHCompositeNode *topNode);

  //! event method
  virtual int
  process_event(PHCompositeNode *topNode);

  //! global termination
  virtual int
  End(PHCompositeNode *topNode);

  int
  dumpHistos(const std::string &filename = "", const std::string &openmode =
      "RECREATE");

  //! Sets the verbosity of this module (0 by default=quiet).
  virtual void
  Verbosity(const int ival);

  /// Gets the verbosity of this module.
  virtual int
  Verbosity() const
  {
    return SubsysReco::Verbosity();
  }

protected:
  //! internal global initialization
  virtual int
  init(PHCompositeNode *topNode);

  //! internal Run initialization
  virtual int
  init_run(PHCompositeNode *topNode);

  //! internal event method
  virtual int
  event(PHCompositeNode *topNode);

  //! internal global termination
  virtual int
  end(PHCompositeNode *topNode);

  bool _initialized;

  //@}

  // --------------------------------------------------------------------------
  //!@name Module IO
  //@{
public:

  //! add module
  int
  RegisterModule(saModuleBase * m, bool is_first_module = false);
  int
  RemoveModule(const std::string &name);
  int
  isModuleRegistered(const std::string &name) const;
  saModuleBase *
  getModule(const std::string &hname) const;
  saModuleBase *
  getModule(const unsigned int iModule) const;
  const char *
  getModuleName(const unsigned int iModule) const;
  unsigned int
  nModules() const;

  saHistManager *
  get_hist_mangager();

  std::string
  get_auto_save_hist() const
  {
    return _auto_save_hist;
  }

  void
  set_auto_save_hist(std::string autoSave)
  {
    _auto_save_hist = autoSave;
  }

  //! register an spin module
  void
  set_mspin(saModuleBase* mspin);

  saModuleBase*
  get_mspin() const
  {
    return _mspin;
  }

protected:

#ifndef __CINT__

  typedef boost::shared_ptr<saModuleBase> sa_module_ptr;
  //  typedef std::map<const std::string, sa_module_ptr> sa_vec_module;
  typedef std::pair<std::string, sa_module_ptr> module_rec;
  typedef std::vector<module_rec> sa_vec_module;

  //! list of modules
  sa_vec_module _modules;

  //! list of histograms
  sa_hist_mangager_ptr _hist_manager;
#endif

  //! file name for auto save of the _hist_manager
  std::string _auto_save_hist;

  //! spin info module
  saModuleBase * _mspin;

  //@}

};

#endif /* SASPINANALYZER_H_ */
