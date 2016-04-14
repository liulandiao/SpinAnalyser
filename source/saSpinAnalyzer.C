// $Id: saSpinAnalyzer.C,v 1.5 2014/04/15 05:21:21 jinhuang Exp $                                                                                             

/*!
 * \file saSpinAnalyzer.C
 * \brief 
 * \author Jin Huang <jhuang@bnl.gov>
 * \version $Revision: 1.5 $
 * \date $Date: 2014/04/15 05:21:21 $
 */

#include <iostream>
#include <cassert>
#include <stdexcept>
#include <vector>
#include <utility>

#include <Fun4AllReturnCodes.h>
#include <PHNodeIterator.h>
#include <TMutNode.h>

#include "saModuleBase.h"
#include "saEventProperty_v1.h"
#include "saModLoadSpinInfo.h"

#include "saSpinAnalyzer.h"

using namespace std;
using namespace boost;

saSpinAnalyzer::saSpinAnalyzer(const std::string &name,
    saModuleBase * mod_loadspin) :
    SubsysReco(name), _initialized(false), //
    _hist_manager(make_shared<saHistManager>()), //
    _auto_save_hist(""), _mspin(NULL)
{
  // load spin info from database

  if (!mod_loadspin)
    {
      cout
          << "saSpinAnalyzer::constructor - Info - use the default spin info module (saModLoadSpinInfo)."
          << endl;
      // use the default spin module if not specified
      mod_loadspin = new saModLoadSpinInfo();

    }

  set_mspin(mod_loadspin);
}

saSpinAnalyzer::~saSpinAnalyzer()
{
  if (Verbosity())
    cout << "saSpinAnalyzer::~saSpinAnalyzer - Destructor. "
        << "Auto deleting modules and histograms through boost::shared_ptr"
        << endl;
}

//! Sets the verbosity of this module (0 by default=quiet).
void
saSpinAnalyzer::Verbosity(const int ival)
{
  SubsysReco::Verbosity(ival);

  if (_hist_manager.get())
    _hist_manager->Verbosity(Verbosity());

  if (_mspin)
    _mspin->Verbosity(Verbosity());

}

//! global initialization
int
saSpinAnalyzer::Init(PHCompositeNode *topNode)
{

  //internal call
    {
      int ret = init(topNode);
      if (ret != EVENT_OK)
        return ret;
    }

  sa_vec_module::iterator hiter;
  for (hiter = _modules.begin(); hiter != _modules.end(); ++hiter)
    {
      const string & name = hiter->first;

      if (Verbosity() >= 1)
        {
          cout << "saSpinAnalyzer::Init start - module " << name << endl;
        }

      try
        {
          int ret = hiter->second.get()->init(topNode, _hist_manager);

          if (ret != EVENT_OK)
            return ret;
        }
      catch (const std::exception & e)
        {

          cout << "saSpinAnalyzer::Init - module " << name << " error : "
              << e.what() << endl;

          return ABORTRUN;

        }

      if (Verbosity() >= 1)
        {
          cout << "saSpinAnalyzer::Init finish - module " << name << endl;
        }
    }

  _initialized = true;

  return EVENT_OK;
}

//! Run initialization
int
saSpinAnalyzer::InitRun(PHCompositeNode *topNode)
{

  //internal call
    {
      int ret = init_run(topNode);
      if (ret != EVENT_OK)
        return ret;
    }

  sa_vec_module::iterator hiter;
  for (hiter = _modules.begin(); hiter != _modules.end(); ++hiter)
    {
      const string & name = hiter->first;

      if (Verbosity() >= 1)
        {
          cout << "saSpinAnalyzer::InitRun start - module " << name << endl;
        }

      try
        {

          int ret = hiter->second.get()->init_run(topNode, _hist_manager);

          if (ret != EVENT_OK)
            return ret;
        }
      catch (const std::exception & e)
        {

          cout << "saSpinAnalyzer::InitRun - module " << name << " error : "
              << e.what() << endl;

          return ABORTRUN;

        }

      if (Verbosity() >= 1)
        {
          cout << "saSpinAnalyzer::InitRun finish - module " << name << endl;
        }
    }

  return EVENT_OK;
}

//! event method
int
saSpinAnalyzer::process_event(PHCompositeNode *topNode)
{

  //internal call
    {
      int ret = event(topNode);
      if (ret != EVENT_OK)
        return ret;
    }

  int ret_global = EVENT_OK;

  sa_vec_module::iterator hiter;
  for (hiter = _modules.begin(); hiter != _modules.end(); ++hiter)
    {
      const string & name = hiter->first;

      if (Verbosity() >= 3)
        {
          cout << "saSpinAnalyzer::process_event start - module " << name
              << endl;
        }

      try
        {

          int ret = hiter->second.get()->event(topNode, _hist_manager);

          switch (ret)
            {
          case EVENT_OK:
            break;

          case DISCARDEVENT:
            ret_global = DISCARDEVENT;
            break;

          default:
            return ret;
            break;

            }
        }
      catch (const std::exception & e)
        {

          cout << "saSpinAnalyzer::InitRun - module " << name << " error : "
              << e.what() << endl;

          return ABORTRUN;

        }

      if (Verbosity() >= 3)
        {
          cout << "saSpinAnalyzer::process_event finish - module " << name
              << endl;
        }
    }

  return ret_global;
}

//! global termination
int
saSpinAnalyzer::End(PHCompositeNode *topNode)
{
  //internal call
    {
      int ret = end(topNode);
      if (ret != EVENT_OK)
        return ret;
    }

  sa_vec_module::iterator hiter;
  for (hiter = _modules.begin(); hiter != _modules.end(); ++hiter)
    {
      const string & name = hiter->first;

      if (Verbosity() >= 1)
        {
          cout << "saSpinAnalyzer::End start - module " << name << endl;
        }

      try
        {

          int ret = hiter->second.get()->end(topNode, _hist_manager);

          if (ret != EVENT_OK)
            return ret;
        }
      catch (const std::exception & e)
        {

          cout << "saSpinAnalyzer::End - module " << name << " error : "
              << e.what() << endl;

          return ABORTRUN;

        }

      if (Verbosity() >= 1)
        {
          cout << "saSpinAnalyzer::End finish - module " << name << endl;
        }
    }

  if (_auto_save_hist.length())
    {
      cout << "saSpinAnalyzer::End - INFO - "
          << "automatically save histograms to " << _auto_save_hist << endl;

      _hist_manager->dumpHistos(_auto_save_hist);
    }

  return EVENT_OK;
}

//! internal global initialization
int
saSpinAnalyzer::init(PHCompositeNode *topNode)
{

  return EVENT_OK;
}

//! internal Run initialization
int
saSpinAnalyzer::init_run(PHCompositeNode *topNode)
{

  return EVENT_OK;
}

//! internal event method
int
saSpinAnalyzer::event(PHCompositeNode *topNode)
{

  return EVENT_OK;
}

//! internal global termination
int
saSpinAnalyzer::end(PHCompositeNode *topNode)
{
  cout << "saSpinAnalyzer::end - ending processing. Current histograms:";
  _hist_manager->Print("");

  return EVENT_OK;
}

int
saSpinAnalyzer::dumpHistos(const string &filename, const string &openmode)
{
  if (Verbosity() >= 1)
    cout << "int saSpinAnalyzer::dumpHistos() dumping histograms" << endl;

  return _hist_manager.get()->dumpHistos( //
      /*const string &*/filename, //
      /*const string &*/openmode);
}

int
saSpinAnalyzer::RegisterModule(saModuleBase * m, bool is_first_module)
{
  if (!m)
    {
      cout
          << "saSpinAnalyzer::RegisterModule - Error - Registering an invalid module!";

      exit(0);
    }

  const string name(m->Name());

  if (_initialized)
    {
      cout
          << "saSpinAnalyzer::RegisterModule - Error - already initialized before register of "
          << name
          << ". Please do saSpinAnalyzer -> RegisterModule before Fun4AllServer -> registerSubsystem. Force exit...";

      exit(0);
    }

  sa_vec_module::const_iterator hiter;
  for (hiter = _modules.begin(); hiter != _modules.end(); ++hiter)
    {
      if (hiter->first == name)
        break;
    }
  if (hiter != _modules.end())
    {
      cerr << "Module " << name << " already registered, I won't overwrite it"
          << endl;
      cerr << "Use a different name and try again" << endl;
      return 0;
    }

  if (Verbosity() >= 1)
    {

      cout << "saSpinAnalyzer::RegisterModule - adding new module " << name
          << endl;

    }

  if (is_first_module)
    {
      _modules.insert(_modules.begin(),
          make_pair<const std::string, sa_module_ptr>(name, sa_module_ptr(m)));
    }
  else
    {
      _modules.push_back(
          make_pair<const std::string, sa_module_ptr>(name, sa_module_ptr(m)));
    }
  return 1;
}

int
saSpinAnalyzer::isModuleRegistered(const std::string &name) const
{

  sa_vec_module::const_iterator hiter;
  for (hiter = _modules.begin(); hiter != _modules.end(); ++hiter)
    {
      if (hiter->first == name)
        break;
    }
  if (hiter != _modules.end())
    {
      return 1;
    }
  return 0;
}

int
saSpinAnalyzer::RemoveModule(const std::string &name)
{
  sa_vec_module::iterator hiter;
  for (hiter = _modules.begin(); hiter != _modules.end(); ++hiter)
    {
      if (hiter->first == name)
        break;
    }
  if (hiter != _modules.end())
    {
      _modules.erase(hiter);
      return 1;
    }
  return 0;
}

saModuleBase *
saSpinAnalyzer::getModule(const unsigned int iModule) const
{

  if (iModule < nModules())
    {
      return _modules[iModule].second.get();
    }
  else
    {
      cout << "saSpinAnalyzer::getModule: ERROR Invalid Modulegram number: "
          << iModule << ", maximum number is " << nModules() << endl;
    }
  return 0;
}

const char *
saSpinAnalyzer::getModuleName(const unsigned int iModule) const
{
  sa_vec_module::const_iterator Moduleiter = _modules.begin();
  unsigned int size = _modules.size();
  if (verbosity > 3)
    {
      cout << "Map contains " << size << " Elements" << endl;
    }
  if (iModule < size)
    {
      for (unsigned int i = 0; i < iModule; i++)
        {
          ++Moduleiter;
        }
      return Moduleiter->first.c_str();
    }
  else
    {
      cout << "saSpinAnalyzer::getModuleName: ERROR Invalid Modulegram number: "
          << iModule << ", maximum number is " << size << endl;
    }
  return 0;
}

saModuleBase *
saSpinAnalyzer::getModule(const string &name) const
{
  sa_vec_module::const_iterator hiter;
  for (hiter = _modules.begin(); hiter != _modules.end(); ++hiter)
    {
      if (hiter->first == name)
        break;
    }
  if (hiter != _modules.end())
    {
      return hiter->second.get();
    }
  cout << "saSpinAnalyzer::getModule: ERROR Unknown Modulegram " << name
      << ", The following are implemented: " << endl;
  Print("ALL");
  return 0;
}

unsigned int
saSpinAnalyzer::nModules() const
{
  return _modules.size();
}

saHistManager *
saSpinAnalyzer::get_hist_mangager()
{
  return _hist_manager.get();
}

//! register an spin module
void
saSpinAnalyzer::set_mspin(saModuleBase* mspin)
{
  if (!mspin)
    {
      cout
          << "saSpinAnalyzer::set_mspin - Error - Registering an invalid module!";

      exit(0);
    }
  std::cout << "saSpinAnalyzer::set_mspin - use spin information module "
      << mspin->Name() << std::endl;

  if (_mspin)
    {
      std::cout
          << "saSpinAnalyzer::set_mspin - removing the old spin information module first:"
          << _mspin->Name() << std::endl;

      RemoveModule(_mspin->Name());
    }

  _mspin = mspin;

  // always insert to the first module
  RegisterModule(_mspin, true);
}
