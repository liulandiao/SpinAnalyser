// $Id: saHistManager.C,v 1.3 2014/04/15 05:19:59 jinhuang Exp $                                                                                             

/*!
 * \file saHistManager.C
 * \brief 
 * \author Jin Huang <jhuang@bnl.gov>
 * \version $Revision: 1.3 $
 * \date $Date: 2014/04/15 05:19:59 $
 */

#include <cassert>
#include <iostream>
#include <Fun4AllServer.h>
#include <Fun4AllHistoManager.h>

#include "saHistManager.h"

using namespace std;
using namespace boost;

saHistManager::saHistManager() :
    Fun4AllBase("saHistManager"), _default_lumi_name("")
{

}

saHistManager::~saHistManager()
{
  while (_hists.begin() != _hists.end())
    {
      if (Verbosity() >= 1)
        cout << "saHistManager::~saHistManager - delete hist "
            << _hists.begin()->first << endl;

      // actual histograms are owned by Fun4AllHistoManager. Not delete here
      _hists.erase(_hists.begin());
    }
}

Fun4AllHistoManager *
saHistManager::getFun4AllHistoManager()
{
  // get the histo manager for fun4all server

  Fun4AllServer* se = Fun4AllServer::instance();

  ostringstream histomanagername;
  histomanagername << se->Name() << "HISTOS";

  static bool once = true;

  if (once)
    {
      once = false;
      cout
          << "saHistManager::getFun4AllHistoManager - use Fun4All server histo mangager with name "
          << histomanagername.str() << endl;

    }

  Fun4AllHistoManager *hm = se->getHistoManager(histomanagername.str());
  if (!hm)
    {
      cout << "saHistManager::getFun4AllHistoManager - create "
          << histomanagername.str() << endl;

      hm = new Fun4AllHistoManager(histomanagername.str());
      se->registerHistoManager(hm);
    }

  assert(hm);

  return hm;
}

int
saHistManager::dumpHistos(const string &filename, const string &openmode)
{
  if (Verbosity() >= 1)
    {
      cout << "int Fun4AllHistoManager::dumpHistos() dumping histograms:"
          << endl;
      Print();
    }
  return getFun4AllHistoManager()->dumpHistos( //
      /*const string &*/filename, //
      /*const string &*/openmode);
}

//! get the histogram of default relative luminorsity
saHist *
saHistManager::getHisto_DefaultLumi() const
{
  if (_default_lumi_name.length() == 0)
    {
      cout << "saHistManager::getHisto_DefaultLumi - Error - "
          << "can not find the default luminosity histogram. Call registerHisto_DefaultLumi through a spin module first"
          << endl;

      return NULL;
    }

  return getHisto(_default_lumi_name);
}

//! register the histogram of default relative luminorsity
bool
saHistManager::registerHisto_DefaultLumi(saHist * h1d)
{

  assert(h1d);

  if (_default_lumi_name.length() > 0)
    {
      cout
          << "saHistManager::registerHisto_DefaultLumi - WARNING - overwrite the default luminorsity histogram from "
          << _default_lumi_name << " to " << h1d->GetName();
    }
  else
    {
      cout << "saHistManager::registerHisto_DefaultLumi - INFO - "
          << "Will use " << h1d->GetName()
          << " for relative luminorsity information" << endl;

      _default_lumi_name = h1d->GetName();
    }

  return registerHisto(h1d);

}

bool
saHistManager::registerHisto(saHist * h)
{
  const string hname(h->GetName());

  map<const string, sa_hist_ptr>::const_iterator histoiter = _hists.find(hname);
  if (histoiter != _hists.end())
    {
      cerr << "Histogram " << hname
          << " already registered, I won't overwrite it" << endl;
      cerr << "Use a different name and try again" << endl;
      return false;
    }

  h->Verbosity(Verbosity());
  _hists[hname] = sa_hist_ptr(h);
  getFun4AllHistoManager()->registerHisto(h);

  return true;
}

int
saHistManager::isHistoRegistered(const std::string &name) const
{

  map<const string, sa_hist_ptr>::const_iterator histoiter = _hists.find(name);
  if (histoiter != _hists.end())
    {
      return 1;
    }
  return 0;
}

saHist *
saHistManager::getHisto(const unsigned int ihisto) const
{
  map<const string, sa_hist_ptr>::const_iterator histoiter = _hists.begin();
  unsigned int size = _hists.size();
  if (Verbosity() > 3)
    {
      cout << "Map contains " << size << " Elements" << endl;
    }
  if (ihisto < size)
    {
      for (unsigned int i = 0; i < ihisto; i++)
        {
          ++histoiter;
        }
      return histoiter->second/*.get()*/;
    }
  else
    {
      cout << "saHistManager::getHisto: ERROR Invalid histogram number: "
          << ihisto << ", maximum number is " << size << endl;
    }
  return 0;
}

const char *
saHistManager::getHistoName(const unsigned int ihisto) const
{
  map<const string, sa_hist_ptr>::const_iterator histoiter = _hists.begin();
  unsigned int size = _hists.size();
  if (verbosity > 3)
    {
      cout << "Map contains " << size << " Elements" << endl;
    }
  if (ihisto < size)
    {
      for (unsigned int i = 0; i < ihisto; i++)
        {
          ++histoiter;
        }
      return histoiter->first.c_str();
    }
  else
    {
      cout << "saHistManager::getHisto: ERROR Invalid histogram number: "
          << ihisto << ", maximum number is " << size << endl;
    }
  return 0;
}

saHist *
saHistManager::getHisto(const string &hname) const
{
  map<const string, sa_hist_ptr>::const_iterator histoiter = _hists.find(hname);
  if (histoiter != _hists.end())
    {
      return histoiter->second/*.get()*/;
    }

  //second try, with auto suffix
  histoiter = _hists.find(hname + "_saHist");
  if (histoiter != _hists.end())
    {
      return histoiter->second/*.get()*/;
    }

  cout << "saHistManager::getHisto: ERROR Unknown Histogram " << hname
      << ", The following are implemented: " << endl;
  Print("ALL");
  return 0;
}

void
saHistManager::Print(const string &what) const
{
  cout << "saHistManager::Print - INFO - " << Name() << " contains "
      << nHistos() << " histograms. ";
  if (_default_lumi_name.length())
    cout << "The default lumi histogram is " << _default_lumi_name << endl;
  else
    cout << "The default lumi histogram was not set." << endl;

  if (what.find("ALL") != string::npos || what.find("HISTOS") != string::npos)
    {
      cout << "List of Histos in saHistManager " << Name() << ":" << endl;

      map<const string, sa_hist_ptr>::const_iterator hiter;
      for (hiter = _hists.begin(); hiter != _hists.end(); ++hiter)
        {
          cout << "\t" << hiter->first << " with " << hiter->second->nHistos()
              << " histograms" << endl;
          if (what.find("EXTEND") != string::npos)
            hiter->second->Print();
        }
      cout << endl;
    }
  return;
}

void
saHistManager::Reset()
{
  map<const string, sa_hist_ptr>::const_iterator hiter;
  for (hiter = _hists.begin(); hiter != _hists.end(); ++hiter)
    {

    }
  return;
}

unsigned int
saHistManager::nHistos() const
{
  return _hists.size();
}

//
//int
//saHistManager::dumpHistos(const string &filename, const string &openmode)
//{
//  int iret = 0;
//  if (!filename.empty())
//    {
//      outfilename = filename;
//    }
//  else
//    {
//      if (outfilename.empty())
//        {
//          recoConsts *rc = recoConsts::instance();
//          ostringstream filnam;
//          int runnumber = -1;
//          if (rc->FlagExist("RUNNUMBER"))
//            {
//              runnumber = rc->get_IntFlag("RUNNUMBER");
//            }
//          // this will set the filename to the name of the manager
//          // add the runnumber in the std 10 digit format and
//          // end it with a .root extension
//          filnam << Name() << "-"
//     << setfill('0') << setw(10)
//     << runnumber << ".root";
//          outfilename = filnam.str();
//        }
//    }
//  cout << "saHistManager::dumpHistos() Writing root file: " << outfilename << endl;
//
//  const int compress = 9;
//  ostringstream creator;
//  creator << "Created by " << Name();
//  TFile hfile(outfilename.c_str(), openmode.c_str(), creator.str().c_str(), compress);
//  if (!hfile.IsOpen())
//    {
//      cout << PHWHERE << " Could not open output file" << outfilename << endl;
//      return -1;
//    }
//
//  map<const string, TNamed *>::const_iterator hiter;
//  for (hiter = Histo.begin(); hiter != Histo.end(); ++hiter)
//    {
//      const std::string & hname = hiter->first;
//      const TNamed*       hptr  = hiter->second;
//      if ( verbosity > 0 )
//        {
//          std::cout << PHWHERE << " Saving histo "
//        << hname
//        << std::endl;
//        }
//
//      //  Decode the string to see if it wants a directory
//      string::size_type pos = hname.find_last_of('/');
//      string dirname, histoname;
//      if ( pos != string::npos ) // string::npos is the result if search unsuccessful
//        {
//          dirname = hname.substr(0, pos);
//          histoname = hname.substr(pos + 1);
//        }
//      else
//        {
//          dirname = "";
//          histoname = hname;
//        }
//
//      if (verbosity)
//        {
//          cout << " Histogram named " << hptr->GetName();
//    cout << " key " << hname;
//    if (dirname.size())
//      {
//        cout << " being saved to directory " << dirname;
//      }
//    cout << endl;
//        }
//
//      if (dirname.size())
//  {
//    TDirectoryHelper::mkdir(&hfile, dirname.c_str());
//    hfile.cd(dirname.c_str());
//  }
//
//      if (hptr)
//        {
//          int byteswritten = hptr->Write();
//    if (!byteswritten)
//      {
//        cout << PHWHERE << "Error saving histogram "
//       << hptr->GetName()
//       << endl;
//        iret = -2;
//      }
//        }
//      else
//        {
//          cout << PHWHERE << "dumpHistos : histogram "
//         << hname << " is a null pointer! Won't be saved."
//         << std::endl;
//        }
//    }
//  hfile.Close();
//  return iret;
//}
