// $Id: saHist.C,v 1.22 2015/08/17 20:09:10 jinhuang Exp $                                                                                             

/*!
 * \file saHist.cpp
 * \brief 
 * \author Jin Huang <jhuang@bnl.gov>
 * \version $Revision: 1.22 $
 * \date $Date: 2015/08/17 20:09:10 $
 */

#include <iostream>
#include <cassert>
#include <stdexcept>
#include <algorithm>    // std::find
#include <vector>       // std::vector
#include <TH1.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TH2.h>
#include <TH3.h>
#include <TClass.h>
#include <TProcessID.h>
#include <TObjArray.h>
#include <TDirectory.h>
#include <TGraphErrors.h>
#include <TFitResultPtr.h>
#include <TFitResult.h>
//#include <TFile.h>

#include <Fun4AllServer.h>
#include <getClass.h>
#include <TMutNode.h>

#include "saEventProperty.h"

#include "saHist.h"

using namespace std;

ClassImp(saHist);

saHist::saHist() :
    verbosity(0), _flags(DEFAULT_FLAG), _name("NONAME")
{
}

saHist::saHist(TH1 * h, const sa_flags_t flag, int verbosity) :
    verbosity(verbosity), _flags(flag), _name("NONAME")
{
  Init(h);
}

saHist::~saHist()
{

  if (Verbosity() >= 1)
    cout << "saHist::" << _name
        << "::destructor - clean sub spin ana histograms" << endl;

  while (_sahists.begin() != _sahists.end())
    {
      TRef & ref = _sahists.begin()->second;

      if (ref.IsValid() && ref.GetObject())
        delete ref.GetObject();

      _sahists.erase(_sahists.begin());
    }

  if (Verbosity() >= 1)
    cout << "saHist::" << _name << "::destructor - clean histograms" << endl;

  while (_hists.begin() != _hists.end())
    {
      TRef & ref = _hists.begin()->second;

      if (ref.IsValid() && ref.GetObject())
        delete ref.GetObject();

      _hists.erase(_hists.begin());
    }
}

void
saHist::Verbosity(int verb, bool recursive)
{
  Verbosity(verb);

  if (recursive)
    for (sa_vec_hist::iterator hiter = _sahists.begin();
        hiter != _sahists.end(); ++hiter)
      {
        TRef & ref = (*hiter).second;
        if (ref.IsValid() && dynamic_cast<saHist *>(ref.GetObject()))
          {
            static_cast<saHist *>(ref.GetObject())->Verbosity(verb);
          }
      }

}

//! insert a clone histogram with suffix into internal hist lists
bool
saHist::insert_clone_hist(const std::string & suffix)
{

  if (!isHistoRegistered(suffix))
    {

      _hists[suffix] = get_clone_hist(make_hist_name(suffix));

      return true;
    }
  else
    return false;

}

//! insert a clone spin analysis histogram with suffix into internal lists
bool
saHist::insert_clone_sahist(const std::string & suffix)
{

  if (!issaHistRegistered(suffix))
    {

      TH1 * new_base =
          static_cast<TH1 *>(get_clone_hist(make_hist_name(suffix)));

      sa_flags_t new_flags = get_flags();

      new_flags &= ~RUN_PROPERTY;
      new_flags &= ~FILL_PROPERTY;

      new_flags |= SUB_HIST;

      saHist * new_hist = new saHist(new_base, new_flags, Verbosity());

      _sahists[suffix] = new_hist;

      return true;
    }
  else
    return false;

}

TObject *
saHist::get_clone_hist(const std::string & name) const
{

  sa_vec_hist::const_iterator histoiter = _hists.find("BASE");
  if (histoiter == _hists.end())
    {
      cout
          << "saHist::get_clone_hist - missing the base histogram for clone to "
          << name << endl;

      throw runtime_error(
          "saHist::get_clone_hist - missing the base histogram");
    }

  const TRef & ref = histoiter->second;
  TH1 * obj_src = dynamic_cast<TH1 *>(ref.GetObject());
  if (!obj_src)
    {
      cout << "saHist::" << _name
          << "::get_clone_hist - invalid base histogram for clone to " << name
          << endl;

      throw runtime_error("saHist::get_clone_hist - invalid base histogram");
    }

  if (Verbosity() >= 3)
    {
      // debug Clone

      cout << "saHist::" << _name << "::get_clone_hist - Directory = "
          << gDirectory->GetName() << endl;
      cout << "saHist::" << _name << "::get_clone_hist src TRef = "
          << ref.GetUniqueID() << endl;
      cout << "saHist::" << _name << "::get_clone_hist obj_src = "
          << obj_src->GetName() << ", id = " << obj_src->GetUniqueID()
          << ", sum = " << static_cast<const TH1 *>(obj_src)->GetSumOfWeights()
          << endl;

      cout << "saHist::" << _name << "::get_clone_hist TProcessID = ";
      TProcessID::GetSessionProcessID()->Print();
      TObjArray * arr = TProcessID::GetSessionProcessID()->GetObjects();
      assert(arr);
      //      arr->Print();
      for (int i = 1; i < arr->GetSize(); i++)
        {
          TObject * obj = arr->UncheckedAt(i);
          if (obj)
            {

              cout << "\t arr @ " << i << " = " << obj->GetName() << ", id = "
                  << obj->GetUniqueID() << endl;

            }
        }
    }

  //  //must reset kIsReferenced before clone,
  //  //otherwise pointer to obj_src in TProcessID table got altered
  //  Bool_t isRef = obj_src->TestBit(kIsReferenced);
  //  obj_src->ResetBit(kIsReferenced);
  //
  //  TObject * obj_clone = obj_src->Clone(name.c_str());
  //  obj_clone->ResetBit(kIsReferenced);
  //  obj_clone->ResetBit(kHasUUID);
  //  obj_clone->ResetBit(kMustCleanup);
  //  obj_clone->SetUniqueID(0);
  //
  //  if (isRef)
  //    obj_src->SetBit(kIsReferenced);

  TH1 * obj_clone = CleanClone(obj_src, name);

  if (Verbosity() >= 3)
    {
      // debug Clone

      cout << "saHist::" << _name << "::get_clone_hist - After clone" << endl;

      cout << "saHist::" << _name << "::get_clone_hist TProcessID = ";
      TProcessID::GetSessionProcessID()->Print();
      TObjArray * arr = TProcessID::GetSessionProcessID()->GetObjects();
      assert(arr);
      //      arr->Print();
      for (int i = 1; i < arr->GetSize(); i++)
        {
          const TObject * obj = arr->UncheckedAt(i);
          if (obj)
            cout << "\t arr @ " << i << " = " << obj->GetName() << ", id = "
                << obj->GetUniqueID() << endl;
        }
      cout << "saHist::" << _name << "::get_clone_hist obj_src = "
          << obj_src->GetName() << ", id = " << obj_src->GetUniqueID()
          << ", sum = " << static_cast<const TH1 *>(obj_src)->GetSumOfWeights()
          << endl;
      cout << "saHist::" << _name << "::get_clone_hist dest = "
          << obj_clone->GetName() << ", id = " << obj_clone->GetUniqueID()
          << ", sum = "
          << static_cast<const TH1 *>(obj_clone)->GetSumOfWeights() << endl;

      cout << "saHist::" << _name << "::get_clone_hist src TRef ID = "
          << ref.GetUniqueID() << endl;
      cout << "saHist::" << _name << "::get_clone_hist src TRef Object = "
          << ref.GetObject()->GetName() << ", id = "
          << ref.GetObject()->GetUniqueID() << ", sum = "
          << static_cast<const TH1 *>(ref.GetObject())->GetSumOfWeights()
          << endl;
      cout << "saHist::" << _name << "::get_clone_hist src TRef ID = "
          << ref.GetUniqueID() << endl;
    }

  return obj_clone;
}

void
saHist::Init(TH1 * h)
{

  if (!h)
    {
      throw runtime_error("saHist::Init - Constructing with invalid histogram");
    }

  cout << "saHist::Init with histogram " << h->GetName() << " in Verbosity "
      << Verbosity() << endl;

  if (isHistoRegistered("BASE"))
    {
      cout << "saHist::" << _name << "::Init - already initalized to " << endl;
      Print();

      throw runtime_error("saHist::Init - already initalized");
    }

  _name = (h->GetName());
  const string name(make_hist_name("saHist"));
  const string title("Histogram container for " + _name);

  SetName(name.c_str());
  SetTitle(title.c_str());

  h->SetName(make_hist_name("BASE").c_str());
  _hists["BASE"] = h;
  if (Verbosity() >= 2)
    {
      TProcessID::GetSessionProcessID()->Print();
      cout << "TProcessID->GetObjectCount() = "
          << TProcessID::GetSessionProcessID()->GetObjectCount() << endl;
      Print();
      gDirectory->Print();
    }

  insert_clone_hist("SUM");

  if (get_flag(CROSSING_CHECKS))
    {
      const string suffix = "CROSSING_CHECKS";

      if (!isHistoRegistered(suffix))
        {
          string t("RHIC crossing check for ");
          t += h->GetTitle();
          t += ";RHIC crossing";

          if (!get_flag(CROSSING_CHECKS_EACH_BIN))
            {

              _hists[suffix] = new TH1F(make_hist_name(suffix).c_str(),
                  t.c_str(), saEventProperty::N_CROSSING, -.5,
                  saEventProperty::N_CROSSING - .5);

            }
          else
            {
              t += ";Histogram Bins";

              int nbins = h->GetNbinsX() * h->GetNbinsY() * h->GetNbinsZ();

              _hists[suffix] = new TH2F(make_hist_name(suffix).c_str(),
                  t.c_str(), saEventProperty::N_CROSSING, -.5,
                  saEventProperty::N_CROSSING - .5, //
                  nbins, .5, nbins + .5);

            }

        }
    }

  if (get_flag(SPIN_DEPENDENT))
    {
      // for spin dependent fillings

      for (saEventProperty::sa_spin_state state = saEventProperty::PP;
          state <= saEventProperty::UNKNOWN; //
          state = static_cast<saEventProperty::sa_spin_state>(state + 1))
        {
          const string sufix = saEventProperty::get_state_name(state);
          insert_clone_hist(sufix);
        }
    }

  if (get_flag(POL_PRODUCTS))
    {
      // produce histograms of polarization products

      //! blue polarization
      insert_clone_hist("Pol_Blue");
      //! blue polarization
      insert_clone_hist("PolStatErr_Blue");

      //! yellow polarization
      insert_clone_hist("Pol_Yellow");
      //! yellow polarization
      insert_clone_hist("PolStatErr_Yellow");

      //! blue polarization *  yellow polarization
      insert_clone_hist("Pol2_Blue_Yellow");

    }

  if (Verbosity() >= 2)
    {
      cout << "saHist::" << _name << "::Init - initalized to ";
      Print();
    }
}

void
saHist::Print(Option_t *option) const
{
  TNamed::Print(option);

  string opt = string(option);
  const bool extend = opt.find("EXTEND") != string::npos;

  string tab = get_flag(SUB_HIST) ? "|--\t" : "";

  cout << tab << "--- Flag List ---" << endl;

  if (get_flag(SUB_HIST))
    cout << tab << "\t" << "This is a sub-histogram of another saHist object"
        << endl;

  if (get_flag(LUMI_HIST))
    cout << tab << "\t" << "This is for storing luminorsity information"
        << endl;

  if (get_flag(SPIN_DEPENDENT))
    cout << tab << "\t" << "Asymmetry analysis mode" << endl;
  else
    cout << tab << "\t" << "Non-Asymmetry analysis mode" << endl;

  if (get_flag(EVENT_PROPERTY))
    cout << tab << "\t" << "analyzed event by event" << endl;

  if (get_flag(RUN_PROPERTY))
    cout << tab << "\t" << "analyzed run by run" << endl;

  if (get_flag(FILL_PROPERTY))
    cout << tab << "\t" << "analyzed fill by fill" << endl;

  if (get_flag(POL_PRODUCTS))
    cout << tab << "\t" << "Produce products of polarizations" << endl;

  if (get_flag(CROSSING_CHECKS))
    cout << tab << "\t"
        << "Produce produce histograms to check RHIC crossing alignment"
        << endl;

  cout << tab << "--- Histogram List ---" << endl;
  for (sa_vec_hist::const_iterator hiter = _hists.begin();
      hiter != _hists.end(); ++hiter)
    {
      const TRef & ref = (*hiter).second;
      if (ref.IsValid() && dynamic_cast<TH1 *>(ref.GetObject()))
        cout << tab << "\t" << "hist[" << (*hiter).first << "] : "
            << ref.GetObject()->GetName() << ", id = "
            << ref.GetObject()->GetUniqueID() << ", sum = "
            << static_cast<const TH1 *>(ref.GetObject())->GetSumOfWeights()
            << endl;
      else
        cout << tab << "hist[" << (*hiter).first << "] is invalid" << endl;
    }

  if (nsaHists())
    {
      cout << tab << "--- Sub Spin Analysis Histogram List ---" << endl;
      for (sa_vec_hist::const_iterator hiter = _sahists.begin();
          hiter != _sahists.end(); ++hiter)
        {
          const TRef & ref = (*hiter).second;
          if (ref.IsValid() && dynamic_cast<saHist *>(ref.GetObject()))
            {
              cout << tab << "\t" << "sahist[" << (*hiter).first << "] : "
                  << ref.GetObject()->GetName() << ", id = "
                  << ref.GetObject()->GetUniqueID() << " containing "
                  << static_cast<const saHist *>(ref.GetObject())->nHistos()
                  << " histograms" << endl;

              if (extend)
                {
                  static_cast<const saHist *>(ref.GetObject())->Print(option);
                }

            }
          else
            cout << tab << "hist[" << (*hiter).first << "] is invalid" << endl;
        }
    }

}

bool
saHist::AutoLoad(int verbosity)
{
  bool all_good = true;

  try
    {

      if (verbosity < 0)
        verbosity = Verbosity();

      string tab = get_flag(SUB_HIST) ? "|--\t" : "";

      if (verbosity >= 2)
        cout << tab << "-------------------------------------" << endl;

      if (verbosity >= 2)
        cout << tab << "saHist::AutoLoad " << GetName() << " , containing "
            << nHistos() << " histograms" << endl;

      for (sa_vec_hist::const_iterator hiter = _hists.begin();
          hiter != _hists.end(); ++hiter)
        {
          if (hiter->second.GetObject())
            continue;

          TH1 * h = dynamic_cast<TH1 *>(gDirectory->Get(
              make_hist_name(hiter->first).c_str()));

          if (h && hiter->second.GetObject())
            {
              if (verbosity >= 2)
                cout << tab << "\t" << "hist[" << (*hiter).first
                    << "] loaded successfully, " << h->GetName() << ", sum = "
                    << h->GetSumOfWeights() << endl;
            }
          else
            {
              cout << tab << "\t" << "hist[" << (*hiter).first
                  << "] cannot be loaded" << endl;
              all_good = false;
            }
        }

      for (sa_vec_hist::iterator hiter = _hists.begin(); hiter != _hists.end();
          ++hiter)
        {
          TRef & ref = (*hiter).second;
          TH1 * obj_src = dynamic_cast<TH1 *>(ref.GetObject());

          if (obj_src)
            {
              if (obj_src->GetDirectory())
                {
                  // decouple histograms and directories

                  //              ref= (NULL);
                  //              TH1 * h = CleanClone(obj_src, obj_src->GetName());
                  //              delete obj_src;
                  //
                  //              ref= (h);

                  ref = (NULL);

                  obj_src->SetDirectory(NULL);
                  obj_src->ResetBit(kIsReferenced);
                  obj_src->ResetBit(kHasUUID);
                  obj_src->ResetBit(kMustCleanup);
                  obj_src->SetUniqueID(0);

                  ref = obj_src;

                  if (verbosity >= 2)
                    cout << tab << "\t" << "hist[" << (*hiter).first
                        << "] refreshed successfully, " << obj_src->GetName()
                        << ", sum = " << obj_src->GetSumOfWeights()
                        //                    <<" - "<<ref.GetObject()
                        << endl;
                }

            }
          else
            {
              cout << tab << "hist[" << (*hiter).first
                  << "] cannot be refreshed" << endl;
              all_good = false;
            }
        }

      if (verbosity >= 2)
        cout << tab << "saHist::AutoLoad " << GetName() << " , containing "
            << nsaHists() << " spin analysis histograms" << endl;

      for (sa_vec_hist::const_iterator hiter = _sahists.begin();
          hiter != _sahists.end(); ++hiter)
        {
          if (hiter->second.GetObject())
            continue;

          saHist * h = dynamic_cast<saHist *>(gDirectory->Get(
              make_hist_name(hiter->first + "_saHist").c_str()));

          if (h && hiter->second.GetObject())
            {
              if (verbosity >= 2)
                cout << tab << "sahist[" << (*hiter).first
                    << "] loaded successfully, " << h->GetName()
                    << ", containing " << h->nHistos() << " histograms and "
                    << h->nsaHists() << " spin analysis histograms" << endl;

              h->AutoLoad(verbosity - 1);
            }
          else
            {
              cout << tab << "sahist[" << (*hiter).first << "] cannot be loaded"
                  << endl;
              all_good = false;
            }
        }

      // rebuild clean reference to sub saHist
      for (sa_vec_hist::iterator hiter = _sahists.begin();
          hiter != _sahists.end(); ++hiter)
        {
          TRef & ref = (*hiter).second;
          saHist * obj_src = dynamic_cast<saHist *>(ref.GetObject());

          if (obj_src)
            {
              ref = (NULL);

              obj_src->ResetBit(kIsReferenced);
              obj_src->ResetBit(kHasUUID);
              obj_src->ResetBit(kMustCleanup);
              obj_src->SetUniqueID(0);

              ref = obj_src;

              if (verbosity >= 2)
                cout << tab << "sahist[" << (*hiter).first
                    << "] refreshed successfully, " << obj_src->GetName()
                    << ", containing " << obj_src->nHistos()
                    << " histograms and " << obj_src->nsaHists()
                    << " spin analysis histograms" << endl;

            }
          else
            {
              cout << tab << "sahist[" << (*hiter).first
                  << "] cannot be refreshed" << endl;
              all_good = false;
            }
        }
    }
  catch (exception& e)
    {
      cout << "saHist::AutoLoad - Error - Failed with exception: " << e.what() << '\n';
      all_good = false;
      return false;
    }

  return all_good;
}

//! clone without dirty association mess with TRef and TDirectory
TH1 *
saHist::CleanClone(TH1 * obj_src, const std::string & name)
{

  if (!obj_src)
    return NULL;

  //must reset kIsReferenced before clone,
  //otherwise pointer to obj_src in TProcessID table got altered
  Bool_t isRef = obj_src->TestBit(kIsReferenced);
  TDirectory * dir = obj_src->GetDirectory();

  obj_src->ResetBit(kIsReferenced);
  obj_src->SetDirectory(NULL);

  TH1 * obj_clone = static_cast<TH1 *>(obj_src->Clone(name.c_str()));
  obj_clone->ResetBit(kIsReferenced);
  obj_clone->ResetBit(kHasUUID);
  obj_clone->ResetBit(kMustCleanup);
  obj_clone->SetUniqueID(0);
  obj_clone->SetDirectory(NULL);

  if (isRef)
    obj_src->SetBit(kIsReferenced);
  obj_src->SetDirectory(dir);

  return obj_clone;
}

Int_t
saHist::Write(const char *name, Int_t option, Int_t bufsize) const
{
  if (name)
    cout << "saHist::" << _name << "Write - WARNING - rename to " << name
        << ", which may not reload correctly from root files." << endl;

  Int_t n_write = 0;

  sa_vec_hist::const_iterator hiter;
  for (hiter = _hists.begin(); hiter != _hists.end(); ++hiter)
    {
      const TRef & ref = (*hiter).second;

      if (ref.IsValid() && ref.GetObject())
        n_write += ref.GetObject()->Write(name, option, bufsize);
    }

  for (hiter = _sahists.begin(); hiter != _sahists.end(); ++hiter)
    {
      const TRef & ref = (*hiter).second;

      if (ref.IsValid() && ref.GetObject())
        n_write += ref.GetObject()->Write(NULL, option, bufsize);
    }

  n_write += TNamed::Write(name, option, bufsize);

  if (Verbosity() >= 1)
    {
      cout << "saHist::" << _name << "::Write - Finished writing with "
          << n_write << "B" << endl;
      //      Print();
      //      gDirectory->Print();
    }

  return n_write;
}

Int_t
saHist::Fill(double v1, double v2, double v3, double v4, double extra_weight)
{

  PHCompositeNode * top_node = Fun4AllServer::instance()->topNode("TOP");

  if (!top_node)
    {
      cout << "saHist::" << _name
          << "::Fill - Cannot find top node in Fun4AllServer" << endl;
      throw std::runtime_error(
      DESCRIPTION("saHist::Fill - Cannot find top node in Fun4AllServer"));
    }

  const saEventProperty * ep = findNode::getClass<saEventProperty>(top_node,
      "saEventProperty");
  if (!ep)
    {
      cout << "saHist::" << _name
          << "::Fill - Cannot find EventProperty node in Top Node" << endl;
      throw std::runtime_error(
      DESCRIPTION("saHist::Fill - Cannot find EventProperty node in Top Node"));
    }

  return Fill(ep, v1, v2, v3, v4, extra_weight);
}

//! fill the histogram
Int_t
saHist::Fill(const saEventProperty * ep, double v1, double v2, double v3,
    double v4, double extra_weight)
{

  if (Verbosity() >= 2)
    {
      cout << "saHist:" << _name << "::Fill - " << "[ " << v1 << ", " << v2
          << ", " << v3 << ", " << v4 << "] in ";
      Print();
    }

  double n_fill = 0;
  double weight = 0;

  TH1 * h_sum = getHisto("SUM");
  if (h_sum)
    weight = Fill(h_sum, v1, v2, v3, v4, extra_weight);
  else
    {
      cout << "saHist::" << _name << "::Fill - Error - missing base histogram."
          << endl;
      return n_fill;
    }
  n_fill += weight;

  if (get_flag(CROSSING_CHECKS))
    {
      if (!ep)
        {
          cout << "saHist::" << _name
              << "::Fill - Error - invalid saEventProperty" << endl;
          return n_fill;
        }

      const int crossing = ep->get_crossing_id_RHIC();

      if (Verbosity() >= 2)
        cout << "saHist::" << _name << "::Fill - INFO - crossing " << crossing
            << " filled with " << weight << endl;

      if (!get_flag(CROSSING_CHECKS_EACH_BIN))
        {
          TH1F * h_spin = dynamic_cast<TH1F *>(getHisto("CROSSING_CHECKS"));
          if (h_spin)
            h_spin->Fill(crossing, weight);
        }
      else
        {

          const int nx = h_sum->GetXaxis()->FindBin(v1);
          const int ny = h_sum->GetYaxis()->FindBin(v2);
          const int nz = h_sum->GetZaxis()->FindBin(v3);

          const int bin = h_sum->GetBin(nx, ny, nz);

          TH2F * h_spin = dynamic_cast<TH2F *>(getHisto("CROSSING_CHECKS"));
          if (h_spin)
            h_spin->Fill(crossing, bin, weight);

        }

    }

  if (get_flag(SPIN_DEPENDENT))
    {
      if (!ep)
        {
          cout << "saHist::" << _name
              << "::Fill - Error - invalid saEventProperty" << endl;
          return n_fill;
        }

      saEventProperty::sa_spin_state state = ep->get_spin_state();
      const string name = saEventProperty::get_state_name(state);

      if (Verbosity() >= 2)
        {
          cout << "saHist:" << _name << "::Fill - spin state " << name << " = "
              << (int) state << " with ";
          ep->identify();
        }

      TH1 * h_spin = getHisto(name);
      if (h_spin)
        {
          n_fill += Fill(h_spin, v1, v2, v3, v4, extra_weight);
        }
      else
        {
          cout << "saHist::" << _name << "::Fill - Error - missing histogram "
              << name << endl;
        }
    }

  if (get_flag(POL_PRODUCTS))
    {
      if (!ep)
        {
          cout << "saHist::" << _name
              << "::Fill - Error - invalid saEventProperty" << endl;
          return n_fill;
        }

      if (Verbosity() >= 2)
        {
          cout << "saHist:" << _name << "::Fill - polarizations "
              //
              << ", Blue = " << ep->get_polarization_blue() << " +/- "
              << ep->get_polarization_stat_err_blue()
              //
              << ", Yellow = " << ep->get_polarization_yellow() << " +/- "
              << ep->get_polarization_stat_err_yellow() << endl;
        }

      n_fill += Fill("Pol_Blue", v1, v2, v3, v4,
          extra_weight * ep->get_polarization_blue());
      n_fill += Fill("Pol_Yellow", v1, v2, v3, v4,
          extra_weight * ep->get_polarization_yellow());
      n_fill += Fill("Pol2_Blue_Yellow", v1, v2, v3, v4,
          extra_weight * ep->get_polarization_blue()
              * ep->get_polarization_yellow());

      n_fill += Fill("PolStatErr_Blue", v1, v2, v3, v4,
          extra_weight * ep->get_polarization_stat_err_blue());
      n_fill += Fill("PolStatErr_Yellow", v1, v2, v3, v4,
          extra_weight * ep->get_polarization_stat_err_yellow());

    }

  if (get_flag(RUN_PROPERTY))
    {
      if (!ep)
        {
          cout << "saHist::" << _name
              << "::Fill - Error - invalid saEventProperty" << endl;
          return n_fill;
        }

      const int run = ep->get_run_number();
      if (run < 0)
        {
          cout << "saHist:" << _name << "::Fill - ERROR - invalid run number "
              << run << endl;

          return n_fill;
        }

      if (Verbosity() >= 2)
        {
          cout << "saHist:" << _name << "::Fill - INFO - run " << run << endl;
        }

      saHist * sahist = getsaHist_run(static_cast<unsigned int>(run), true);
      if (sahist)
        sahist->Fill(ep, v1, v2, v3, v4, extra_weight);
    }

  if (get_flag(FILL_PROPERTY))
    {
      if (!ep)
        {
          cout << "saHist::" << _name
              << "::Fill - Error - invalid saEventProperty" << endl;
          return n_fill;
        }

      const int fill = ep->get_fill_number();
      if (fill < 0)
        {
          cout << "saHist:" << _name << "::Fill - ERROR - invalid fill number "
              << fill << endl;

          return n_fill;
        }

      if (Verbosity() >= 2)
        {
          cout << "saHist:" << _name << "::Fill - INFO - fill " << fill << endl;
        }

      saHist * sahist = getsaHist_fill(static_cast<unsigned int>(fill), true);
      if (sahist)
        sahist->Fill(ep, v1, v2, v3, v4, extra_weight);
    }

  return n_fill;
}

double
saHist::Fill(const std::string &hname, double v1, double v2, double v3,
    double v4, double extra_weight) const
{
  TH1 * h = saHist::getHisto(hname);

  if (!h)
    {
      cout << "saHist::" << _name << "::Fill - Error - can not find histogram "
          << hname << endl;

      return 0;
    }

  return Fill(h, v1, v2, v3, v4, extra_weight);
}

double
saHist::Fill(TH1 * h, double v1, double v2, double v3, double v4,
    double extra_weight)
{

  if (!h)
    {
      cout << "saHist::Fill - Error - Empty histogram" << endl;
      return 0;
    }

  static const char * negative_weight_warning = //
      "You are filling a negative weight to this histogram. \n"
          "In this case, the distribution for the content of this histogram is not Poission anymore.\n"
          "The default SpinAnalyzer error calculation did assumed independent Poission error distribution for the histograms."
          "Therefore, Please avoid using the built-in error calculations, e.g. saHist::CalcAsymmetry().\n"
          "And you are responsible for handle the error propagation now. \n";

  if (h->GetDimension() == 3)
    {
      (static_cast<TH3 *>(h))->Fill(v1, v2, v3, v4 * extra_weight);

      static bool once = true;
      if (once && v4 * extra_weight < 0)
        {
          once = false;

          cout << "saHist::" << h->GetName() << "::Fill - WARNING - "
              << negative_weight_warning << endl;
        }

      return v4 * extra_weight;
    }
  else if (h->GetDimension() == 2)
    {
      (static_cast<TH2 *>(h))->Fill(v1, v2, v3 * extra_weight);

      static bool once = true;
      if (once && v3 * extra_weight < 0)
        {
          once = false;

          cout << "saHist::" << h->GetName() << "::Fill - WARNING - "
              << negative_weight_warning << endl;
        }

      return v3 * extra_weight;
    }
  else if (h->GetDimension() == 1)
    {
      (static_cast<TH1 *>(h))->Fill(v1, v2 * extra_weight);

      static bool once = true;
      if (once && v2 * extra_weight < 0)
        {
          once = false;

          cout << "saHist::" << h->GetName() << "::Fill - WARNING - "
              << negative_weight_warning << endl;
        }

      return v2 * extra_weight;
    }
  else
    {
      cout << "saHist::" << h->GetName() << "::Fill - Error - Zombie object:";
      h->Class()->Print();
      cout << endl;

      return 0;
    }

  assert(0);
  //should never reach here;
  return 0;
}

int
saHist::isHistoRegistered(const std::string &name) const
{
  sa_vec_hist::const_iterator histoiter = _hists.find(name);
  if (histoiter != _hists.end())
    {
      return 1;
    }
  return 0;
}

TH1 *
saHist::getHisto(const unsigned int ihisto) const
{
  sa_vec_hist::const_iterator histoiter = _hists.begin();
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
      return static_cast<TH1 *>(histoiter->second.GetObject());
    }
  else
    {
      cout << "saHist::" << _name
          << "::getHisto: ERROR Invalid histogram number: " << ihisto
          << ", maximum number is " << size << endl;
    }
  return 0;
}

const char *
saHist::getHistoName(const unsigned int ihisto) const
{
  sa_vec_hist::const_iterator histoiter = _hists.begin();
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
      cout << "saHist::" << _name
          << "::getHisto: ERROR Invalid histogram number: " << ihisto
          << ", maximum number is " << size << endl;
    }
  return 0;
}

TH1 *
saHist::getHisto(const string &suffix) const
{
  sa_vec_hist::const_iterator histoiter = _hists.find(suffix);
  if (histoiter != _hists.end())
    {
      return static_cast<TH1 *>(histoiter->second.GetObject());
    }

  if (Verbosity() > 1)
    {
      cout << "saHist::" << _name << "::getHisto: ERROR Unknown Histogram "
          << suffix << ", The following are implemented: " << endl;
      Print("ALL");
    }

  return 0;
}

//! remove a histogram from the list
bool
saHist::removeHisto(const std::string& suffix)
{

  sa_vec_hist::iterator histoiter = _hists.find(suffix);
  if (histoiter != _hists.end())
    {

      //      TRef & ref = histoiter->second;

      if (Verbosity() > 1)
        {
          cout << "saHist::" << _name << "::removeHisto: erase Histogram "
              << suffix << endl;
        }

      //      if (ref.IsValid() && ref.GetObject())
      //        delete ref.GetObject();

      _hists.erase(histoiter);

      if (Verbosity() > 2)
        {
          cout << "saHist::" << _name << "::removeHisto: remaining histos: ";
          for (sa_vec_hist::const_iterator hiter = _hists.begin();
              hiter != _hists.end(); ++hiter)
            {
              const string name = (*hiter).first;

              cout << name << ", ";

            }
          cout << endl;
        }

      return true;
    }

  if (Verbosity() > 1)
    {
      cout << "saHist::" << _name << "::removeHisto: ERROR Unknown Histogram "
          << suffix << endl;
    }
  return false;
}

TH1 *
saHist::getHisto(const saEventProperty::sa_spin_state state) const
{
  if (!get_flag(SPIN_DEPENDENT))
    {
      cout << "saHist::" << _name
          << "::getHisto - Error - try to access spin dependent histogram for non-spin-dependent object "
          << _name << endl;

      return NULL;
    }

  return getHisto(saEventProperty::get_state_name(state));
}

unsigned int
saHist::nHistos() const
{
  return _hists.size();
}

int
saHist::issaHistRegistered(const std::string &name) const
{
  sa_vec_hist::const_iterator histoiter = _sahists.find(name);
  if (histoiter != _sahists.end())
    {
      return 1;
    }
  return 0;
}

saHist *
saHist::getsaHist(const unsigned int ihisto) const
{
  sa_vec_hist::const_iterator histoiter = _sahists.begin();
  unsigned int size = _sahists.size();
  if (Verbosity() > 3)
    {
      cout << "saHist::" << _name << "::getsaHist:Map contains " << size
          << " Elements" << endl;
    }
  if (ihisto < size)
    {
      for (unsigned int i = 0; i < ihisto; i++)
        {
          ++histoiter;
        }
      return static_cast<saHist *>(histoiter->second.GetObject());
    }
  else
    {
      cout << "saHist::" << _name
          << "::getsaHist: ERROR Invalid histogram number: " << ihisto
          << ", maximum number is " << size << endl;
    }
  return 0;
}

const char *
saHist::getsaHistName(const unsigned int ihisto) const
{
  sa_vec_hist::const_iterator histoiter = _sahists.begin();
  unsigned int size = _sahists.size();
  if (verbosity > 3)
    {
      cout << "saHist::" << _name << "::getsaHist:Map contains " << size
          << " Elements" << endl;
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
      cout << "saHist::" << _name
          << "::getsaHist: ERROR Invalid histogram number: " << ihisto
          << ", maximum number is " << size << endl;
    }
  return 0;
}

saHist *
saHist::getsaHist(const string &suffix, bool add)
{
  sa_vec_hist::const_iterator histoiter = _sahists.find(suffix);
  if (histoiter != _sahists.end())
    {
      return static_cast<saHist *>(histoiter->second.GetObject());
    }

  if (add)
    {
      cout << "saHist::" << _name
          << "::getsaHist - INFO - add new sub spin analysis histogram with suffix "
          << suffix << endl;

      insert_clone_sahist(suffix);

      return getsaHist(suffix, false);
    }
  else
    {
      if (Verbosity() > 1)
        {
          cout << "saHist::" << _name << "::getsaHist: ERROR Unknown Histogram "
              << suffix << ", The following are implemented: " << endl;
          Print("ALL");
        }
    }

  return 0;
}

unsigned int
saHist::nsaHists() const
{
  return _sahists.size();
}

saHist::n_list
saHist::get_run_list(void)
{

  n_list list;

  for (sa_vec_hist::const_iterator hiter = _sahists.begin();
      hiter != _sahists.end(); ++hiter)
    {
      const string name = (*hiter).first;

      int number = 0;

      if (sscanf(name.c_str(), "RUN_%d", &number))
        {
          list.push_back(number);
        }
    }

  return list;
}

saHist::n_list
saHist::get_fill_list(void)
{

  n_list list;

  for (sa_vec_hist::const_iterator hiter = _sahists.begin();
      hiter != _sahists.end(); ++hiter)
    {
      const string name = (*hiter).first;

      int number = 0;

      if (sscanf(name.c_str(), "FILL_%d", &number))
        {
          list.push_back(number);
        }
    }

  return list;
}

void
saHist::RemoveAsymmetryHisto()
{

  if (Verbosity() >= 2)
    cout << "saHist::" << _name
        << "::RemoveAsymmetryHisto - INFO - remove asymmetry histograms for"
        << GetName() << endl;

  removeHisto("A_LL");
  removeHisto("A_L_Blue");
  removeHisto("A_L_Yellow");
  removeHisto("YIELD");

  // process sub histogram
  for (sa_vec_hist::const_iterator hiter = _sahists.begin();
      hiter != _sahists.end(); ++hiter)
    {

      const TRef & ref = hiter->second;

      saHist * sub_hist = dynamic_cast<saHist *>(ref.GetObject());

      if (ref.IsValid() && sub_hist)
        {

          sub_hist->RemoveAsymmetryHisto();

        }
    }
}

//! Calculate asymmetry according to luminorsity histogram h_lumi
bool
saHist::CalcAsymmetry(saHist * h_lumi, enu_lumi_id lumi_type)
{
  RemoveAsymmetryHisto();

  bool status = true;

  status &= CalcAsymmetry_Work(h_lumi, lumi_type);

  // process sub histogram
  for (sa_vec_hist::const_iterator hiter = _sahists.begin();
      hiter != _sahists.end(); ++hiter)
    {

      const string & suffix = hiter->first;
      const TRef & ref = hiter->second;

      saHist * sub_hist = dynamic_cast<saHist *>(ref.GetObject());

      if (ref.IsValid() && sub_hist)
        {
          if (Verbosity() >= 1)
            cout << "saHist::" << _name
                << "::CalcAsymmetry - INFO - process sub histogram "
                << ref.GetObject()->GetName() << endl;

          saHist * sub_lumi = h_lumi->getsaHist(suffix);

          if (!sub_lumi)
            {
              cout << "saHist::" << _name
                  << "::CalcAsymmetry - ERROR - can not find lumi info for "
                  << suffix << endl;

              status = false;

              continue;
            }

          status &= sub_hist->CalcAsymmetry_Work(sub_lumi, lumi_type);

        }
    }

  return status;
}

//! Calculate asymmetry according to luminorsity histogram h_lumi, DO NOT iterate through sub-histograms
bool
saHist::CalcAsymmetry_Work(saHist * h_lumi, enu_lumi_id lumi_type)
{

  if (!get_flag(SPIN_DEPENDENT))
    {
      cout << "saHist::" << _name << "::CalcAsymmetry - Error - " << _name
          << " is not set up for spin analysis." << endl;
      return false;
    }

  if ( //
  !(this->getHisto("BASE")) or //
      !(this->getHisto("PP")) or //
      !(this->getHisto("MP")) or //
      !(this->getHisto("PM")) or //
      !(this->getHisto("MM")) //
      )
    {
      cout << "saHist::" << _name
          << "::CalcAsymmetry - Error - invalid spin info." << endl;
      this->Print();
      return false;
    }

  if (!h_lumi)
    {
      cout << "saHist::" << _name
          << "::CalcAsymmetry - Error - missing luminorsity info." << endl;
      return false;
    }

  if ( //
  !(h_lumi->getHisto("BASE")) or //
      !(h_lumi->getHisto("PP")) or //
      !(h_lumi->getHisto("MP")) or //
      !(h_lumi->getHisto("PM")) or //
      !(h_lumi->getHisto("MM")) or //
      !(h_lumi->getHisto("Pol2_Blue_Yellow")) or //
      !(h_lumi->getHisto("Pol_Blue")) or //
      !(h_lumi->getHisto("Pol_Yellow")) //
      )
    {
      cout << "saHist::" << _name
          << "::CalcAsymmetry - Error - invalid luminorsity info:" << endl;
      h_lumi->Print();
      return false;
    }

  if (h_lumi->getHisto("BASE")->GetNbinsX() != LUMI_COUNT)
    {
      cout << "saHist::" << _name
          << "::CalcAsymmetry - Error - invalid luminorsity format" << endl;
      return false;
    }

  cout << "saHist::" << _name << "::CalcAsymmetry with lumi"
      << h_lumi->get_name() << "type = " << lumi_type << endl;

  insert_clone_hist("A_LL");
  insert_clone_hist("A_L_Blue");
  insert_clone_hist("A_L_Yellow");
  insert_clone_hist("YIELD");

  const TH1 * h_NPP = getHisto("PP");
  const TH1 * h_NPM = getHisto("PM");
  const TH1 * h_NMP = getHisto("MP");
  const TH1 * h_NMM = getHisto("MM");

  const TH1 * h_RPP = h_lumi->getHisto("PP");
  const TH1 * h_RPM = h_lumi->getHisto("PM");
  const TH1 * h_RMP = h_lumi->getHisto("MP");
  const TH1 * h_RMM = h_lumi->getHisto("MM");
  const TH1 * h_RSUM = h_lumi->getHisto("SUM");

  const double RPP = h_RPP->GetBinContent(lumi_type);
  const double RMP = h_RMP->GetBinContent(lumi_type);
  const double RPM = h_RPM->GetBinContent(lumi_type);
  const double RMM = h_RMM->GetBinContent(lumi_type);
  const double RSUM = h_RSUM->GetBinContent(lumi_type);

  const TH1 * h_Pol2_Blue_Yellow = h_lumi->getHisto("Pol2_Blue_Yellow");
  const TH1 * h_Pol_Blue = h_lumi->getHisto("Pol_Blue");
  const TH1 * h_Pol_Yellow = h_lumi->getHisto("Pol_Yellow");

  const double POL2_BLUE_YELLOW = h_Pol2_Blue_Yellow->GetBinContent(lumi_type)
      / RSUM;
  const double POL_BLUE = h_Pol_Blue->GetBinContent(lumi_type) / RSUM;
  const double POL_YELLOW = h_Pol_Yellow->GetBinContent(lumi_type) / RSUM;

  cout << "saHist::" << _name << "::CalcAsymmetry polarization (B, Y, B*Y) = " //
      << POL_BLUE << ", " //
      << POL_YELLOW << ", " //
      << POL2_BLUE_YELLOW << ")" << endl;

  TH1 * h_A_LL = getHisto("A_LL");
  TH1 * h_A_L_Blue = getHisto("A_L_Blue");
  TH1 * h_A_L_Yellow = getHisto("A_L_Yellow");
  TH1 * h_YIELD = getHisto("YIELD");

  const int n_bin_x = getHisto("BASE")->GetNbinsX();
  const int n_bin_y = getHisto("BASE")->GetNbinsY();
  const int n_bin_z = getHisto("BASE")->GetNbinsZ();

  for (int x = 1; x <= n_bin_x; x++)
    for (int y = 1; y <= n_bin_y; y++)
      for (int z = 1; z <= n_bin_z; z++)
        {
          /*
           * put this model in mathematica, and copy whatever I get.
           *
           NPP == RPP*Y*(1 + AB + ALL + AY)
           NMP == RMP*Y*(1 - AB - ALL + AY)
           NPM == RPM*Y*(1 + AB - ALL - AY)
           NMM == RMM*Y*(1 - AB + ALL - AY)

           you are welcome to check result,..
           */

          const double NPP = h_NPP->GetBinContent(x, y, z);
          const double NMP = h_NMP->GetBinContent(x, y, z);
          const double NPM = h_NPM->GetBinContent(x, y, z);
          const double NMM = h_NMM->GetBinContent(x, y, z);

          double A_LL = ((NPP * RMM * RMP * RPM - NPM * RMM * RMP * RPP
              - NMP * RMM * RPM * RPP + NMM * RMP * RPM * RPP)
              / (NPP * RMM * RMP * RPM + NPM * RMM * RMP * RPP
                  + NMP * RMM * RPM * RPP + NMM * RMP * RPM * RPP));

          double A_L_Blue = ((NPP * RMM * RMP * RPM + NPM * RMM * RMP * RPP
              - NMP * RMM * RPM * RPP - NMM * RMP * RPM * RPP)
              / (NPP * RMM * RMP * RPM + NPM * RMM * RMP * RPP
                  + NMP * RMM * RPM * RPP + NMM * RMP * RPM * RPP));

          double A_L_Yellow = ((NPP * RMM * RMP * RPM - NPM * RMM * RMP * RPP
              + NMP * RMM * RPM * RPP - NMM * RMP * RPM * RPP)
              / (NPP * RMM * RMP * RPM + NPM * RMM * RMP * RPP
                  + NMP * RMM * RPM * RPP + NMM * RMP * RPM * RPP));

          double YIELD = ((NPP * RMM * RMP * RPM + NPM * RMM * RMP * RPP
              + NMP * RMM * RPM * RPP + NMM * RMP * RPM * RPP)
              / (4. * RMM * RMP * RPM * RPP));

          double A_LL_ERR2 = (4 * Power(RMM, 2) * Power(RMP, 2) * Power(RPM, 2)
              * Power(RPP, 2)
              * (NPP * Power(RMM, 2)
                  * (NPM * (NPM + NPP) * Power(RMP, 2)
                      + 2 * NMP * NPM * RMP * RPM
                      + NMP * (NMP + NPP) * Power(RPM, 2))
                  + 2 * NMM * NPP * RMM
                      * (NPM * Power(RMP, 2) + NMP * Power(RPM, 2)) * RPP
                  + NMM
                      * (Power(NPM * RMP + NMP * RPM, 2)
                          + NMM * (NPM * Power(RMP, 2) + NMP * Power(RPM, 2)))
                      * Power(RPP, 2)))
              / Power(
                  NPP * RMM * RMP * RPM
                      + (NPM * RMM * RMP + NMP * RMM * RPM + NMM * RMP * RPM)
                          * RPP, 4);

          double A_L_Blue_ERR2 = (4 * Power(RMM, 2) * Power(RMP, 2)
              * Power(RPM, 2) * Power(RPP, 2)
              * (NPP
                  * (NMP * (NMP + NPP) * Power(RMM, 2)
                      + 2 * NMM * NMP * RMM * RMP
                      + NMM * (NMM + NPP) * Power(RMP, 2)) * Power(RPM, 2)
                  + 2 * NPM * NPP * (NMP * Power(RMM, 2) + NMM * Power(RMP, 2))
                      * RPM * RPP
                  + NPM
                      * (NMP * (NMP + NPM) * Power(RMM, 2)
                          + 2 * NMM * NMP * RMM * RMP
                          + NMM * (NMM + NPM) * Power(RMP, 2)) * Power(RPP, 2)))
              / Power(
                  NPP * RMM * RMP * RPM
                      + (NPM * RMM * RMP + NMP * RMM * RPM + NMM * RMP * RPM)
                          * RPP, 4);

          double A_L_Yellow_ERR2 = (4 * Power(RMM, 2) * Power(RMP, 2)
              * Power(RPM, 2) * Power(RPP, 2)
              * (NPP * Power(RMP, 2)
                  * (NPM * (NPM + NPP) * Power(RMM, 2)
                      + 2 * NMM * NPM * RMM * RPM
                      + NMM * (NMM + NPP) * Power(RPM, 2))
                  + 2 * NMP * NPP * RMP
                      * (NPM * Power(RMM, 2) + NMM * Power(RPM, 2)) * RPP
                  + NMP
                      * (Power(NPM * RMM + NMM * RPM, 2)
                          + NMP * (NPM * Power(RMM, 2) + NMM * Power(RPM, 2)))
                      * Power(RPP, 2)))
              / Power(
                  NPP * RMM * RMP * RPM
                      + (NPM * RMM * RMP + NMP * RMM * RPM + NMM * RMP * RPM)
                          * RPP, 4);

          double YIELD_ERR2 = (NMM / Power(RMM, 2) + NMP / Power(RMP, 2)
              + NPM / Power(RPM, 2) + NPP / Power(RPP, 2)) / 16.;

          //average polarization scale
          A_LL /= POL2_BLUE_YELLOW;
          A_L_Blue /= POL_BLUE;
          A_L_Yellow /= POL_YELLOW;

          A_LL_ERR2 /= POL2_BLUE_YELLOW * POL2_BLUE_YELLOW;
          A_L_Blue_ERR2 /= POL_BLUE * POL_BLUE;
          A_L_Yellow_ERR2 /= POL_YELLOW * POL_YELLOW;

          //nan numbers handling
          if (isnan(A_LL) or isnan(A_LL_ERR2))
            {
              A_LL = 0;
              A_LL_ERR2 = 1;
            }
          if (isnan(A_L_Blue) or isnan(A_L_Blue_ERR2))
            {
              A_L_Blue = 0;
              A_L_Blue_ERR2 = 1;
            }
          if (isnan(A_L_Yellow) or isnan(A_L_Yellow_ERR2))
            {
              A_L_Yellow = 0;
              A_L_Yellow_ERR2 = 1;
            }
          if (isnan(YIELD) or isnan(YIELD_ERR2))
            {
              YIELD = 0;
              YIELD_ERR2 = 0;
            }

          // save results
          h_A_LL->SetBinContent(x, y, z, A_LL);
          h_A_L_Blue->SetBinContent(x, y, z, A_L_Blue);
          h_A_L_Yellow->SetBinContent(x, y, z, A_L_Yellow);
          h_YIELD->SetBinContent(x, y, z, YIELD);

          h_A_LL->SetBinError(x, y, z, sqrt(A_LL_ERR2));
          h_A_L_Blue->SetBinError(x, y, z, sqrt(A_L_Blue_ERR2));
          h_A_L_Yellow->SetBinError(x, y, z, sqrt(A_L_Yellow_ERR2));
          h_YIELD->SetBinError(x, y, z, sqrt(YIELD_ERR2));

        }

  return true;
}

//! Calculate final asymmetry by fit run-by-run or fill-by-fill asymmetry results from the sub saHits
//! \param[in] run_or_fill true : use the run list, false : use the fill list
//! \return return a new sub saHist, storing the new asymmetry fit results
saHist *
saHist::CalcAsymmetry_Chi2Fit(bool run_or_fill)
{
  string suffix_sahist = run_or_fill ? "RUN_FIT" : "FILL_FIT";

  // build new saHist to save result
  saHist * sahAsym = getsaHist(suffix_sahist, false);
  if (!sahAsym)
    {

      TH1 * new_base = static_cast<TH1 *>(get_clone_hist(
          make_hist_name(suffix_sahist)));

      sa_flags_t new_flags = get_flags();

      new_flags &= ~RUN_PROPERTY;
      new_flags &= ~FILL_PROPERTY;
      new_flags &= ~EVENT_PROPERTY;
      new_flags &= ~SPIN_DEPENDENT;
      new_flags &= ~CROSSING_CHECKS;

      new_flags |= SUB_HIST;

      sahAsym = new saHist(new_base, new_flags, Verbosity());
      _sahists[suffix_sahist] = sahAsym;
    }

  assert(sahAsym);

  CalcAsymmetry_Chi2Fit_Work(run_or_fill, sahAsym, "A_LL");
  CalcAsymmetry_Chi2Fit_Work(run_or_fill, sahAsym, "A_L_Blue");
  CalcAsymmetry_Chi2Fit_Work(run_or_fill, sahAsym, "A_L_Yellow");
  CalcAsymmetry_Chi2Fit_Work(run_or_fill, sahAsym, "YIELD");

  if (Verbosity())
    {

      cout << "saHist::" << get_name()
          << "::CalcAsymmetry_Chi2Fit - Rebuild asymmetry result using "
          << (run_or_fill ? "run-by-run" : "fill-by-fill") //
          << " fits. Results saved to:" << endl;

      sahAsym->Print();
    }

  return sahAsym;
}

bool
saHist::CalcAsymmetry_Chi2Fit_Work(bool run_or_fill, saHist * sahAsym,
    TString hist_type)
{
  assert(sahAsym);

  n_list src_list;

  if (run_or_fill)
    src_list = get_run_list();
  else
    src_list = get_fill_list();

  TGraphErrors * ge_type = new TGraphErrors(src_list.size());

  sahAsym->insert_clone_hist((hist_type + "_FitPar").Data());
  sahAsym->insert_clone_hist((hist_type + "_Chi2").Data());
  sahAsym->insert_clone_hist((hist_type + "_Ndf").Data());

  TH1 *hFitPar = sahAsym->getHisto((hist_type + "_FitPar").Data());
  TH1 *hChi2 = sahAsym->getHisto((hist_type + "_Chi2").Data());
  TH1 *hNdf = sahAsym->getHisto((hist_type + "_Ndf").Data());

  //TH1 *hFitPar = CleanClone(getHisto(hist_type.Data()),hist_type.Data());
  //TH1 *hChi2 = CleanClone(getHisto(hist_type.Data()),hist_type.Data());
  //TH1 *hNdf = CleanClone(getHisto(hist_type.Data()),hist_type.Data());

  const int n_bin_x = getHisto("BASE")->GetNbinsX();
  const int n_bin_y = getHisto("BASE")->GetNbinsY();
  const int n_bin_z = getHisto("BASE")->GetNbinsZ();

  for (int x = 1; x <= n_bin_x; x++)
    for (int y = 1; y <= n_bin_y; y++)
      for (int z = 1; z <= n_bin_z; z++)
        {
          int i_ge_bin = 0;
          for (n_list::const_iterator iter = src_list.begin();
              iter != src_list.end(); ++iter)
            {
              const int n = (*iter);

              string suffix;
              if (run_or_fill)
                suffix = get_suffix_run(n);
              else
                suffix = get_suffix_fill(n);

              sa_vec_hist::iterator histoiter = (_sahists).find(suffix);
              if (histoiter != (_sahists).end())
                {

                  saHist * sub_hist =
                      static_cast<saHist *>(histoiter->second.GetObject());

                  if (!sub_hist)
                    {
                      cout << "saHist::" << get_name()
                          << "::CalcAsymmetry_Chi2Fit - Error - can not load sub saHist "
                          << histoiter->first << endl;
                      continue;
                    }

                  // read the saHist here
                  TH1 *h_type = sub_hist->getHisto(hist_type.Data());

                  (ge_type->GetX())[i_ge_bin] = n;
                  (ge_type->GetY())[i_ge_bin] = h_type->GetBinContent(x, y, z);
                  (ge_type->GetEY())[i_ge_bin] = h_type->GetBinError(x, y, z);
                  i_ge_bin++;
                }
            }

          TFitResultPtr t_fit_result_ptr = ge_type->Fit("pol0", "MSQ");

          hFitPar->SetBinContent(x, y, z, t_fit_result_ptr->Parameter(0));
          hFitPar->SetBinError(x, y, z, t_fit_result_ptr->ParError(0));
          hChi2->SetBinContent(x, y, z, t_fit_result_ptr->Chi2());
          hNdf->SetBinContent(x, y, z, t_fit_result_ptr->Ndf());
        }

  if (ge_type)
    delete ge_type;

  return true;
}

//! use this saHist as template to make an unfilled new saHist
saHist *
saHist::MakeTemplate(std::string new_name)
{
  if (new_name.length() == 0)
    new_name = get_name();

  if (Verbosity())
    cout << "saHist::" << get_name()
        << "::MakeTemplate - INFO - make an unfilled saHist :" << new_name
        << endl;

  TH1 * new_base = static_cast<TH1 *>(get_clone_hist(new_name));
  sa_flags_t new_flags = get_flags();

  saHist * new_hist = new saHist(new_base, new_flags, Verbosity());

  for (sa_vec_hist::const_iterator hiter = _hists.begin();
      hiter != _hists.end(); ++hiter)
    {
      const string & name = (*hiter).first;

      new_hist->insert_clone_hist(name);
    }

  if (Verbosity())
    new_hist->Print();

  return new_hist;
}

//! Merge count in hist to this object
void
saHist::MergeHistCount(saHist * hist)
{

  if (!hist)
    {

      cout << "saHist::" << get_name()
          << "::MergeHistCount - Error - invalid input sub histogram " << endl;
      return;
    }

  //consistency check
  TH1 * dest = getHisto("BASE");
  TH1 * src = hist->getHisto("BASE");

  assert(dest);
  assert(src);

  if (dest->GetNbinsX() != src->GetNbinsX()
      or dest->GetNbinsY() != src->GetNbinsY()
      or dest->GetNbinsZ() != src->GetNbinsZ())
    {

      cout << "saHist::" << get_name()
          << "::MergeHistCount - Error - inconsistent input sub histogram dimensions:"
          << hist->GetName() << endl;
    }

  //merging
  for (sa_vec_hist::const_iterator hiter = (hist->_hists).begin();
      hiter != (hist->_hists).end(); ++hiter)
    {
      const string & suffix = (*hiter).first;
      const TRef & ref = (*hiter).second;
      const TH1 * h_src = dynamic_cast<const TH1 *>(ref.GetObject());

      if (!(ref.IsValid() && h_src))
        {

          cout << "saHist::" << get_name()
              << "::MergeHistCount - Error - invalid TH1 histogram :" << suffix
              << endl;

          continue;

        }

      insert_clone_hist(suffix);
      sa_vec_hist::const_iterator hiter_dets = _hists.find(suffix);
      TH1 * h_dest = dynamic_cast<TH1 *>(((*hiter_dets).second).GetObject());

      h_dest->Add(h_src);

    }

}

//! Add a saHist as sub histogram, merge its counts to the master histograms
//!  and possess the ownership to this object
void
saHist::MergeSubHist(const std::string & suffix, saHist * sub_hist)
{

  if (!sub_hist)
    {

      cout << "saHist::" << get_name()
          << "::MergeSubHist - Error - invalid input sub histogram " << endl;
      return;
    }

  saHist * dest_hist = getsaHist(suffix);
  if (!dest_hist)
    {

      if (Verbosity())
        cout << "saHist::" << get_name()
            << "::MergeSubHist - adopting and merge sub hist [" << suffix
            << "] = " << sub_hist->GetName() << endl;

      _sahists[suffix] = sub_hist;

      MergeHistCount(sub_hist);

    }
  else
    {

      if (!get_flag(LUMI_HIST))
        {

          if (Verbosity())
            cout << "saHist::" << get_name()
                << "::MergeSubHist - already exist sub hist [" << suffix
                << "]. Merge counts..." << endl;

          dest_hist->MergeHistCount(sub_hist);
          MergeHistCount(sub_hist);
        }
      else
        {
          if (Verbosity())
            cout << "saHist::" << get_name()
                << "::MergeSubHist - already exist luminorsity sub hist ["
                << suffix << "]. DO NOT Merge counts." << endl;
        }

      if (Verbosity())
        cout << "saHist::" << get_name()
            << "::MergeSubHist - done and delete the source."
            << sub_hist->GetName() << endl;
      delete sub_hist;
    }

}

//! adopt sub saHist of a source master saHist to this object.
//! Support use of a good run list. Those sub saHist will be removed from the source.
void
saHist::AdoptSubHist(saHist * master_hist, bool run_or_fill, bool use_good_list,
    const saHist::n_list & good_list)
{
  if (!master_hist)
    {

      cout << "saHist::" << get_name()
          << "::AdoptSubHist - Error - invalid input master histogram " << endl;
      return;
    }

  if (Verbosity())
    cout << "saHist::" << get_name()
        << "::AdoptSubHist - adopting sub hist from " << master_hist->GetName()
        << endl;

  n_list src_list;

  if (run_or_fill)
    src_list = master_hist->get_run_list();
  else
    src_list = master_hist->get_fill_list();

  if (src_list.size()==0)
    {
      cout << "saHist::" << get_name()
          << "::AdoptSubHist - WARNING - No sub hist found in " << master_hist->GetName()
          << endl;
    }

  for (n_list::const_iterator iter = src_list.begin(); iter != src_list.end();
      ++iter)
    {
      const int n = (*iter);

      if (use_good_list)
        {
          n_list::const_iterator iter_good = find(good_list.begin(),
              good_list.end(), n);
          if (iter_good == good_list.end())
            {
              if (Verbosity())
                cout << "saHist::" << get_name() << "::AdoptSubHist - " << n
                    << " is a bad fill/run. Ignore it" << endl;

              continue;
            }
        }

      string suffix;
      if (run_or_fill)
        suffix = master_hist->get_suffix_run(n);
      else
        suffix = master_hist->get_suffix_fill(n);

      sa_vec_hist::iterator histoiter = (master_hist->_sahists).find(suffix);
      if (histoiter != (master_hist->_sahists).end())
        {

          saHist * sub_hist =
              static_cast<saHist *>(histoiter->second.GetObject());

          (master_hist->_sahists).erase(histoiter);

          MergeSubHist(suffix, sub_hist);

        }
    }

}

//---------------------------------------------------------------------------------------------------------------
void
saHist::MergeRunSubHist(const std::string & suffix, saHist * sub_hist)
{

  if (!sub_hist)
    {

      cout << "saHist::" << get_name()
          << "::MergeRunSubHist - Error - invalid input sub histogram " << endl;
      return;
    }

  saHist * dest_hist = getsaHist(suffix);
  if (!dest_hist)
    {

      if (Verbosity())
        cout << "saHist::" << get_name()
            << "::MergeRunSubHist - adopting and merge sub hist [" << suffix
            << "] = " << sub_hist->GetName() << endl;

      _sahists[suffix] = sub_hist;
    }
  else
    {

      if (!get_flag(LUMI_HIST))
        {

          if (Verbosity())
            cout << "saHist::" << get_name()
                << "::MergeRunSubHist - already exist sub hist [" << suffix
                << "]. Merge counts..." << endl;

          dest_hist->MergeHistCount(sub_hist);
        }
      else
        {
          if (Verbosity())
            cout << "saHist::" << get_name()
                << "::MergeRunSubHist - already exist luminorsity sub hist ["
                << suffix << "]. DO NOT Merge counts." << endl;
        }

      if (Verbosity())
        cout << "saHist::" << get_name()
            << "::MergeRunSubHist - done and delete the source."
            << sub_hist->GetName() << endl;
      delete sub_hist;
    }

}

void
saHist::AdoptRunSubHist(saHist * master_hist, bool use_good_list,
    const saHist::n_list & good_list)
{
  if (!master_hist)
    {

      cout << "saHist::" << get_name()
          << "::AdoptRunSubHist - Error - invalid input master histogram "
          << endl;
      return;
    }

  if (Verbosity())
    cout << "saHist::" << get_name()
        << "::AdoptRunSubHist - adopting run sub hist from "
        << master_hist->GetName() << endl;

  n_list src_list = master_hist->get_run_list();

  for (n_list::const_iterator iter = src_list.begin(); iter != src_list.end();
      ++iter)
    {
      const int n = (*iter);

      if (use_good_list)
        {
          n_list::const_iterator iter_good = find(good_list.begin(),
              good_list.end(), n);
          if (iter_good == good_list.end())
            {
              if (Verbosity())
                cout << "saHist::" << get_name() << "::AdoptRunSubHist - " << n
                    << " is a bad run. Ignore it" << endl;

              continue;
            }
        }

      string suffix = master_hist->get_suffix_run(n);

      sa_vec_hist::iterator histoiter = (master_hist->_sahists).find(suffix);
      if (histoiter != (master_hist->_sahists).end())
        {

          saHist * sub_hist =
              static_cast<saHist *>(histoiter->second.GetObject());

          (master_hist->_sahists).erase(histoiter);

          MergeRunSubHist(suffix, sub_hist);

        }
    }

}

void
saHist::MergeFillSubHist(const std::string & suffix, saHist * sub_hist)
{

  if (!sub_hist)
    {

      cout << "saHist::" << get_name()
          << "::MergeFillSubHist - Error - invalid input sub histogram "
          << endl;
      return;
    }

  saHist * dest_hist = getsaHist(suffix);
  if (!dest_hist)
    {

      if (Verbosity())
        cout << "saHist::" << get_name()
            << "::MergeFillSubHist - adopting and merge fill sub hist ["
            << suffix << "] = " << sub_hist->GetName() << endl;

      _sahists[suffix] = sub_hist;
    }
  else
    {

      if (Verbosity())
        cout << "saHist::" << get_name()
            << "::MergeFillSubHist - already exist fill sub hist [" << suffix
            << "]. Merge counts..." << endl;

      dest_hist->MergeHistCount(sub_hist);

      if (Verbosity())
        cout << "saHist::" << get_name()
            << "::MergeFillSubHist - done and delete the source."
            << sub_hist->GetName() << endl;
      delete sub_hist;
    }

}

void
saHist::AdoptFillSubHist(saHist * master_hist, bool use_good_list,
    const saHist::n_list & good_list)
{
  if (!master_hist)
    {

      cout << "saHist::" << get_name()
          << "::AdoptFillSubHist - Error - invalid input master histogram "
          << endl;
      return;
    }

  if (Verbosity())
    cout << "saHist::" << get_name()
        << "::AdoptFillSubHist - adopting fill sub hist from "
        << master_hist->GetName() << endl;

  n_list src_list = master_hist->get_fill_list();

  for (n_list::const_iterator iter = src_list.begin(); iter != src_list.end();
      ++iter)
    {
      const int n = (*iter);

      if (use_good_list)
        {
          n_list::const_iterator iter_good = find(good_list.begin(),
              good_list.end(), n);
          if (iter_good == good_list.end())
            {
              if (Verbosity())
                cout << "saHist::" << get_name() << "::AdoptFillSubHist - " << n
                    << " is a bad fill. Ignore it" << endl;

              continue;
            }
        }

      string suffix = master_hist->get_suffix_fill(n);

      sa_vec_hist::iterator histoiter = (master_hist->_sahists).find(suffix);
      if (histoiter != (master_hist->_sahists).end())
        {

          saHist * sub_hist =
              static_cast<saHist *>(histoiter->second.GetObject());

          (master_hist->_sahists).erase(histoiter);

          MergeFillSubHist(suffix, sub_hist);

        }
    }

}

//! remove a sahist from the list
bool
saHist::removeFillsaHist(const std::string& suffix)
{
  if (get_flag(SUB_HIST))
    {
      if (Verbosity() > 1)
        {
          cout << "saHist::" << _name
              << "::removeFillsaHist ERROR - it is a SUB_HIST!" << endl;
        }
      return false;
    }
  else if (!get_flag(FILL_PROPERTY))
    {
      if (Verbosity() > 1)
        {
          cout << "saHist::" << _name
              << "::removeFillsaHist ERROR - no FILL_PROPERTY flag!" << endl;
        }
      return false;
    }

  sa_vec_hist::iterator iter = _sahists.find(suffix);

  if (iter != _sahists.end())
    {
      if (Verbosity() > 1)
        {
          cout << "saHist::" << _name << "::removeFillsaHist: erase saHist "
              << suffix << endl;
        }

      //      if (ref.IsValid() && ref.GetObject())
      //        delete ref.GetObject();

      _sahists.erase(iter);

      if (Verbosity() > 2)
        {
          cout << "saHist::" << _name
              << "::removeFillsaHist: remaining saHists: ";
          for (sa_vec_hist::const_iterator hiter = _sahists.begin();
              hiter != _sahists.end(); ++hiter)
            {
              const string name = (*hiter).first;

              cout << name << ", ";

            }
          cout << endl;
        }

      return true;

    }

  if (Verbosity() > 1)
    {
      cout << "saHist::" << _name
          << "::removeFillsaHist: ERROR Unknown Fill saHist " << suffix << endl;
    }
  return false;

}

bool
saHist::removeRunsaHist(const std::string& suffix)
{
  if (get_flag(SUB_HIST))
    {
      if (Verbosity() > 1)
        {
          cout << "saHist::" << _name
              << "::removeRunsaHist ERROR - it is a SUB_HIST!" << endl;
        }
      return false;
    }
  else if (!get_flag(RUN_PROPERTY))
    {
      if (Verbosity() > 1)
        {
          cout << "saHist::" << _name
              << "::removeRunsaHist ERROR - no RUN_PROPERTY flag!" << endl;
        }
      return false;
    }

  sa_vec_hist::iterator iter = _sahists.find(suffix);

  if (iter != _sahists.end())
    {
      if (Verbosity() > 1)
        {
          cout << "saHist::" << _name << "::removeRunsaHist: erase saHist "
              << suffix << endl;
        }

      //      if (ref.IsValid() && ref.GetObject())
      //        delete ref.GetObject();

      _sahists.erase(iter);

      if (Verbosity() > 2)
        {
          cout << "saHist::" << _name << "::removeRunsaHist: remaining saHists: ";
          for (sa_vec_hist::const_iterator hiter = _sahists.begin();
              hiter != _sahists.end(); ++hiter)
            {
              const string name = (*hiter).first;

              cout << name << ", ";

            }
          cout << endl;
        }

      return true;      

    }

  if (Verbosity() > 1)
    {
      cout << "saHist::" << _name << "::removeRunsaHist: ERROR Unknown Run saHist "
          << suffix << endl;
    }
  return false;

}
