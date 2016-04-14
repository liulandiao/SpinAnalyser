// $Id: saHistManager.h,v 1.3 2013/05/12 01:33:45 jinhuang Exp $                                                                                             

/*!
 * \file saHistManager.h
 * \brief 
 * \author Jin Huang <jhuang@bnl.gov>
 * \version $Revision: 1.3 $
 * \date $Date: 2013/05/12 01:33:45 $
 */

#ifndef SAHISTMANAGER_H_
#define SAHISTMANAGER_H_

#ifndef __CINT__
#include <boost/smart_ptr.hpp>
#endif

#include <map>
#include <string>

#include <Fun4AllBase.h>

#include "saHist.h"

class Fun4AllHistoManager;
/*!
 * \brief saHistManager
 */
class saHistManager : public Fun4AllBase
{
public:
  saHistManager();
  virtual
  ~saHistManager();

  //! what do I got?
  void
  Print(const std::string &what = "ALL") const;

  //! add a histogram
  bool
  registerHisto(saHist * h1d);

  //! is some histogram already there?
  int
  isHistoRegistered(const std::string &name) const;

  //! get the histogram of default relative luminorsity
  saHist *
  getHisto_DefaultLumi() const;

  //! register the histogram of default relative luminorsity
  bool
  registerHisto_DefaultLumi(saHist * h1d);

  saHist *
  getHisto(const std::string &hname) const;

  saHist *
  getHisto(const unsigned int ihisto) const;

  const char *
  getHistoName(const unsigned int ihisto) const;

  unsigned int
  nHistos() const;

  void
  Reset();

  static Fun4AllHistoManager *
  getFun4AllHistoManager();

  int
  dumpHistos(const std::string &filename = "", const std::string &openmode =
      "RECREATE");

protected:

#ifndef __CINT__

//  typedef boost::shared_ptr<saHist> sa_hist_ptr;

#endif

  typedef saHist * sa_hist_ptr;
  typedef std::map<const std::string, sa_hist_ptr> sa_vec_hist;

  sa_vec_hist _hists;

  std::string _default_lumi_name;

};

#ifndef __CINT__
typedef boost::shared_ptr<saHistManager> sa_hist_mangager_ptr;
#else
typedef saHistManager * sa_hist_mangager_ptr;
#endif

typedef saHistManager & sa_hist_mangager_ref;

#endif /* SAHISTMANAGER_H_ */
