// $Id: saModuleJPsi.h,v 1.2 2013/11/05 16:23:32 keyaaron Exp $                                                                                             

/*!
 * \file saModuleJPsi.h
 * \brief 
 * \author Aaron Key <keyaaron@unm.edu>
 * \version $Revision: 1.2 $
 * \date $Date: 2013/11/05 16:23:32 $
 */

#ifndef SAMODULEJPSI_H_
#define SAMODULEJPSI_H_

#include <vector>

#include "saModuleBase.h"

class DiMuon;
class saEventProperty;
class saHist;
class TH1;

/*!
 * \brief saModuleJPsi
 */
class saModuleJPsi : public saModuleBase
{
public:
  saModuleJPsi(const std::string &name = "JPsi");
  virtual
  ~saModuleJPsi();

#ifndef __CINT__

  //! global initialization
  virtual int
  init(PHCompositeNode *topNode, sa_hist_mangager_ptr hm);

  //! Run initialization
  virtual int
  init_run(PHCompositeNode *topNode, sa_hist_mangager_ptr hm);

  //! event method
  virtual int
  event(PHCompositeNode *topNode, sa_hist_mangager_ptr hm);

  //! global termination
  virtual int
  end(PHCompositeNode *topNode, sa_hist_mangager_ptr hm);

#endif
  saHist *_h_mass_pT,*_h_pT_like;
  saHist *_h_pT,*_h_pT_m,*_h_pT_p;
  saHist *_h_pT_side,*_h_pT_wide;

  saHist *_h_mass_pT_FVTX,*_h_pT_like_FVTX;
  saHist *_h_pT_FVTX,*_h_pT_m_FVTX,*_h_pT_p_FVTX;
  saHist *_h_pT_side_FVTX,*_h_pT_wide_FVTX;

  saHist *_h_mass_pT_sel,*_h_pT_like_sel;
  saHist *_h_pT_sel,*_h_pT_m_sel,*_h_pT_p_sel;
  saHist *_h_pT_side_sel,*_h_pT_wide_sel;

private:
  saHist* make_saHist(sa_hist_mangager_ptr hm, TH1* h);

  bool passesCuts() const;
  bool passesFVTXCuts() const;

  void selectOutputHistos(bool has_fvtx);
  void fillHistos();

  std::string dimuonType() const;

  DiMuon* _dimuon;
  const saEventProperty* _ep;
  std::vector<saHist*> _histos;
};

#endif /* SAMODULEJPSI_H_ */
