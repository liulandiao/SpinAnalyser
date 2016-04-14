// $Id: saModuleDimuonMing.h,v 1.2 2013/08/17 05:42:48 mxliu Exp $                                                                                             
// 
// modefied from Jin's saModuleSimpleDimuon example   7/12/2013  MXL
//
// for Run12+ dimuons analyses: J/Psi, DY and B2B dimuons
//
//////////////////////////////////////////////////////////////////////

#ifndef SAMODULEDIMUONMING_H_
#define SAMODULEDIMUONMING_H_


#include "saModuleBase.h"

/*!
 * \brief saModuleSimpleDimuon
 */
class saModuleDimuonMing : public saModuleBase
{
public:
  saModuleDimuonMing(const std::string &name = "DiMuonMing");
  virtual
  ~saModuleDimuonMing();

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

  // -- Same Arm dimuons
  saHist * _h_mass;
  saHist * _h_mass_n;
  saHist * _h_mass_s;

  saHist * _h_mass_p;    // ++ charged dimuons
  saHist * _h_mass_p_n;  // ++ charged dimuons, north
  saHist * _h_mass_p_s;  // ++ charged dimuons, south

  saHist * _h_mass_m;    // -- charged dimuons
  saHist * _h_mass_m_n;  // -- charged dimuons, north
  saHist * _h_mass_m_s;  // -- charged dimuons, south


  saHist * _h_pT;
  saHist * _h_pT_n;      // north arm
  saHist * _h_pT_s;      // south arm

  saHist * _h_pT_p;
  saHist * _h_pT_p_n;      // north arm
  saHist * _h_pT_p_s;      // south arm

  saHist * _h_pT_m;
  saHist * _h_pT_m_n;      // north arm
  saHist * _h_pT_m_s;      // south arm


  // -- B2B dimuons
  saHist * _h_b2b;
  saHist * _h_b2b_p;  // ++ charged dimuons
  saHist * _h_b2b_m;  // -- charged dimuons


};

#endif /* SAMODULEDIMUONMING_H_ */
