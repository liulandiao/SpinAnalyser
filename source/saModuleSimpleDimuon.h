// $Id: saModuleSimpleDimuon.h,v 1.1 2013/05/10 15:59:48 jinhuang Exp $                                                                                             

/*!
 * \file saModuleSimpleDimuon.h
 * \brief 
 * \author Jin Huang <jhuang@bnl.gov>
 * \version $Revision: 1.1 $
 * \date $Date: 2013/05/10 15:59:48 $
 */

#ifndef SAMODULESIMPLEDIMUON_H_
#define SAMODULESIMPLEDIMUON_H_


#include "saModuleBase.h"

/*!
 * \brief saModuleSimpleDimuon
 */
class saModuleSimpleDimuon : public saModuleBase
{
public:
  saModuleSimpleDimuon(const std::string &name = "SimpleDiMuon");
  virtual
  ~saModuleSimpleDimuon();

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


  saHist * _h_mass;

};

#endif /* SAMODULESIMPLEDIMUON_H_ */
