// $Id: saModuleBase.h,v 1.1 2013/05/10 15:59:47 jinhuang Exp $                                                                                             

/*!
 * \file saModuleBase.h
 * \brief 
 * \author Jin Huang <jhuang@bnl.gov>
 * \version $Revision: 1.1 $
 * \date $Date: 2013/05/10 15:59:47 $
 */

#ifndef SAMODULEBASE_H_
#define SAMODULEBASE_H_

#include <SubsysReco.h>
#include <Fun4AllBase.h>
#include <string>

#include <Fun4AllReturnCodes.h>

#include "saHistManager.h"


/*!
 * \brief saModuleBase
 */
class saModuleBase : public Fun4AllBase
{
public:
  saModuleBase(const std::string &name);

  virtual
  ~saModuleBase();

#ifndef __CINT__

  //! global initialization
  virtual int
  init(PHCompositeNode *topNode, sa_hist_mangager_ptr hm) {return EVENT_OK;}

  //! Run initialization
  virtual int
  init_run(PHCompositeNode *topNode, sa_hist_mangager_ptr hm)  {return EVENT_OK;}

  //! event method
  virtual int
  event(PHCompositeNode *topNode, sa_hist_mangager_ptr hm)  {return EVENT_OK;}

  //! global termination
  virtual int
  end(PHCompositeNode *topNode, sa_hist_mangager_ptr hm)  {return EVENT_OK;}

#endif

private:

};

#endif /* SAMODULEBASE_H_ */
