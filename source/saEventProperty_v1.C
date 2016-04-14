// $Id: saEventProperty_v1.C,v 1.2 2013/05/12 05:48:22 jinhuang Exp $                                                                                             

/*!
 * \file saEventProperty_v1.C
 * \brief 
 * \author Jin Huang <jhuang@bnl.gov>
 * \version $Revision: 1.2 $
 * \date $Date: 2013/05/12 05:48:22 $
 */

#include "saEventProperty_v1.h"

ClassImp(saEventProperty_v1);

saEventProperty_v1::saEventProperty_v1()
{
  Reset();
}

saEventProperty_v1::~saEventProperty_v1()
{
  // TODO Auto-generated destructor stub
}

//! Clear Event
void
saEventProperty_v1::Reset()
{

  _fill_number = -1;
  _run_number = -1;
  _event_number = -1;

  _crossing_id = -1;
  _crossing_id_RHIC = -1;
  _spin_state = UNKNOWN;

  _polarization_blue = -1;
  _polarization_yellow = -1;

  _polarization_stat_err_blue = 1;

  _polarization_stat_err_yellow = 1;

  _polarization_sys_err_blue = 1;

  _polarization_sys_err_yellow = 1;

}
