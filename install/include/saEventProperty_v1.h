// $Id: saEventProperty_v1.h,v 1.2 2013/05/12 01:33:45 jinhuang Exp $                                                                                             

/*!
 * \file saEventProperty_v1.h
 * \brief 
 * \author Jin Huang <jhuang@bnl.gov>
 * \version $Revision: 1.2 $
 * \date $Date: 2013/05/12 01:33:45 $
 */

#ifndef saEventProperty_v1_H_
#define saEventProperty_v1_H_

#include "saEventProperty.h"

/*!
 * \brief saEventProperty_v1
 */
class saEventProperty_v1 : public saEventProperty
{
public:
  saEventProperty_v1();
  virtual
  ~saEventProperty_v1();

  //! isValid returns non zero if object contains vailid data
  virtual int
  isValid() const
  {
    return 1;
  }

  //! Clear Event
  virtual void
  Reset();

public:

  short
  get_crossing_id() const
  {
    return _crossing_id;
  }

  void
  set_crossing_id(short crossingId)
  {
    _crossing_id = crossingId;
  }

  short
  get_crossing_id_RHIC() const
  {
    return _crossing_id_RHIC;
  }

  void
  set_crossing_id_RHIC(short crossingIdRhic)
  {
    _crossing_id_RHIC = crossingIdRhic;
  }

  int
  get_event_number() const
  {
    return _event_number;
  }

  void
  set_event_number(int eventNumber)
  {
    _event_number = eventNumber;
  }

  int
  get_fill_number() const
  {
    return _fill_number;
  }

  void
  set_fill_number(int fillNumber)
  {
    _fill_number = fillNumber;
  }

  float
  get_polarization_blue() const
  {
    return _polarization_blue;
  }

  void
  set_polarization_blue(float polarization)
  {
    _polarization_blue = polarization;
  }

  float
  get_polarization_yellow() const
  {
    return _polarization_yellow;
  }

  void
  set_polarization_yellow(float polarization)
  {
    _polarization_yellow = polarization;
  }

  float
  get_polarization_stat_err_blue() const
  {
    return _polarization_stat_err_blue;
  }

  void
  set_polarization_stat_err_blue(float polarizationStatErrBlue)
  {
    _polarization_stat_err_blue = polarizationStatErrBlue;
  }

  float
  get_polarization_stat_err_yellow() const
  {
    return _polarization_stat_err_yellow;
  }

  void
  set_polarization_stat_err_yellow(float polarizationStatErrYellow)
  {
    _polarization_stat_err_yellow = polarizationStatErrYellow;
  }

  float
  get_polarization_sys_err_blue() const
  {
    return _polarization_sys_err_blue;
  }

  void
  set_polarization_sys_err_blue(float polarizationSysErrBlue)
  {
    _polarization_sys_err_blue = polarizationSysErrBlue;
  }

  float
  get_polarization_sys_err_yellow() const
  {
    return _polarization_sys_err_yellow;
  }

  void
  set_polarization_sys_err_yellow(float polarizationSysErrYellow)
  {
    _polarization_sys_err_yellow = polarizationSysErrYellow;
  }

  int
  get_run_number() const
  {
    return _run_number;
  }

  void
  set_run_number(int runNumber)
  {
    _run_number = runNumber;
  }

  sa_spin_state
  get_spin_state() const
  {
    return _spin_state;
  }

  void
  set_spin_state(sa_spin_state spinState)
  {
    _spin_state = spinState;
  }

protected:

  int _fill_number;

  int _run_number;

  int _event_number;

  //! PHENIX crossing id = TrigLvl1::get_lvl1_clock_cross()
  short _crossing_id;

  //! RHIC crossing id = (crossing + SpinDBContent::cross_shift)%120
  short _crossing_id_RHIC;

  sa_spin_state _spin_state;

  float _polarization_blue;

  float _polarization_yellow;

  float _polarization_stat_err_blue;

  float _polarization_stat_err_yellow;

  float _polarization_sys_err_blue;

  float _polarization_sys_err_yellow;

ClassDef(saEventProperty_v1, 1)

};

#endif /* saEventProperty_v1_H_ */
