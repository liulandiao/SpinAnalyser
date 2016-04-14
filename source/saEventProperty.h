// $Id: saEventProperty.h,v 1.2 2013/05/12 01:33:45 jinhuang Exp $                                                                                             

/*!
 * \file saEventProperty.h
 * \brief 
 * \author Jin Huang <jhuang@bnl.gov>
 * \version $Revision: 1.2 $
 * \date $Date: 2013/05/12 01:33:45 $
 */

#ifndef SAEVENTPROPERTY_H_
#define SAEVENTPROPERTY_H_

#include "PHObject.h"

/*!
 * \brief saEventProperty
 */
class saEventProperty : public PHObject
{
public:
  saEventProperty();
  virtual
  ~saEventProperty();

  /** identify Function from PHObject
   @param os Output Stream
   */
  virtual void
  identify(std::ostream& os = std::cout) const
  {
    std::cout << "This is " << this->ClassName() << " with run = "
        << get_run_number() //
        << " spin state = " << (int) get_spin_state() << std::endl;
  }

public:

  //! spin states , in order of Blue and Yellow spin sign
  enum sa_spin_state
  {
    PP, PM, MP, MM, UNKNOWN, SUM
  };

  //! number of RHIC crossings
  static const int N_CROSSING = 120;

  static const std::string
  get_state_name(sa_spin_state state)
  {
    switch (state)
      {
    case PP:
      return std::string("PP");
    case PM:
      return std::string("PM");
    case MP:
      return std::string("MP");
    case MM:
      return std::string("MM");
    case SUM:
      return std::string("SUM");
    default:
      return std::string("UNKNOWN");
      }
    return std::string("UNKNOWN");
  }

  virtual short
  get_crossing_id() const
  {
    return 0;
  }

  virtual void
  set_crossing_id(short crossingId)
  {
  }

  virtual
  short
  get_crossing_id_RHIC() const
  {
    return 0;
  }

  virtual
  void
  set_crossing_id_RHIC(short crossingIdRhic)
  {
  }

  virtual int
  get_event_number() const
  {
    return 0;
  }

  virtual void
  set_event_number(int eventNumber)
  {
  }

  virtual int
  get_fill_number() const
  {
    return 0;
  }

  virtual void
  set_fill_number(int fillNumber)
  {
  }

  virtual float
  get_polarization_blue() const
  {
    return 0;
  }

  virtual void
  set_polarization_blue(float polarization)
  {
  }

  virtual float
  get_polarization_yellow() const
  {
    return 0;
  }

  virtual void
  set_polarization_yellow(float polarization)
  {
  }

  virtual
  float
  get_polarization_stat_err_blue() const
  {
    return 0;
  }

  virtual
  void
  set_polarization_stat_err_blue(float polarizationStatErrBlue)
  {
  }

  virtual
  float
  get_polarization_stat_err_yellow() const
  {
    return 0;
  }

  virtual
  void
  set_polarization_stat_err_yellow(float polarizationStatErrYellow)
  {
  }

  virtual
  float
  get_polarization_sys_err_blue() const
  {
    return 0;
  }

  virtual
  void
  set_polarization_sys_err_blue(float polarizationSysErrBlue)
  {
  }

  virtual
  float
  get_polarization_sys_err_yellow() const
  {
    return 0;
  }

  virtual
  void
  set_polarization_sys_err_yellow(float polarizationSysErrYellow)
  {
  }

  virtual int
  get_run_number() const
  {
    return 0;
  }

  virtual
  void
  set_run_number(int runNumber)
  {
  }

  virtual sa_spin_state
  get_spin_state() const
  {
    return UNKNOWN;
  }

  virtual
  void
  set_spin_state(sa_spin_state spinState)
  {
  }

protected:

ClassDef(saEventProperty, 1)

};

#endif /* SAEVENTPROPERTY_H_ */
