// $Id: saFlag.h,v 1.1 2013/07/23 04:40:09 jinhuang Exp $

/*!
 * \file saFlag.h
 * \brief 
 * \author Jin Huang <jhuang@bnl.gov>
 * \version $Revision: 1.1 $
 * \date $Date: 2013/07/23 04:40:09 $
 */

#ifndef saFlag_H_
#define saFlag_H_

#include "PHObject.h"

/*!
 * \brief saFlag save an array of flag which is synchronized with picodst_objects
 */

template<class TArray_t>
  class saFlag : public PHObject
  {
  public:
    saFlag()
    {
    }

    virtual
    ~saFlag()
    {
    }

    /** identify Function from PHObject
     @param os Output Stream
     */
    virtual void
    identify(std::ostream& os = std::cout) const
    {
      std::cout << "This is " << this->ClassName() << std::endl;
    }

    //! Clear Event
    virtual void
    Reset(){flags.Reset(0);}

  public:

    TArray_t flags;

  ClassDef(saFlag, 1)

  };

templateClassImp(saFlag)

#include <TArrayI.h>

typedef saFlag<TArrayI> saFlagI;

#include <TArrayS.h>

typedef saFlag<TArrayS> saFlagS;

#include <TArrayC.h>

typedef saFlag<TArrayC> saFlagC;

#endif /* saFlag_H_ */
