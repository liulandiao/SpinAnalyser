// $Id: saModSimpleHist.h,v 1.5 2016/02/08 02:19:28 jinhuang Exp $                                                                                             

/*!
 * \file saModSimpleHist.h
 * \brief 
 * \author Jin Huang <jhuang@bnl.gov>
 * \version $Revision: 1.5 $
 * \date $Date: 2016/02/08 02:19:28 $
 */

#ifndef SAMODSIMPLEHIST_H_
#define SAMODSIMPLEHIST_H_

#include <string>
#include <map>
#include <utility>

#include "saModuleBase.h"

class SingleMuonContainer;
class DiMuonContainer;

/*!
 * \brief example module to use spin in-dependent histograms
 */
class saModSimpleHist : public saModuleBase
{
  // --------------------------------------------------------------------------
  //!@name Interface
  //@{
public:
  saModSimpleHist(const std::string &name = "SimpleHist");
  virtual
  ~saModSimpleHist();

#ifndef __CINT__
  //! global initialization
  virtual int
  init(PHCompositeNode* topNode, sa_hist_mangager_ptr hm);
  //! Run initialization
  virtual int
  init_run(PHCompositeNode* topNode, sa_hist_mangager_ptr hm);
  //! event method
  virtual int
  event(PHCompositeNode* topNode, sa_hist_mangager_ptr hm);
  //! global termination
  virtual int
  end(PHCompositeNode* topNode, sa_hist_mangager_ptr hm);


#endif

  void set_usv_svx_only(int i) {_usv_svx_only = i;}
private:
  int _run_num;
  int _usv_svx_only;

  //@}

  // --------------------------------------------------------------------------
  //!@name Vertex Fit
  //@{

public:
  saHist*
  get_vtx() const
  {
    return _h_vtx;
  }

  void
  set_vtx(saHist* vtx)
  {
    _h_vtx = vtx;
  }

  std::string
  get_vertex_name() const
  {
    return _vertex_name;
  }

  void
  set_vertex_name(std::string vertexName)
  {
    _vertex_name = vertexName;
  }

private:

  std::string _vertex_name;

  saHist * _h_vtx;

  saHist * _h_vtx_xz;
  saHist * _h_vtx_yz;

  //@}

  // --------------------------------------------------------------------------
  //!@name DCA
  //@{

public:

  //! call this function to load beam position data
  //! in the order of file >> run >> x >> y >> ex >> ey;
  void
  load_beam_xy_data(std::string data_file = "BeamPos.dat")
  {
    _beam_pos_xy.load_data(data_file, Verbosity());
  }

  //! get precise beam X-Y position from off-line fit of VTX vertex
  class BeamPosXY
  {
  public:
    //! load data files
    void
    load_data(std::string data_file = "BeamPos.dat", int verbose = 0);

    //! have this run in data file?
    bool
    have_run(int run);

    //! get beam position in x
    double
    get_x(int run);
    //! get beam position in x
    double
    get_y(int run);

    //! get beam position in y
    double
    get_dx(int run);

    //! get beam position in y
    double
    get_dy(int run);

    //! position in x y
    typedef std::pair<double,double> pos;
    typedef std::map<int, pos> record;

    //! map to run -> pos
    record position;
    //! map to run -> weight = 1/error^2
    record weight;
    //! map to run -> slope
    record dposition;
    //! map to run -> weight = 1/error^2
    record dweight;

  };

private:

  void
  DCA(SingleMuonContainer * smc, DiMuonContainer * dmc);

  void
  DCA_FVTX(SingleMuonContainer * smc, DiMuonContainer * dmc);

  saHist * _h_DCA_Phi_Arm;
  saHist * _h_DCA_Pz_Arm;

  saHist * _h_DCA_Phi_z_Arm0;
  saHist * _h_DCA_Phi_z_Arm1;

  saHist * _h_fvtx_DCA_Phi_Arm;
  saHist * _h_fvtx_DCA_Pz_Arm;
  saHist * _h_Normalization;

  BeamPosXY _beam_pos_xy;

  double _VTX_Avg_X;
  double _VTX_Avg_Y;
  double _VTX_Avg_dX;
  double _VTX_Avg_dY;

  std::vector<double> _event_vertex;

  //@}

};

#endif /* SAMODSIMPLEHIST_H_ */
