// $Id: saModSimpleHist.C,v 1.6 2016/02/08 02:19:27 jinhuang Exp $                                                                                             

/*!
 * \file saModSimpleHist.C
 * \brief 
 * \author Jin Huang <jhuang@bnl.gov>
 * \version $Revision: 1.6 $
 * \date $Date: 2016/02/08 02:19:27 $
 */

#include <TH2F.h>
#include <TH3F.h>
#include <TMath.h>

#include <EventHeader.h>
#include <RunHeader.h>
#include <SyncObject.h>
#include <TMutNode.h>
#include <recoConsts.h>
#include <VtxOut.h>
#include <PHPoint.h>
#include <getClass.h>

#include <cmath>
#include <iostream>
#include <fstream>
#include <boost/algorithm/string.hpp>

#include <DiMuonContainer.h>
#include <SingleMuonContainer.h>

#include "saModSimpleHist.h"

using namespace std;

saModSimpleHist::saModSimpleHist(const std::string &name) :
    saModuleBase(name), _vertex_name("SVX_PRECISE"), //
    _h_vtx(NULL), //
    _h_DCA_Phi_Arm(NULL), //
    _h_DCA_Pz_Arm(NULL), //
    _h_fvtx_DCA_Phi_Arm(NULL), //
    _h_fvtx_DCA_Pz_Arm(NULL), _h_Normalization(NULL), //
    _event_vertex(3) //
{
  // TODO Auto-generated constructor stub

  _usv_svx_only = 0;

  _VTX_Avg_X = 0;
  _VTX_Avg_Y = 0;
  _VTX_Avg_dX = 0;
  _VTX_Avg_dY = 0;
}

saModSimpleHist::~saModSimpleHist()
{
  // TODO Auto-generated destructor stub
}

//! global initialization
int
saModSimpleHist::init(PHCompositeNode *topNode, sa_hist_mangager_ptr hm)
{

  TH2F * h_vtx = new TH2F("VTX_XY",
      "Vertex Transverse Distribution;VTX X (cm);VTX Y (cm)", //
      300, -.3, .3, 300, -.3, .3);
  _h_vtx = new saHist(h_vtx, //the template histogram, any TH1 and derivatives is accepted
      saHist::EVENT_PROPERTY | saHist::RUN_PROPERTY | saHist::FILL_PROPERTY, // flags
      Verbosity() // verbosity
      );
  hm.get()->registerHisto(_h_vtx);

  TH2F * h_vtx_xz = new TH2F("VTX_XZ",
      "Vertex Distribution;VTX Z (cm);VTX X (cm)", //
      300, -15, 15, 300, -.3, .3);
  _h_vtx_xz = new saHist(h_vtx_xz, //the template histogram, any TH1 and derivatives is accepted
      saHist::EVENT_PROPERTY | saHist::RUN_PROPERTY | saHist::FILL_PROPERTY, // flags
      Verbosity() // verbosity
      );
  hm.get()->registerHisto(_h_vtx_xz);

  TH2F * h_vtx_yz = new TH2F("VTX_YZ",
      "Vertex Distribution;VTX Z (cm);VTX Y (cm)", //
      300, -15, 15, 300, -.3, .3);
  _h_vtx_yz = new saHist(h_vtx_yz, //the template histogram, any TH1 and derivatives is accepted
      saHist::EVENT_PROPERTY | saHist::RUN_PROPERTY | saHist::FILL_PROPERTY, // flags
      Verbosity() // verbosity
      );
  hm.get()->registerHisto(_h_vtx_yz);

  static const double DCA_lim = 1;
  static const int DCA_bin = 500;
  static const double Phi_lim = TMath::Pi();
  static const int Phi_bin = 64;

  //
  TH3F * DCA_Phi_z_Arm0 = new TH3F("DCA_Phi_z_Arm0",
      "DCA VS Phi;DCA (cm);Phi (rad);z (cm)", //
      DCA_bin, -DCA_lim, DCA_lim, //
      Phi_bin, -Phi_lim, Phi_lim, //
      4, -10, 10);
  _h_DCA_Phi_z_Arm0 = new saHist(DCA_Phi_z_Arm0, //the template histogram, any TH1 and derivatives is accepted
      saHist::EVENT_PROPERTY | saHist::FILL_PROPERTY, // flags
      Verbosity() // verbosity
      );
  hm.get()->registerHisto(_h_DCA_Phi_z_Arm0);

  //
  TH3F * DCA_Phi_z_Arm1 = new TH3F("DCA_Phi_z_Arm1",
      "DCA VS Phi;DCA (cm);Phi (rad);z (cm)", //
      DCA_bin, -DCA_lim, DCA_lim, //
      Phi_bin, -Phi_lim, Phi_lim, //
      4, -10, 10);
  _h_DCA_Phi_z_Arm1 = new saHist(DCA_Phi_z_Arm1, //the template histogram, any TH1 and derivatives is accepted
      saHist::EVENT_PROPERTY | saHist::FILL_PROPERTY, // flags
      Verbosity() // verbosity
      );
  hm.get()->registerHisto(_h_DCA_Phi_z_Arm1);

  //
  TH3F * h_DCA_Pz_Arm = new TH3F("DCA_Pz_Arm",
      "DCA VS Phi;DCA (cm);Pz (GeV/c);Arm", //
      DCA_bin, -DCA_lim, DCA_lim, //
      40, 2, 100, //
      4, -.5, 3.5);
  _h_DCA_Pz_Arm = new saHist(h_DCA_Pz_Arm, //the template histogram, any TH1 and derivatives is accepted
      saHist::EVENT_PROPERTY | saHist::FILL_PROPERTY, // flags
      Verbosity() // verbosity
      );
  hm.get()->registerHisto(_h_DCA_Pz_Arm);

  //
  TH3F * h_DCA_Phi_Arm = new TH3F("DCA_Phi_Arm",
      "DCA VS Phi;DCA (cm);Phi (rad);Arm", //
      DCA_bin, -DCA_lim, DCA_lim, //
      Phi_bin, -Phi_lim, Phi_lim, //
      4, -.5, 3.5);
  _h_DCA_Phi_Arm = new saHist(h_DCA_Phi_Arm, //the template histogram, any TH1 and derivatives is accepted
      saHist::EVENT_PROPERTY | saHist::FILL_PROPERTY, // flags
      Verbosity() // verbosity
      );
  hm.get()->registerHisto(_h_DCA_Phi_Arm);
  //
  TH3F * h_FVTX_DCA_Phi_Arm = new TH3F("FVTX_DCA_Phi_Arm",
      "DCA VS Phi;DCA (cm);Phi (rad);Arm", //
      DCA_bin, -DCA_lim, DCA_lim, //
      Phi_bin, -Phi_lim, Phi_lim, //
      4, -.5, 3.5);
  _h_fvtx_DCA_Phi_Arm = new saHist(h_FVTX_DCA_Phi_Arm, //the template histogram, any TH1 and derivatives is accepted
      saHist::EVENT_PROPERTY | saHist::FILL_PROPERTY, // flags
      Verbosity() // verbosity
      );
  hm.get()->registerHisto(_h_fvtx_DCA_Phi_Arm);
  //
  h_DCA_Pz_Arm = new TH3F("FVTX_DCA_Pz_Arm",
      "DCA VS Phi;DCA (cm);Pz (GeV/c);Arm", //
      DCA_bin, -DCA_lim, DCA_lim, //
      40, 2, 100, //
      4, -.5, 3.5);
  _h_fvtx_DCA_Pz_Arm = new saHist(h_DCA_Pz_Arm, //the template histogram, any TH1 and derivatives is accepted
      saHist::EVENT_PROPERTY | saHist::FILL_PROPERTY, // flags
      Verbosity() // verbosity
      );
  hm.get()->registerHisto(_h_fvtx_DCA_Pz_Arm);

  //
  TH1D * Normalization = new TH1D("Normalization", "Normalization;Item", 10,
      0.5, 10.5);
  Normalization->GetXaxis()->SetBinLabel(1, "Event");
  Normalization->GetXaxis()->SetBinLabel(2, "Vertex Cut");
  Normalization->GetXaxis()->SetBinLabel(3, "Single Muon");
  Normalization->GetXaxis()->SetBinLabel(4, "Single Muon J/Psi");
  Normalization->GetXaxis()->SetBinLabel(5, "Single Muon J/Psi good_muID_muon");
  Normalization->GetXaxis()->SetBinLabel(6, "Single Muon pz4 good_muID_muon");
//  Normalization->GetXaxis()->SetBinLabel(7, "Dimuon");
  Normalization->GetXaxis()->LabelsOption("v");
  _h_Normalization = new saHist(Normalization, //the template histogram, any TH1 and derivatives is accepted
      saHist::EVENT_PROPERTY | saHist::RUN_PROPERTY | saHist::FILL_PROPERTY, // flags
      Verbosity() // verbosity
      );

  hm.get()->registerHisto(_h_Normalization);

  if (Verbosity())
    cout << "saModSimpleHist::init - N beam position records = "
        << _beam_pos_xy.position.size() << endl;

  return EVENT_OK;
}

//! Run initialization
int
saModSimpleHist::init_run(PHCompositeNode *top_node, sa_hist_mangager_ptr hm)
{
  if (Verbosity())
    cout << "saModSimpleHist::init_run - N beam position records = "
        << _beam_pos_xy.position.size() << endl;

  RunHeader* _run_header;
  try
    {
      _run_header = TMutNode<RunHeader>::find_io_node(top_node, "RunHeader");
    }
  catch (...)
    {
      _run_header = NULL;
    }
  if (_run_header)
    {
      _run_num = _run_header->get_RunNumber();

      if (Verbosity() >= 1)
        {
          cout << "saModSimpleHist::init_run - from RunHeader : run number = "
              << _run_num << endl;
        }

    }
  else
    {
      recoConsts *rc = recoConsts::instance();

      _run_num = rc->get_IntFlag("RUNNUMBER");

      if (Verbosity() >= 1)
        {
          cout << "saModSimpleHist::init_run - from recoConsts : run number = "
              << _run_num << endl;
        }
    }

  // Get offsets
  if (!_beam_pos_xy.have_run(_run_num))
    {
      cout
          << "saModSimpleHist::InitRun - Error - Do not have average beam XY for run "
          << _run_num << ". Use 0 cm instead. N records = "
          << _beam_pos_xy.position.size() << endl;

      _VTX_Avg_X = 0;
      _VTX_Avg_Y = 0;
    }
  else
    {

      _VTX_Avg_X = _beam_pos_xy.get_x(_run_num);
      _VTX_Avg_Y = _beam_pos_xy.get_y(_run_num);
      _VTX_Avg_dX = _beam_pos_xy.get_dx(_run_num);
      _VTX_Avg_dY = _beam_pos_xy.get_dy(_run_num);

      cout << "saModSimpleHist::InitRun - INFO - Use vertex [" << _VTX_Avg_X
          << ", " << _VTX_Avg_Y << "] cm + z * [" << _VTX_Avg_dX << ", "
          << _VTX_Avg_dY << "]" << endl;
    }

  return EVENT_OK;
}

//! event method
int
saModSimpleHist::event(PHCompositeNode *topNode, sa_hist_mangager_ptr hm)
{
  _h_Normalization->Fill(1);

  // VTX
  VtxOut* vtx = findNode::getClass<VtxOut>(topNode, "VtxOut");
  if (!vtx)
    {
      static bool once = true;

      if (once)
        {
          cout << "saModSimpleHist::event - Error - VtxOut not in Node Tree"
              << endl;
          once = false;
        }
    }
  else
    {

      const PHPoint _vertex = vtx->get_Vertex(_vertex_name.c_str());
      const bool bad_vertex = isnan(_vertex.getX()) || isnan(_vertex.getY())
          || isnan(_vertex.getZ()) || abs(_vertex.getZ()) > 15;

      const PHPoint _vertex_err = vtx->get_VertexError(_vertex_name.c_str());
      const bool precise_vertex = _vertex_err.getX() < 300e-4
          and _vertex_err.getY() < 300e-4 and _vertex_err.getZ() < 600e-4;

      if (!bad_vertex and precise_vertex)
        {
          _h_Normalization->Fill(2);
          _h_vtx->Fill(_vertex.getX(), _vertex.getY());
          _h_vtx_xz->Fill(_vertex.getZ(), _vertex.getX());
          _h_vtx_yz->Fill(_vertex.getZ(), _vertex.getY());
        }
//      else
//        return DISCARDEVENT;

      _event_vertex[2] = _vertex.getZ();
      _event_vertex[0] = _VTX_Avg_X + _VTX_Avg_dX * _event_vertex[2];
      _event_vertex[1] = _VTX_Avg_Y + _VTX_Avg_dY * _event_vertex[2];

    }

  // VTX
  SingleMuonContainer* SingleMuons = findNode::getClass<SingleMuonContainer>(
      topNode, "SingleMuonContainer");
  DiMuonContainer* DiMuons = findNode::getClass<DiMuonContainer>(topNode,
      "DiMuonContainer");
  if (!SingleMuons)
    {
      static bool once = true;

      if (once)
        {
          cout
              << "saModSimpleHist::event - Error - SingleMuons not in Node Tree"
              << endl;
          once = false;
        }
    }
  else if (!DiMuons)
    {
      static bool once = true;

      if (once)
        {
          cout
              << "saModSimpleHist::event - Error - DiMuonContainer not in Node Tree"
              << endl;
          once = false;
        }
    }
  else if (_VTX_Avg_X == 0 && _VTX_Avg_Y == 0)
    {
      static bool once = true;

      if (once)
        {
          cout << "saModSimpleHist::event - Error - no average vertex found"
              << endl;
          once = false;
        }
    }
  else
    {

      DCA(SingleMuons, DiMuons);
      DCA_FVTX(SingleMuons, DiMuons);
    }

  return EVENT_OK;
}

//! global termination
int
saModSimpleHist::end(PHCompositeNode *topNode, sa_hist_mangager_ptr hm)
{

  return EVENT_OK;
}

template<typename T>
  int
  sign(T val)
  {
    return (T(0) < val) - (val < T(0));
  }

void
saModSimpleHist::DCA(SingleMuonContainer * smc, DiMuonContainer * dmc)
{
  assert(smc);
  assert(dmc);

//  T->SetAlias("SingleMuons_arm", "1*(SingleMuons.pz>0)");
//  T->SetAlias("SingleMuons_arm_sign", "1*sign(SingleMuons.pz)");
//
//  T->SetAlias("disp_x_kalman",
//      "SingleMuons.x0_fvtxmutr + SingleMuons.px_fvtxmutr/SingleMuons.pz_fvtxmutr*(vertex_z - SingleMuons.z0_fvtxmutr) - vertex_x");
//  T->SetAlias("disp_y_kalman",
//      "SingleMuons.y0_fvtxmutr + SingleMuons.py_fvtxmutr/SingleMuons.pz_fvtxmutr*(vertex_z - SingleMuons.z0_fvtxmutr) - vertex_y");
//  T->SetAlias("phi",
//      "1*atan2(SingleMuons.py_fvtxmutr,SingleMuons.px_fvtxmutr)");
//
//  T->SetAlias("disp_pt_vtx_kalman",
//      "disp_x_kalman * cos(phi) + disp_y_kalman * sin(phi)");
//  T->SetAlias("disp_lateral_vtx_kalman",
//      "disp_x_kalman * cos(phi+pi/2) + disp_y_kalman * sin(phi+pi/2)");
//
//  // DCA corrections - FVTX workshop Oct 2012
//  T->SetAlias("disp_pt_vtx_kalman_cor",
//      "disp_pt_vtx_kalman + 0.1201 * disp_lateral_vtx_kalman * disp_lateral_vtx_kalman");

  //  const double vertex_x = smc->get_Evt_vtxX();
  //  const double vertex_y = smc->get_Evt_vtxY();
  const double vertex_z = _event_vertex[2];

  const double vertex_x = _event_vertex[0];
  const double vertex_y = _event_vertex[1];

  const double pi = TMath::Pi();

  //Cut
  if (isnan(vertex_z) || abs(vertex_z) > 15 || vertex_x == 0 || vertex_y == 0)
    return;

  for (unsigned int i = 0; i < smc->get_nSingleMuons(); i++)
    {
      SingleMuon * sm = smc->get_SingleMuon(i);
      assert(sm);
      _h_Normalization->Fill(3);

      // quality cuts
      const bool good_muID_muon = fabs(sm->get_DG0()) < 15
          and fabs(sm->get_DDG0()) < 10 and sm->get_lastgap() >= 4
          and sm->get_nidhits() > 5;

      double invmass = 0;
      for (unsigned int j = 0; j < dmc->get_nDiMuons(); j++)
        {
          DiMuon* dm = dmc->get_DiMuon(j);

          if (dm->get_charge() != 0)
            continue;

          if (abs(dm->get_Tr0_pz() - sm->get_pz()) < .1
              || abs(dm->get_Tr1_pz() - sm->get_pz()) < .1)
            invmass = dm->get_mass();
        }

      //Cut
      if (sm->get_dr_fvtx() <= 0)
        continue;
      if ((sm->get_hit_pattern() & 0xFF00) != 0 && _usv_svx_only == 0)
        continue;
      if ((sm->get_hit_pattern() & 0xFF00) == 0 && _usv_svx_only > 0)
        continue;

      const double SingleMuons_arm = 1 * (sm->get_pz() > 0);
//      const double SingleMuons_arm_sign = sign<float>(sm->get_pz());

//      const double proj_z = sm->get_pz() > 0 ? 30 : -30;
      const double proj_z = sm->get_pz() > 0 ? 40 : -40;

      const double disp_x_kalman = sm->get_x0_fvtxmutr()
          + sm->get_px_fvtxmutr() / sm->get_pz_fvtxmutr()
              * (vertex_z - sm->get_z0_fvtxmutr()) - vertex_x;
      const double disp_y_kalman = sm->get_y0_fvtxmutr()
          + sm->get_py_fvtxmutr() / sm->get_pz_fvtxmutr()
              * (vertex_z - sm->get_z0_fvtxmutr()) - vertex_y;

      const double proj_x_kalman = sm->get_x0_fvtxmutr()
          + sm->get_px_fvtxmutr() / sm->get_pz_fvtxmutr()
              * (proj_z - sm->get_z0_fvtxmutr());
      const double proj_y_kalman = sm->get_y0_fvtxmutr()
          + sm->get_py_fvtxmutr() / sm->get_pz_fvtxmutr()
              * (proj_z - sm->get_z0_fvtxmutr());

      const double phi = atan2(proj_y_kalman, proj_x_kalman);

      const double disp_pt_vtx_kalman = disp_x_kalman * cos(phi)
          + disp_y_kalman * sin(phi);
      const double disp_lateral_vtx_kalman = disp_x_kalman * cos(phi + pi / 2)
          + disp_y_kalman * sin(phi + pi / 2);

      // DCA corrections - FVTX workshop Oct 2012
//      const double disp_pt_vtx_kalman_cor = disp_pt_vtx_kalman
//          + 0.1201 * disp_lateral_vtx_kalman * disp_lateral_vtx_kalman;

      if (abs(invmass - 3.11) < .3)
        {

          _h_Normalization -> Fill(4);

          if (good_muID_muon)
            {

              _h_Normalization -> Fill(5);

              // JPsi muons
              _h_DCA_Pz_Arm->Fill(disp_pt_vtx_kalman,
                  abs(sm->get_pz_fvtxmutr()), SingleMuons_arm);
              _h_DCA_Pz_Arm->Fill(disp_lateral_vtx_kalman,
                  abs(sm->get_pz_fvtxmutr()), SingleMuons_arm + 2);

              if (SingleMuons_arm == 0)
                {
                  _h_DCA_Phi_z_Arm0->Fill(disp_pt_vtx_kalman, phi, vertex_z);
                }
              else
                {
                  _h_DCA_Phi_z_Arm1->Fill(disp_pt_vtx_kalman, phi, vertex_z);
                }

              if (abs(sm->get_pz_fvtxmutr()) > 3)
                {
                  _h_DCA_Phi_Arm->Fill(disp_pt_vtx_kalman, phi, SingleMuons_arm);
                  _h_DCA_Phi_Arm->Fill(disp_lateral_vtx_kalman, phi,
                      SingleMuons_arm + 2);
                }
            }

        }

      if (good_muID_muon and abs(sm->get_pz_fvtxmutr()) > 6)
        {

          _h_Normalization -> Fill(6);
          if (SingleMuons_arm == 0)
            {
              _h_DCA_Phi_z_Arm0->Fill(disp_pt_vtx_kalman, phi, vertex_z);
            }
          else
            {
              _h_DCA_Phi_z_Arm1->Fill(disp_pt_vtx_kalman, phi, vertex_z);
            }

        }

    }

}

void
saModSimpleHist::DCA_FVTX(SingleMuonContainer * smc, DiMuonContainer * dmc)
{
  assert(smc);
  assert(dmc);

//  T->SetAlias("SingleMuons_arm", "1*(SingleMuons.pz>0)");
//  T->SetAlias("SingleMuons_arm_sign", "1*sign(SingleMuons.pz)");
//
//  T->SetAlias("disp_x_kalman",
//      "SingleMuons.x0_fvtxmutr + SingleMuons.px_fvtxmutr/SingleMuons.pz_fvtxmutr*(vertex_z - SingleMuons.z0_fvtxmutr) - vertex_x");
//  T->SetAlias("disp_y_kalman",
//      "SingleMuons.y0_fvtxmutr + SingleMuons.py_fvtxmutr/SingleMuons.pz_fvtxmutr*(vertex_z - SingleMuons.z0_fvtxmutr) - vertex_y");
//  T->SetAlias("phi",
//      "1*atan2(SingleMuons.py_fvtxmutr,SingleMuons.px_fvtxmutr)");
//
//  T->SetAlias("disp_pt_vtx_kalman",
//      "disp_x_kalman * cos(phi) + disp_y_kalman * sin(phi)");
//  T->SetAlias("disp_lateral_vtx_kalman",
//      "disp_x_kalman * cos(phi+pi/2) + disp_y_kalman * sin(phi+pi/2)");
//
//  // DCA corrections - FVTX workshop Oct 2012
//  T->SetAlias("disp_pt_vtx_kalman_cor",
//      "disp_pt_vtx_kalman + 0.1201 * disp_lateral_vtx_kalman * disp_lateral_vtx_kalman");

  const double vertex_z = _event_vertex[2];

  const double vertex_x = _event_vertex[0];
  const double vertex_y = _event_vertex[1];

  const double pi = TMath::Pi();

  //Cut
  if (isnan(vertex_z) || abs(vertex_z) > 15 || vertex_x == 0 || vertex_y == 0)
    return;

  for (unsigned int i = 0; i < smc->get_nSingleMuons(); i++)
    {
      SingleMuon * sm = smc->get_SingleMuon(i);

      double invmass = 0;
      for (unsigned int j = 0; j < dmc->get_nDiMuons(); j++)
        {
          DiMuon* dm = dmc->get_DiMuon(j);

          if (dm->get_charge() != 0)
            continue;

          if (abs(dm->get_Tr0_pz() - sm->get_pz()) < .1
              || abs(dm->get_Tr1_pz() - sm->get_pz()) < .1)
            invmass = dm->get_mass();
        }
      //Cut
      if (sm->get_dr_fvtx() <= 0)
        continue;
      if ((sm->get_hit_pattern() & 0xFF00) != 0 && _usv_svx_only == 0)
        continue;
      if ((sm->get_hit_pattern() & 0xFF00) == 0 && _usv_svx_only > 0)
        continue;
      if (sm->get_x0_fvtx() == 0 && sm->get_y0_fvtx() == 0)
        continue;

      const double SingleMuons_arm = 1 * (sm->get_pz() > 0);
//      const double SingleMuons_arm_sign = sign<float>(sm->get_pz());

      const double disp_x_kalman = sm->get_x0_fvtx()
          + sm->get_px_fvtx() / sm->get_pz_fvtx()
              * (vertex_z - sm->get_z0_fvtx()) - vertex_x;
      const double disp_y_kalman = sm->get_y0_fvtx()
          + sm->get_py_fvtx() / sm->get_pz_fvtx()
              * (vertex_z - sm->get_z0_fvtx()) - vertex_y;
      const double phi = atan2(sm->get_py_fvtx(), sm->get_px_fvtx());

      const double disp_pt_vtx_kalman = disp_x_kalman * cos(phi)
          + disp_y_kalman * sin(phi);
      const double disp_lateral_vtx_kalman = disp_x_kalman * cos(phi + pi / 2)
          + disp_y_kalman * sin(phi + pi / 2);

      // DCA corrections - FVTX workshop Oct 2012
      const double disp_pt_vtx_kalman_cor = disp_pt_vtx_kalman
          + 0.1201 * disp_lateral_vtx_kalman * disp_lateral_vtx_kalman;

      if (abs(invmass - 3.11) < .2)
        if ((sm->get_lastgap()) >= 4 && sm->get_DG0() < 20
            && sm->get_DDG0() < 9)
          {
            _h_fvtx_DCA_Pz_Arm->Fill(disp_pt_vtx_kalman_cor,
                abs(sm->get_pz_fvtxmutr()), SingleMuons_arm);
            _h_fvtx_DCA_Pz_Arm->Fill(disp_lateral_vtx_kalman,
                abs(sm->get_pz_fvtxmutr()), SingleMuons_arm + 2);
          }

      if (abs(sm->get_pz_fvtxmutr()) > 3)
        {
          _h_fvtx_DCA_Phi_Arm->Fill(disp_pt_vtx_kalman_cor, phi,
              SingleMuons_arm);
          _h_fvtx_DCA_Phi_Arm->Fill(disp_lateral_vtx_kalman, phi,
              SingleMuons_arm + 2);
        }
    }

}

//______________________________________________________
void
saModSimpleHist::BeamPosXY::load_data(string data_file, int verbose)
{
  cout << "saModSimpleHist::BeamPosXY::load_data - Info - loading data file "
      << data_file << endl;

  fstream fdata;
  try
    {
      fdata.open(data_file.c_str(), ios_base::in);
    }
  catch (const std::exception& e)
    {
      cout << "saModSimpleHist::BeamPosXY::load_data - Error - cannot open "
          << data_file << e.what() << endl;
      return;
    }

  if (!fdata.is_open())
    {
      cout
          << "saModSimpleHist::BeamPosXY::load_data - Error - cannot open data file "
          << data_file << endl;
      return;
    }

  while (!fdata.eof())
    {
      int run = -9999;
      double x = -9999, y = -9999, ex = -9999, ey = -9999;
      double dx = -0, dy = -0, edx = -9999, edy = -9999;

      fdata >> run >> x >> y >> ex >> ey;
      fdata >> dx >> dy >> edx >> edy;
      if (verbose >= 2)
        cout
            << "saModSimpleHist::BeamPosXY::load_data - Info - Single record for beam positon of run "
            << run << " = [" << x << " , " << y << "] +/- " << " [" << ex
            << " , " << ey << "] cm" << endl;

      if (ex <= 0)
        ex = 1e-6;
      if (ey <= 0)
        ey = 1e-6;

      const double wx = 1 / ex / ex;
      const double wy = 1 / ey / ey;
      const double wdx = 1 / edx / edx;
      const double wdy = 1 / edy / edy;

      if (position.count(run))
        {
          const double oldx = position[run].first;
          const double oldy = position[run].second;
          const double oldwx = weight[run].first;
          const double oldwy = weight[run].second;

          position[run] = pos((x * wx + oldx * oldwx) / (wx + oldwx),
              (y * wy + oldy * oldwy) / (wy + oldwy));
          weight[run] = pos(wx + oldwx, wy + oldwy);
        }
      else
        {
          position[run] = pos(x, y);
          weight[run] = pos(wx, wy);
        }

      if (dposition.count(run))
        {
          const double oldx = dposition[run].first;
          const double oldy = dposition[run].second;
          const double oldwx = dweight[run].first;
          const double oldwy = dweight[run].second;

          dposition[run] = pos((dx * wdx + oldx * oldwx) / (wdx + oldwx),
              (dy * wdy + oldy * oldwy) / (wdy + oldwy));
          dweight[run] = pos(wdx + oldwx, wdy + oldwy);
        }
      else
        {
          position[run] = pos(x, y);
          weight[run] = pos(wx, wy);
          dposition[run] = pos(dx, dy);
          dweight[run] = pos(wdx, wdy);
        }
    }
  fdata.close();

  for (record::const_iterator iter = position.begin(); iter != position.end();
      iter++)
    {
      const int run = iter->first;
      pos xy = iter->second;
      pos xy_weight = weight[run];

      if (verbose)
        cout
            << "saModSimpleHist::BeamPosXY::load_data - Info - beam positon for run "
            << run << " = [" << xy.first << " , " << xy.second << "] +/- "
            << " [" << 1 / sqrt(xy_weight.first) << " , "
            << 1 / sqrt(xy_weight.second) << "] cm " << have_run(run) << endl;
    }
}

//________________________________________
bool
saModSimpleHist::BeamPosXY::have_run(int run)
{
  return position.count(run) > 0;
}

//________________________________________
double
saModSimpleHist::BeamPosXY::get_x(int run)
{
  return position[run].first;
}

//________________________________________
double
saModSimpleHist::BeamPosXY::get_y(int run)
{
  return position[run].second;
}

//________________________________________
double
saModSimpleHist::BeamPosXY::get_dx(int run)
{
  return dposition[run].first;
}

//________________________________________
double
saModSimpleHist::BeamPosXY::get_dy(int run)
{
  return dposition[run].second;
}

