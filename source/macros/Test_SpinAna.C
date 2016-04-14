
// $Id: Test_SpinAna.C,v 1.2 2013/10/09 22:31:26 jinhuang Exp $

/*!
 * \file Test_SpinAna.C
 * \brief Example how to draw the result of spin analyzer
 * \author Jin Huang <jhuang@bnl.gov>
 * \version $Revision: 1.2 $
 * \date $Date: 2013/10/09 22:31:26 $
 */

#include "SaveCanvas.C"
#include "SetOKStyle.C"
#include <cassert>
using namespace std;

void
Test_SpinAna()
{
  SetOKStyle();
  gStyle->SetOptStat(0);
  gStyle->SetOptFit(1111);
  TVirtualFitter::SetDefaultFitter("Minuit2");

  gSystem->Load("libspin_analyzer");

  ReadTest();
}

void
ReadTest()
{

  //  TFile *_file0 = TFile::Open("data/366636_0.lst_hist.root");
  //  TFile *_file0 = TFile::Open("data/taxi_1264.filelist_hist.root");
  //  TFile *_file0 = TFile::Open("data/taxi_1264_100.filelist_hist.root");
  TFile *_file0 = TFile::Open("data/taxi_1264_2.filelist_hist.root");

  saHist * InvMass_saHist = (saHist *) _file0->GetObjectChecked(
      "InvMass_saHist", "saHist");
  assert(InvMass_saHist);
  saHist * RelLumi_saHist = (saHist *) _file0->GetObjectChecked(
      "RelLumi_saHist", "saHist");
  assert(RelLumi_saHist);

  InvMass_saHist->AutoLoad();

  RelLumi_saHist->AutoLoad();

  ///////////////////////

  TCanvas *c1 = new TCanvas("Test_SpinAna_ReadTest", "Test_SpinAna_ReadTest",
      1000, 900);
  c1->Divide(2, 2);
  int idx = 1;
  TPad * p;

  p = (TPad *) c1->cd(idx++);
  c1->Update();

  p->SetLogy();

  InvMass_saHist->getHisto("YIELD")->Draw();
  InvMass_saHist->getHisto("YIELD")->SetYTitle("Yield");

  p = (TPad *) c1->cd(idx++);
  c1->Update();

  InvMass_saHist->getHisto("A_LL")->Draw();
  InvMass_saHist->getHisto("A_LL")->SetYTitle("Double spin asymmetry");

  p = (TPad *) c1->cd(idx++);
  c1->Update();

  InvMass_saHist->getHisto("A_L_Blue")->Draw();
  InvMass_saHist->getHisto("A_L_Blue")->SetYTitle(
      "Single spin asymmetry (Blue)");

  p = (TPad *) c1->cd(idx++);
  c1->Update();

  InvMass_saHist->getHisto("A_L_Yellow")->Draw();
  InvMass_saHist->getHisto("A_L_Yellow")->SetYTitle(
      "Single spin asymmetry (Yelllow)");

  c1->Paint();

  SaveCanvas(c1,
      TString(_file0->GetName()) + TString("_") + TString(c1->GetName()),
      kFALSE);

}
