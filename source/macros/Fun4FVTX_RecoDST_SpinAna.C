/*!
 * \file Fun4FVTX_RecoDST_SpinAna.C
 * \brief example of how to use spin analyzer
 * \author Jin Huang <jhuang@bnl.gov>
 * \version $Revision: 1.4 $
 * \date $Date: 2013/07/24 07:43:27 $
 */

#include <exception>
using namespace std;

void
Fun4FVTX_RecoDST_SpinAna(int nEvents = 0, //
//    char *input_file = "data/366636_0.lst", //
//    char *input_file = "data/taxi_1264.filelist", //
//        char *input_file = "data/taxi_1264_100.filelist", //
    char *input_file = "data/taxi_1264_3.filelist", //
    char *dst_file = "dst_out.root" //
    )
{

  // load libraries
  gSystem->Load("libmutoo_subsysreco.so");
  gSystem->Load("libfun4all.so");
  gSystem->Load("librecal.so");
  gSystem->Load("libfun4allfuncs.so");
  gSystem->Load("liblpc.so");
  gSystem->Load("libcompactCNT.so");
  gSystem->Load("libfun4allfuncs_muons");
  gSystem->Load("libMWGOO");
  gSystem->Load("libmutrg");
  gSystem->Load("librpc_subsysreco");
  gSystem->Load("libSvxDstQA.so");
  gSystem->Load("libpicodst_object");
  gSystem->Load("libspin_analyzer");

  ///////////////////////////////////////////
  // Make the Server
  //////////////////////////////////////////

  Fun4AllServer *se = Fun4AllServer::instance();
  se->Verbosity(0);

  //! load a QA text file marking good run and good corssings
  saModLoadSpinInfo * mspin = new saModLoadSpinInfo("SpinInfoWithAdditionalQA");
  mspin -> load_crossing_flag("spinQA.txt");

  saSpinAnalyzer * sa = new saSpinAnalyzer("SpinAnalyzer", mspin);

  sa->set_auto_save_hist(
      string(input_file) + string("_") + string("hist.root"));
  sa->Verbosity(1);

  sa->RegisterModule(new saModuleSimpleDimuon());
//  saModSimpleHist * sh = new saModSimpleHist();
//  sh -> Verbosity(1);
////  sh -> load_beam_xy_data("data/taxi_1627.list.Smooth.GetBeamPos");
//  sh -> load_beam_xy_data("data/taxi_1627.list.Smooth_Shifted.GetBeamPos");
//  sa->RegisterModule(sh);

  se->registerSubsystem(sa);

  ///////////////////////////////////////////
  // DST
  //////////////////////////////////////////
  if (dst_file)
    {
      std::cout << "registering Fun4AllDstOutputManager" << std::endl;

      Fun4AllDstOutputManager *dstManager = new Fun4AllDstOutputManager(
          "DSTOUT", string(input_file) + string("_") + string(dst_file));

      dstManager->AddNode("saEventProperty");
      dstManager->AddNode("DiMuonContainer");
      dstManager->AddNode("SimpleDimuonFlag");
      dstManager->AddNode("Sync");
      dstManager->AddNode("TrigLvl1");


      dstManager->AddEventSelector("SpinAnalyzer");
      se->registerOutputManager(dstManager);
    }
  ///////////////////////////////////////////
  // Analyze the Data.
  //////////////////////////////////////////

  gSystem->ListLibraries();

  Fun4AllDstInputManager *in = new Fun4AllDstInputManager("DSTin");
//  in->fileopen(input_file);
  in->AddListFile(input_file);
  se->registerInputManager(in);

  se->run(nEvents);
  se->End();

  delete se;

  cout << "Completed reconstruction." << endl;
}
