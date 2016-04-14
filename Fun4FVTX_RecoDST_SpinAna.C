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
Fun4FVTX_RecoDST_SpinAna(
		int nEvents = 0, //
//		char *input_file = "picoDST_test.list",
		char *input_file = "picoDST.list", //
		bool doControlData = true,
		bool use_bbc_cut = true,
		char *dst_file = NULL //
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
	gSystem->Load("libSvxDstQA.so");
	gSystem->Load("libpicodst_object.so");
	gSystem->Load("/gpfs/mnt/gpfs02/phenix/spin/spin3/liuld/codefile/spinAnalyser/install/lib/libspin_analyzer.so");
		
	ifstream readfile;
	string file_name;
	string s1,s2;
	readfile.open(input_file,ios::in);
	getline(readfile,file_name);
	size_t ipos = file_name.find(".");
	s1 = file_name.substr(ipos-6,6);
	cout <<" ========================="<<s1<<endl;

	///////////////////////////////////////////
	// Make the Server
	//////////////////////////////////////////

	Fun4AllServer *se = Fun4AllServer::instance();
	se->Verbosity(0);

	//! load a QA text file marking good run and good corssings
	//saModLoadSpinInfo * mspin = new saModLoadSpinInfo("SpinInfoWithAdditionalQA");
	//mspin -> load_crossing_flag("spinQA.txt");
	//saSpinAnalyzer * sa = new saSpinAnalyzer("SpinAnalyzer", mspin);

	//! use the database
	saSpinAnalyzer * sa = new saSpinAnalyzer("SpinAnalyzer");

	sa->set_auto_save_hist(
			string(s1) + string("_") + string("hist.root"));
	sa->Verbosity(0);

	saModuleJpsiAN *samodule = new saModuleJpsiAN("saModuleJpsiAN",doControlData);
	samodule->set_use_bbc_cut(use_bbc_cut);
	samodule->verbosity = 0;
	sa->RegisterModule(samodule);

	//saModuleDimuonJpsiHaiwang *samodule = new saModuleDimuonJpsiHaiwang("saModuleDimuonJpsiHaiwang",doControlData);
	//samodule->set_use_bbc_cut(use_bbc_cut);
	//samodule->verbosity = 0;
	//sa->RegisterModule(samodule);

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
		dstManager->AddNode("SimpleDimuonFlag");
		dstManager->AddNode("Sync");
		dstManager->AddNode("TrigLvl1");

		dstManager->AddNode("DiMuonContainer");
		dstManager->AddNode("RunHeader");

		dstManager->AddEventSelector("SpinAnalyzer");
		se->registerOutputManager(dstManager);
	}
	///////////////////////////////////////////
	// Analyze the Data.
	//////////////////////////////////////////

	gSystem->ListLibraries();

	Fun4AllDstInputManager *in = new Fun4AllDstInputManager("DSTin");
	//in->fileopen(input_file);
	in->AddListFile(input_file);
//	in->AddListFile(dst_file);
	se->registerInputManager(in);

	se->run(nEvents);
	se->End();

	delete se;

	cout << "Completed reconstruction." << endl;
}
