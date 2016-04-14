using namespace std;
typedef std::vector<int> n_list;


void
Draw_vs_FILL_RUN_work(
		TFile *output_file
		,TFile *input_file
		,std::string sahBase_name
		,bool run_or_fill
		,std::string hist_type = "PP"
		)
{
	input_file->cd();
	saHist * sahBase = (saHist *) input_file->GetObjectChecked(sahBase_name.data(),"saHist");
	if (!sahBase){
		cout<<__FILE__<<" : "<<__LINE__<<endl;
		continue;
	}

	sahBase->AutoLoad(0);

	//assert(sahBase);

	n_list src_list;

	if (run_or_fill)
		src_list = sahBase->get_run_list();
	else
		src_list = sahBase->get_fill_list();

	TGraphErrors * ge_type = new TGraphErrors(src_list.size());

	const int n_bin_x = sahBase->getHisto("BASE")->GetNbinsX();
	const int n_bin_y = sahBase->getHisto("BASE")->GetNbinsY();
	const int n_bin_z = sahBase->getHisto("BASE")->GetNbinsZ();

	for (int x = 1; x <= n_bin_x; x++)
	{
		for (int y = 1; y <= n_bin_y; y++)
		{
			for (int z = 1; z <= n_bin_z; z++)
			{
				int i_ge_bin = 0;
				for (n_list::const_iterator iter = src_list.begin();
						iter != src_list.end(); ++iter)
				{
					const int n = (*iter);

					saHist *sahSub = NULL;
					if (run_or_fill)
						sahSub = sahBase->getsaHist_run(n);
					else
						sahSub = sahBase->getsaHist_fill(n);

					if (!sahSub){
						cout 
							<< "saHist::" << sahBase->get_name()
							<< "::Draw_vs_FILL_RUN_work - Error - can not load sub saHist "
							<< endl;
						continue;
					}	else {
						// read the saHist here
						TH1 *h_type = sahSub->getHisto(hist_type.data());

						(ge_type->GetX())[i_ge_bin] = n;
						(ge_type->GetY())[i_ge_bin] = h_type->GetBinContent(x, y, z);
						(ge_type->GetEY())[i_ge_bin] = h_type->GetBinError(x, y, z);
						i_ge_bin++;
					}
				}

				TFitResultPtr t_fit_result_ptr = ge_type->Fit("pol0", "MSQ");
				cout<<t_fit_result_ptr->Parameter(0)<<" +- "<<t_fit_result_ptr->ParError(0)<<endl;

				std::string name_base = sahBase->get_name();
				std::string name_suffix = run_or_fill ? "RUN_FIT" : "FILL_FIT";
				ge_type->SetName(Form("%s_%d_%d_%d",(name_base+"_"+hist_type+"_"+name_suffix).data(),x,y,z));

				output_file->cd();
				ge_type->Write();
			}
		}
	}
}

void
Draw_vs_FILL_RUN(
		std::string input_name = "picoDST.list_hist.root_A_LL.root"
		, bool run_or_fill = false)
{
	TFile *input_file = new TFile(input_name.data(),"read");
	TFile *output_file = new TFile("Draw_vs_FILL_RUN.root","recreate");

	std::vector<std::string> v_saHist_names;
	//v_saHist_names.push_back("pT_cb_bkg_n_saHist");
	//v_saHist_names.push_back("pT_cb_bkg_s_saHist");
	//v_saHist_names.push_back("pT_os_lsb_n_saHist");
	//v_saHist_names.push_back("pT_os_lsb_s_saHist");
	//v_saHist_names.push_back("pT_os_sig_n_saHist");
	//v_saHist_names.push_back("pT_os_sig_s_saHist");
	//v_saHist_names.push_back("pT_ss_all_n_saHist");
	//v_saHist_names.push_back("pT_ss_all_s_saHist");
	//v_saHist_names.push_back("RelLumi_saHist");

	Draw_vs_FILL_RUN_work(output_file, input_file, "pT_os_sig_n_saHist", run_or_fill, "A_LL");
	//Draw_vs_FILL_RUN_work(output_file, input_file, "pT_os_sig_n_saHist", run_or_fill, "PP");
	//Draw_vs_FILL_RUN_work(output_file, input_file, "RelLumi_saHist", run_or_fill, "PP");
	//Draw_vs_FILL_RUN_work(output_file, input_file, "RelLumi_saHist", run_or_fill, "Pol2_Blue_Yellow");

	input_file->Close();
	output_file->Close();
}
