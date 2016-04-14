// $Id: saModuleDimuonJpsiHaiwang.h,v 1.4 2014/12/11 22:33:08 yuhw Exp $                                                                                             

/*!
 * \file saModuleDimuonJpsiHaiwang.h
 * \brief saModuleSimpleDimuon 
 *		modified from Ming's saModuleDimuonMing. for DiMuon Jpsi Analysis 
 * \author Haiwang Yu <yuhw@rcf.rhic.bnl.gov>
 * \version $Revision: 1.4 $
 * \date $Date: 2014/12/11 22:33:08 $
 */

#ifndef SAMODULEDIMUONJPSIHAIWANG_H_
#define SAMODULEDIMUONJPSIHAIWANG_H_

#include <iostream>
#include <iomanip>
#include <fstream>

#include "saModuleBase.h"

class DiMuon;
class RpcDiMuon;

class saModuleDimuonJpsiHaiwang  : public saModuleBase
{
	public:
		saModuleDimuonJpsiHaiwang(const std::string &name = "DimuonJpsiHaiwang", bool doControlData = false);
		virtual
			~saModuleDimuonJpsiHaiwang();

		int verbosity;

		bool _use_fixed_mass_win;
		bool _use_bbc_cut;
		bool _use_rpc_cut;

		void set_use_bbc_cut(const bool a){_use_bbc_cut = a;}
		void set_use_rpc_cut(const bool a){_use_rpc_cut = a;}

		std::ofstream fout;

		bool _doControlData;
		int _FillCounter;

		float JPSI_SIG_MASS_MIN_N;
		float JPSI_SIG_MASS_MAX_N;
		float JPSI_SIG_MASS_MIN_S;
		float JPSI_SIG_MASS_MAX_S;
		float JPSI_LSB_MASS_MIN;
		float JPSI_LSB_MASS_MAX;
		float JPSI_ALL_MASS_MIN;
		float JPSI_ALL_MASS_MAX;

		float _pT_MinJpsiMs[2][3];
		float _pT_MaxJpsiMs[2][3];
		double _pT_bins[4];
		std::ifstream ifs_fitparam;
		void set_fitparam_file_name(const std::string a){ifs_fitparam.open(a.data());}
		bool set_fitparam(int nSigma = 2);

		int MIN_FILL_NUM;
		int MAX_FILL_NUM;

		enum EVENT_TYPE {OS_SIG_N, OS_SIG_S, OS_LSB_N, OS_LSB_S, SS_ALL_N, SS_ALL_S, BAD_EVENT};
		enum ARM_TYPE {NORTH, SOUTH};

#ifndef __CINT__

		//! global initialization
		virtual int
			init(PHCompositeNode *topNode, sa_hist_mangager_ptr hm);

		//! Run initialization
		virtual int
			init_run(PHCompositeNode *topNode, sa_hist_mangager_ptr hm);

		//! event method
		virtual int
			event(PHCompositeNode *topNode, sa_hist_mangager_ptr hm);

		//! global termination
		virtual int
			end(PHCompositeNode *topNode, sa_hist_mangager_ptr hm);

		//! dimuon cut
		virtual bool
			pass_dimuon_cut(const DiMuon *dimuon, const RpcDiMuon *rpcdimuon);

#endif

		// -- Same Arm dimuons
		saHist * _h_Xing_Check;// for Xing check
		// TH1D bined by pT
		saHist * _h_pT_os_sig_n; // opposite sign, signal region, north
		saHist * _h_pT_os_sig_s; // opposite sign, signal region, south

		saHist * _h_pT_os_lsb_n; // opposite sign, lower sideband, north
		saHist * _h_pT_os_lsb_s; // opposite sign, lower sideband, south

		saHist * _h_pT_ss_all_n; // same sign, all region, north
		saHist * _h_pT_ss_all_s; // same sign, all region, south

		saHist * _h_pT_cb_bkg_n; // combined bkg, north
		saHist * _h_pT_cb_bkg_s; // combined bkg, south
};

#endif /* SAMODULEDIMUONJPSIHAIWANG_H_ */
