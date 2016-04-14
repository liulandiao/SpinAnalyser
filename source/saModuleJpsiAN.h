// $Id: saModuleJpsiAN.h,v 1.3 2015/09/03 20:45:38 yuhw Exp $                                                                                             

/*!
 * \file saModuleJpsiAN.h
 * \brief saModuleSimpleDimuon 
 *		modified from Ming's saModuleDimuonMing. for DiMuon Jpsi Analysis 
 * \author Haiwang Yu <yuhw@rcf.rhic.bnl.gov>
 * \version $Revision: 1.3 $
 * \date $Date: 2015/09/03 20:45:38 $
 */

#ifndef _SAMODULEJPSIAN_H_
#define _SAMODULEJPSIAN_H_

//STL
#include <iostream>
#include <iomanip>
#include <fstream>
#include <vector>

#include "saModuleBase.h"

class DiMuon;

class saModuleJpsiAN  : public saModuleBase
{
	public:
		saModuleJpsiAN(const std::string &name = "DimuonJpsiHaiwang", bool doControlData = false);
		virtual
			~saModuleJpsiAN();

		int verbosity;

		bool _use_fixed_mass_win;
		bool _use_bbc_cut;

		void set_use_bbc_cut(const bool a){_use_bbc_cut = a;}

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

		float _pT_MinJpsiMs[2][2];
		float _pT_MaxJpsiMs[2][2];
		double _pT_bins[2];
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
			pass_dimuon_cut(const DiMuon *dimuon);

		TH1 *makeTHist_pT_phi(
				std::string name,
				std::string title
				);
		saHist *makesaHist_pT_phi(
				sa_hist_mangager_ptr hm,
				TH1 *h
				);

#endif

		// TH2D bined by pT phi
		saHist * _sah_pT_os_sig_n; // opposite sign, signal region, north
		saHist * _sah_pT_os_sig_s; // opposite sign, signal region, south

		saHist * _sah_pT_os_lsb_n; // opposite sign, lower sideband, north
		saHist * _sah_pT_os_lsb_s; // opposite sign, lower sideband, south

		std::vector<saHist*> _v_sahist;
};

#endif /* _SAMODULEJPSIAN_H_ */
