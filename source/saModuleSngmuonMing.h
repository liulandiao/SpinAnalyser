// $Id: saModuleSngmuonMing.h,v 1.3 2013/08/17 05:41:42 mxliu Exp $                                                                                             
// 
// modefied from Jin's saModuleSimpleDimuon example   7/12/2013  MXL
//
// for Run12+ single muons analyses: light hadron, heavy flavor, W
//
/////////////////////////////////////////////////////////////////////

#ifndef SAMODULESNGMUONMING_H_
#define SAMODULESNGMUONMING_H_


#include "saModuleBase.h"

/*!
 * \brief saModuleSimpleDimuon
 */
class saModuleSngmuonMing : public saModuleBase
{
public:
  saModuleSngmuonMing(const std::string &name = "SngMuonMing");
  virtual
  ~saModuleSngmuonMing();

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

#endif
  // -- general event/track histograms 

  saHist * _h_ST1Z_N;         // St1K multiple scattering angle 
  saHist * _h_ST1Z_S;         // St1 hit Z

  saHist * _h_dA;         // St1 multiple scattering angle 
  saHist * _h2_dA;        // St1 multiple scattering angle vs pZ 

  saHist * _h_phi_1N;         // St1 phi angle 
  saHist * _h_phi_1S;         // St1 phi angle 

  saHist * _h_phi_2N;         // St2 phi angle 
  saHist * _h_phi_2S;         // St2 phi angle 

  saHist * _h_phi_3N;         // St3 phi angle 
  saHist * _h_phi_3S;         // St3 phi angle 

  saHist * _h_phi_13N;         // St1-3 phi angle 
  saHist * _h_phi_13S;         // St1-3 phi angle 

  saHist * _h2_phi_13N;         // St1-3 phi angle 
  saHist * _h2_phi_13S;         // St1-3 phi angle 

  saHist * _h_phi_23N;         // St1-3 phi angle 
  saHist * _h_phi_23S;         // St1-3 phi angle 

  saHist * _h2_phi_23N;         // St1-3 phi angle 
  saHist * _h2_phi_23S;         // St1-3 phi angle 

  // -- RPC stuff 
  saHist * _h_Rpc1DCA_N; // RPC1 DCA north 
  saHist * _h_Rpc1DCA_S; // RPC1 DCA south 

  saHist * _h_Rpc3DCA_N; // RPC3 DCA north 
  saHist * _h_Rpc3DCA_S; // RPC3 DCA south 


  //for comparision with run11/12 style W->muon analysis 
  saHist * _h_dw23N_p;         // charge=+; St2-3 rescaled phi angle == pT*sin(theta)*dPhi23
  saHist * _h_dw23N_m;         // charge=-;St2-3 rescaled phi angle == pT*sin(theta)*dPhi23
  saHist * _h_dw23S_p;         // St2-3 rescaled phi angle 
  saHist * _h_dw23S_m;         // St2-3 rescaled phi angle 

  saHist * _h_dw23N_PT10_p;         // charge=+; St2-3 rescaled phi angle == pT*sin(theta)*dPhi23
  saHist * _h_dw23N_PT10_m;         // charge=-;St2-3 rescaled phi angle == pT*sin(theta)*dPhi23
  saHist * _h_dw23S_PT10_p;         // St2-3 rescaled phi angle 
  saHist * _h_dw23S_PT10_m;         // St2-3 rescaled phi angle 

  saHist * _h_dw23N_PT15_p;         // charge=+; St2-3 rescaled phi angle == pT*sin(theta)*dPhi23
  saHist * _h_dw23N_PT15_m;         // charge=-;St2-3 rescaled phi angle == pT*sin(theta)*dPhi23
  saHist * _h_dw23S_PT15_p;         // St2-3 rescaled phi angle 
  saHist * _h_dw23S_PT15_m;         // St2-3 rescaled phi angle 

  saHist * _h_dw23N_PT20_p;         // charge=+; St2-3 rescaled phi angle == pT*sin(theta)*dPhi23
  saHist * _h_dw23N_PT20_m;         // charge=-;St2-3 rescaled phi angle == pT*sin(theta)*dPhi23
  saHist * _h_dw23S_PT20_p;         // St2-3 rescaled phi angle 
  saHist * _h_dw23S_PT20_m;         // St2-3 rescaled phi angle 

  // -- eta for W/Z analysis, for distribution check only, with fine bins --  
  saHist * _hx_etaW10;
  saHist * _hx_etaW10_p;      // Q+ high pT>10 tracks
  saHist * _hx_etaW10_m;      // Q- high pT>10 tracks

  saHist * _hx_etaW15;
  saHist * _hx_etaW15_p;      // Q+ high pT>15 tracks
  saHist * _hx_etaW15_m;      // Q- high pT>15 tracks

  saHist * _hx_etaW20;
  saHist * _hx_etaW20_p;      // Q+ high pT>20 tracks
  saHist * _hx_etaW20_m;      // Q- high pT>20 tracks


  //
  // --- the following histograms are used to calculate the spin asymmetry ---
  //

  // -- Single Deep muons, lastGap =>4

  saHist * _h_pT;
  saHist * _h_pT_n;      // north arm
  saHist * _h_pT_s;      // south arm

  saHist * _h_pT_p;      // Q+ tracks
  saHist * _h_pT_p_n;    // north arm
  saHist * _h_pT_p_s;    // south arm

  saHist * _h_pT_m;      // Q- tracks
  saHist * _h_pT_m_n;    // north arm
  saHist * _h_pT_m_s;    // south arm


  // -- Single Stooped hadrons, lastGap <=3

  saHist * _h_pTH;
  saHist * _h_pTH_n;      // north arm
  saHist * _h_pTH_s;      // south arm

  saHist * _h_pTH_p;      // Q+ tracks
  saHist * _h_pTH_p_n;    // north arm
  saHist * _h_pTH_p_s;    // south arm

  saHist * _h_pTH_m;      // Q- tracks
  saHist * _h_pTH_m_n;    // north arm
  saHist * _h_pTH_m_s;    // south arm


  // -- eta for heavy flavor --  

  saHist * _h_eta;
  saHist * _h_eta_p;      // Q+ tracks
  saHist * _h_eta_m;      // Q- tracks

  // -- eta for stopped hadrons --  
  saHist * _h_etaH;
  saHist * _h_etaH_p;      // Q+ tracks
  saHist * _h_etaH_m;      // Q- tracks

  // -- eta for W/Z analysis --  
  saHist * _h_etaW10;
  saHist * _h_etaW10_p;      // Q+ high pT>10 tracks
  saHist * _h_etaW10_m;      // Q- high pT>10 tracks

  saHist * _h_etaW15;
  saHist * _h_etaW15_p;      // Q+ high pT>15 tracks
  saHist * _h_etaW15_m;      // Q- high pT>15 tracks

  saHist * _h_etaW20;
  saHist * _h_etaW20_p;      // Q+ high pT>20 tracks
  saHist * _h_etaW20_m;      // Q- high pT>20 tracks



};

#endif /* SAMODULESNGMUONMING_H_ */
