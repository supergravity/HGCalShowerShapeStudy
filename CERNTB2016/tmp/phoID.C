////high pT photon ID
float PhotonIDHighPt(TreeReader &data, int i)
{

  // load necessary tree branches
  float* phoEt            = data.GetPtrFloat("phoEt");
  float* phoSCEta         = data.GetPtrFloat("phoSCEta");
  Int_t* phoEleVeto       = data.GetPtrInt("phoEleVeto");
  float* phoHoverE12      = data.GetPtrFloat("phoHoverE");
  float* phoSigmaIEtaIEta = data.GetPtrFloat("phoSigmaIEtaIEtaFull5x5");
  float* phoPFChIso       = data.GetPtrFloat("phoPFChIso");
  float* phoPFNeuIso      = data.GetPtrFloat("phoPFNeuIso");
  float* phoPFPhoIso      = data.GetPtrFloat("phoPFPhoIso");
  float  rho2012          = data.GetFloat("rho");

  // common ID cuts
  if (phoEleVeto[i] == 0) return 0;


  // to improve the code readability
  float absEta = fabs(phoSCEta[i]);

  // 0=ECAL barrel or 1=ECAL endcaps
  int iBE = (absEta < 1.4442) ? 0 : 1;

  // NOTE: the two-dimensional arrays below are [loose,medium,tight][EB,EE]

  ///H/E
  float HoverEcut[2] = {0.05, 0.05};
  
  if (phoHoverE12[i] > HoverEcut[iBE]) return 0;

  // sigmaIEtaIEta cut
  float sigmaIEIECut[2] = {0.0105, 0.028};
  if (phoSigmaIEtaIEta[i] > sigmaIEIECut[iBE]) return 0;

  // find eta bin which determines effective areas for the rho corrections to
  // isolations
  int bin;
  if (absEta < 0.9) bin = 0;
  else if (absEta >= .9 && absEta < 1.4442) bin = 1;
  else if (absEta >= 1.566 && absEta < 2.0) bin = 2;
  else if (absEta >= 2.0 && absEta < 2.2) bin = 3;
  else if (absEta >= 2.2 && absEta < 2.5) bin = 4;
  else bin = 4;
  
  // rho-corrected PF photon isolation
  double A[] = {0.17, 0.14, 0.11, 0.14, 0.22};
  double alpha = 2.5;
  double kappa[] = {4.53e-3, 4.53e-3,4.53e-3, 3e-3, 3e-3};

  double chIsoCut = 5;
  double phoIsoCut[2] = {2.75, 2.0};

  
  double corrphoIso = alpha + phoPFPhoIso[i] - rho2012 * A[bin] - kappa[bin] * phoEt[i];


  if ( corrphoIso > phoIsoCut[iBE]) return 0;
  if ( phoPFChIso[i] > chIsoCut ) return 0;
  
  
  // passed all ID and isolation cuts
  return 1;
}

