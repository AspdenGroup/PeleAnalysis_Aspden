
namespace amrex {

const char* buildInfoGetBuildDate() {

  static const char BUILD_DATE[] = "2021-03-26 20:22:12.082239";
  return BUILD_DATE;
}

const char* buildInfoGetBuildDir() {

  static const char BUILD_DIR[] = "/home/thomas/src/PeleAnalysis_Aspden/Src/ModelSpecificAnalysis";
  return BUILD_DIR;
}

const char* buildInfoGetBuildMachine() {

  static const char BUILD_MACHINE[] = "Linux 19S-STB-89575.campus.ncl.ac.uk 5.4.0-67-generic #75~18.04.1-Ubuntu SMP Tue Feb 23 19:17:50 UTC 2021 x86_64 x86_64 x86_64 GNU/Linux";
  return BUILD_MACHINE;
}

const char* buildInfoGetAMReXDir() {

  static const char AMREX_DIR[] = "../../amrex";
  return AMREX_DIR;
}

const char* buildInfoGetComp() {

  static const char COMP[] = "gcc";
  return COMP;
}

const char* buildInfoGetCompVersion() {

  static const char COMP_VERSION[] = "7.5.0";
  return COMP_VERSION;
}

// deprecated
const char* buildInfoGetFcomp() {

  static const char FCOMP[] = "";
  return FCOMP;
}

// deprecated
const char* buildInfoGetFcompVersion() {

  static const char FCOMP_VERSION[] = "";
  return FCOMP_VERSION;
}

const char* buildInfoGetCXXName() {

  static const char CXX_comp_name[] = "";
  return CXX_comp_name;
}

const char* buildInfoGetFName() {

  static const char F_comp_name[] = "";
  return F_comp_name;
}

const char* buildInfoGetCXXFlags() {

  static const char CXX_flags[] = "";
  return CXX_flags;
}

const char* buildInfoGetFFlags() {

  static const char F_flags[] = "";
  return F_flags;
}

const char* buildInfoGetLinkFlags() {

  static const char link_flags[] = "";
  return link_flags;
}

const char* buildInfoGetLibraries() {

  static const char libraries[] = "";
  return libraries;
}

const char* buildInfoGetAux(int i) {

  //static const char AUX1[] = "${AUX[1]}";

  static const char EMPT[] = "";

  switch(i)
  {

    default: return EMPT;
  }
}

int buildInfoGetNumModules() {
  // int num_modules = X;
  int num_modules = 4;

  return num_modules;
}

const char* buildInfoGetModuleName(int i) {

  //static const char MNAME1[] = "${MNAME[1]}";
  static const char AUX1[] = "EOS";
  static const char AUX2[] = "REACTIONS";
  static const char AUX3[] = "CHEMISTRY";
  static const char AUX4[] = "TRANSPORT";

  static const char EMPT[] = "";

  switch(i)
  {
    case 1: return AUX1;
    case 2: return AUX2;
    case 3: return AUX3;
    case 4: return AUX4;

    default: return EMPT;
  }
}

const char* buildInfoGetModuleVal(int i) {

  //static const char MVAL1[] = "${MVAL[1]}";
  static const char AUX1[] = "../../../PeleProduction_ISC2022/Submodules/PelePhysics/Eos/Fuego";
  static const char AUX2[] = "../../../PeleProduction_ISC2022/Submodules/PelePhysics/Reactions/Fuego";
  static const char AUX3[] = "BurkeDryer_mod_noArHe";
  static const char AUX4[] = "../../../PeleProduction_ISC2022/Submodules/PelePhysics/Transport";

  static const char EMPT[] = "";

  switch(i)
  {
    case 1: return AUX1;
    case 2: return AUX2;
    case 3: return AUX3;
    case 4: return AUX4;

    default: return EMPT;
  }
}

const char* buildInfoGetGitHash(int i) {

  //static const char HASH1[] = "${GIT[1]}";
  static const char HASH1[] = "5313462-dirty";
  static const char HASH2[] = "21.03-63-g9f297b9bc";
  static const char HASH3[] = "";

  static const char EMPT[] = "";

  switch(i)
  {
    case 1: return HASH1;
    case 2: return HASH2;
    case 3: return HASH3;

    default: return EMPT;
  }
}

const char* buildInfoGetBuildGitHash() {

  //static const char HASH[] = "${GIT}";
  static const char HASH[] = "";


  return HASH;
}

const char* buildInfoGetBuildGitName() {

  //static const char NAME[] = "";
  static const char NAME[] = "";


  return NAME;
}

}
