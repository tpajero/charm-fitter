#ifndef CHARM_COMBINER
#define CHARM_COMBINER

#include <stdlib.h>
#include "GammaComboEngine.h"

#include "ParametersCharmCombo.h"
#include "PDF_AcpHH_LHCb_Run12.h"
#include "PDF_BES.h"
#include "PDF_BinFlip.h"
#include "PDF_Cleo.h"
#include "PDF_DY.h"
#include "PDF_DY_RS.h"
#include "PDF_Kpipi0.h"
#include "PDF_Kshh.h"
#include "PDF_RM.h"
#include "PDF_scan_DY_RS.h"
#include "PDF_XY.h"
#include "PDF_yCP.h"
#include "PDF_yCP_minus_yCP_RS.h"
#include "PDF_yCP_minus_yCP_KP.h"
#include "PDF_yCP_plus_yCP_RS.h"
#include "PDF_WS.h"
#include "PDF_WS_NoCPV.h"


void run_combiner(std::string name, theory_config th_cfg, int argc, char* argv[]);

#endif
