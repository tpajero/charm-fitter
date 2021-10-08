/**
 * Charm Combination
 * Author: Tommaso Pajero, tommaso.pajero@cern.ch
 * Date: October 2021
 **/

#ifndef PDF_scan_DY_RS_h
#define PDF_scan_DY_RS_h

#include "PDF_Abs.h"
#include "ParametersCharmCombo.h"

using namespace RooFit;
using namespace std;

class PDF_scan_DY_RS : public PDF_Abs
{
    public:
        PDF_scan_DY_RS(const theory_config& th_cf);
        ~PDF_scan_DY_RS();
        void buildPdf();
        void initParameters(const theory_config& th_cf);
        void initObservables();
        virtual void initRelations(const theory_config& th_cf);
        void setCorrelations();
        void setObservables();
        void setUncertainties();
};

#endif

