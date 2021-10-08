/**
 * Charm Combination
 * Author: Tommaso Pajero, tommaso.pajero@cern.ch
 * Date: October 2021
 **/

#ifndef PDF_yCP_minus_yCP_KP_h
#define PDF_yCP_minus_yCP_KP_h

#include "PDF_Abs.h"
#include "ParametersCharmCombo.h"

using namespace RooFit;
using namespace std;

class PDF_yCP_minus_yCP_KP : public PDF_Abs
{
    public:
        PDF_yCP_minus_yCP_KP(TString cObs, TString cErr, TString cCor,
                 const theory_config& th_cf);
        ~PDF_yCP_minus_yCP_KP();
        void buildPdf();
        void initObservables(const TString& setName);
        virtual void initParameters(const theory_config& th_cf);
        virtual void initRelations(const theory_config& th_cf);
        void setCorrelations(TString c);
        void setObservables(TString c);
        void setUncertainties(TString c);
};

#endif

