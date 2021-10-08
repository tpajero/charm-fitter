/**
 * Charm Combination
 * Author: Tommaso Pajero, tommaso.pajero@cern.ch
 * Date: October 2021
 **/

#ifndef PDF_DY_h
#define PDF_DY_h

#include "PDF_Abs.h"
#include "ParametersCharmCombo.h"

using namespace RooFit;
using namespace std;

class PDF_DY : public PDF_Abs
{
    public:
        PDF_DY(TString cObs, TString cErr, TString cCor,
                   const theory_config& th_cf);
        ~PDF_DY();
        void buildPdf();
        void initParameters(const theory_config& th_cf);
        void initObservables(const TString& setName);
        virtual void initRelations(const theory_config& th_cf);
        void setCorrelations(TString c);
        void setObservables(TString c);
        void setUncertainties(TString c);
};

#endif

