/**
 * Charm Combination
 * Author: Tommaso Pajero, tommaso.pajero@cern.ch
 * Date: October 2021
 **/

#ifndef PDF_BES_h
#define PDF_BES_h

#include "PDF_Abs.h"
#include "ParametersCharmCombo.h"

using namespace RooFit;
using namespace std;
using namespace Utils;

class PDF_BES : public PDF_Abs
{
    public:
        PDF_BES(TString cObs, TString cErr, TString cCor,
                const theory_config& th_cf);
        ~PDF_BES();
        void buildPdf();
        void initObservables(const TString& setName);
        virtual void initParameters(const theory_config& th_cf);
        virtual void initRelations(const theory_config& th_cf);
        void setCorrelations(TString c);
        void setObservables(TString c);
        void setUncertainties(TString c);
};

#endif

