/**
 * Charm Combination
 * Author: Tommaso Pajero, tommaso.pajero@cern.ch
 * Date: October 2021
 **/

#ifndef PDF_BinFlip_h
#define PDF_BinFlip_h

#include "PDF_Abs.h"
#include "ParametersCharmCombo.h"

using namespace RooFit;
using namespace std;
using namespace Utils;

class PDF_BinFlip : public PDF_Abs
{
    public:
        PDF_BinFlip(TString cObs, TString cErr, TString cCor,
                    const theory_config& th_cf);
        ~PDF_BinFlip();
        void buildPdf();
        void initObservables(const TString& setName);
        virtual void initParameters(const theory_config& th_cf);
        virtual void initRelations(const theory_config& th_cf);
        void setCorrelations(TString c);
        void setObservables(TString c);
        void setUncertainties(TString c);
};

#endif

