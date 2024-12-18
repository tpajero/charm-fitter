/**
 * Charm Combination
 * Author: Tommaso Pajero, tommaso.pajero@cern.ch
 * Date: October 2021
 **/

#ifndef PDF_WS_h
#define PDF_WS_h

#include <TString.h>

#include <PDF_Abs.h>

#include "CharmUtils.h"

enum class WS_parametrisation { raxy, rrxy, ccprime };

class PDF_WS : public PDF_Abs {
 public:
  PDF_WS(TString measurement_id, const theory_config& th_cfg, WS_parametrisation p = WS_parametrisation::rrxy);
  PDF_WS(TString val, TString err, const theory_config& th_cfg);
  ~PDF_WS();
  void buildPdf() override;
  void initObservables(const TString& setName);
  virtual void initParameters() override;
  virtual void initRelations() override;
  void setCorrelations(TString measurement_id) override;
  void setObservables(TString measurement_id) override;
  void setUncertainties(TString measurement_id) override;

 private:
  void initRelationsCCPrime();
  void initRelationsRAXY();
  void initRelationsRRXY();
  void initRelationsXYM(RooArgList* theory);

  const theory_config th_cfg;
  const WS_parametrisation ws_param = WS_parametrisation::rrxy;
};

#endif
