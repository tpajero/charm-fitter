/**
 * Charm Combination
 * Author: Tommaso Pajero, tommaso.pajero@cern.ch
 * Date: October 2021
 **/

#ifndef PDF_WS_h
#define PDF_WS_h

#include "CharmUtils.h"

#include <PDF_Charm.h>

#include <TString.h>

#include <set>
#include <string>

class RooArgList;

namespace parametrisations {
  enum class kpi { raxy, rrxy, ccprime };
}  // namespace parametrisations

class PDF_WS : public PDF_Charm {
 public:
  PDF_WS(TString measurement_id, parametrisations::mix mix_param,
         parametrisations::kpi p = parametrisations::kpi::rrxy);
  PDF_WS(TString val, TString err, parametrisations::mix mix_param);
  void buildPdf() override;
  void initObservables() override;
  void initRelations() override;
  void setCorrelations(TString measurement_id) override;
  void setObservables(TString measurement_id) override;
  void setUncertainties(TString measurement_id) override;

 private:
  void initRelationsCCPrime();
  void initRelationsRAXY();
  void initRelationsRRXY();
  void initRelationsXYM(RooArgList* theory);
  std::set<std::string> getParameterNames() const override;

  const parametrisations::mix mix_param;
  const parametrisations::kpi ws_param = parametrisations::kpi::rrxy;
  TString label;
};

#endif
