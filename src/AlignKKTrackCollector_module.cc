// Ryunosuke O'Neil, 2019
// roneil@fnal.gov
// ryunoneil@gmail.com

// A module to collect Cosmic NoField tracks and write out 'Mille' data files used as input to a
// Millepede-II alignment fit.

// Consult README.md for more information
#include <algorithm>
#include <cmath>
#include <cstddef>
#include <exception>
#include <fstream>
#include <iostream>
#include <memory>
#include <ostream>
#include <stdint.h>
#include <string>
#include <unordered_map>
#include <utility>
#include <vector>

#include "Offline/DbTables/inc/TrkAlignElement.hh"
#include "Offline/DbTables/inc/TrkAlignStraw.hh"
#include "Offline/GeneralUtilities/inc/BitMap.hh"
#include "Offline/GeneralUtilities/inc/HepTransform.hh"
#include "Offline/GeometryService/inc/GeomHandle.hh"
#include "Minuit2/MnUserCovariance.h"
#include "Offline/Mu2eUtilities/inc/TwoLinePCA.hh"

#include "boost/math/distributions/chi_squared.hpp"
#include "boost/math/distributions/normal.hpp"

#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Core/ProducerTable.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "art_root_io/TFileService.h"

#include "Offline/TrackerConditions/inc/StrawResponse.hh"
#include "Offline/TrackerConditions/inc/TrackerStatus.hh"
#include "Offline/TrackerGeom/inc/Panel.hh"
#include "Offline/TrackerGeom/inc/Plane.hh"
#include "Offline/TrackerGeom/inc/Straw.hh"
#include "Offline/TrackerGeom/inc/Tracker.hh"

#include "Offline/DbService/inc/DbHandle.hh"
#include "Offline/ProditionsService/inc/ProditionsHandle.hh"

#include "Offline/RecoDataProducts/inc/ComboHit.hh"
#include "Offline/RecoDataProducts/inc/KalSeed.hh"
#include "Offline/RecoDataProducts/inc/TrkFitFlag.hh"

#include "Offline/DataProducts/inc/StrawId.hh"

#include "TAxis.h"
#include "TH1F.h"
#include "TMatrixDSym.h"
#include "TTree.h"

#include "canvas/Utilities/Exception.h"
#include "canvas/Utilities/InputTag.h"

#include "CLHEP/Vector/ThreeVector.h"

#include "cetlib_except/exception.h"

#include "fhiclcpp/types/Atom.h"
#include "fhiclcpp/types/Comment.h"
#include "fhiclcpp/types/Name.h"
#include "fhiclcpp/types/Table.h"
#include "fhiclcpp/types/detail/validationException.h"

#include "TrackerAlignment/inc/MilleDataWriter.hh"
#include "TrackerAlignment/inc/AlignmentUtilities.hh"



namespace mu2e {

class AlignKKTrackCollector : public art::EDAnalyzer {
public:
  size_t _dof_per_plane = 6; // dx, dy, dz, a, b, g (translation, rotation)
  size_t _dof_per_panel = 6; // dx, dy, dz, a, b, g (translation, rotation)

  struct Config {
    using Name = fhicl::Name;
    using Comment = fhicl::Comment;

    fhicl::Atom<int> diaglvl{Name("diagLevel"), Comment("diagnostic level")};

    fhicl::Atom<art::InputTag> kktag{    Name("KalSeedCollection"), Comment("tag for cosmic track seed collection")};
    fhicl::Atom<bool> millefilegzip{      Name("GzipCompression"), Comment("Enable gzip compression for millepede output file")};
    fhicl::Atom<std::string> millefile{   Name("MilleFile"), Comment("Output filename for Millepede track data file")};
    fhicl::Atom<std::string> steerfile{   Name("SteerFile"), Comment("Output filename for Millepede steering file")};
    fhicl::Atom<std::string> paramfile{   Name("ParamFile"), Comment("Output filename for Millepede parameters file")};
    fhicl::Atom<std::string> extrafile{   Name("ExtraFile"), Comment("Output filename for Millepede parameters file")};
    fhicl::Atom<std::string> constrfile{  Name("ConstrFile"), Comment("Output filename for Millepede constraints file")};
    fhicl::Sequence<std::string> mpsteers{Name("SteeringOpts"), Comment("Additional configuration options to add to generated steering file.")};

    fhicl::Sequence<int> fixplane{        Name("FixPlane"), Comment("Planes to fix. The parameters are fixed to the proditions value.")};
    fhicl::Sequence<int> fixpanel{        Name("FixPanel"), Comment("Panels to fix. The parameters are fixed to the proditions value.")};
    fhicl::Atom<bool> enableplanetranslation{ Name("EnablePlaneTranslationDOF"), Comment("Fit for plane translations")};
    fhicl::Atom<bool> enableplanerotation{ Name("EnablePlaneRotationDOF"), Comment("Fit for plane rotations")};
    fhicl::Atom<bool> enablepaneltranslation{ Name("EnablePanelTranslationDOF"), Comment("Fit for panel translations")};
    fhicl::Atom<bool> enablepanelrotation{ Name("EnablePanelRotationDOF"), Comment("Fit for panel rotations")};
    fhicl::Atom<std::string> weakconstraints{ Name("WeakConstraints"), Comment("Which weak constraint strategy to use. Either 'None', 'Fix', or 'Measurement'")};
    fhicl::Sequence<float> weakvalues{Name("WeakValues"), Comment("Constrain weak modes at these values"), std::vector<float>{0,0,0,0,0,0}};
    fhicl::Sequence<float> weaksigmas{Name("WeakSigmas"), Comment("Constrain weak modes with these sigmas"), std::vector<float>{0,0,0,0,0,0}};
    fhicl::Atom<bool> fixpanelperplane{ Name("FixPanelPerPlane"), Comment("Fix single panel per plane to remove redundant DOF")};
    fhicl::Atom<bool> panelconstraints{       Name("PanelConstraints"), Comment("Whether to constrain each panel DOF")};
    fhicl::Sequence<float> panelvaluesx{Name("PanelValuesX"), Comment("Constrain panel x at these values"), std::vector<float>{}};
    fhicl::Sequence<float> panelvaluesy{Name("PanelValuesY"), Comment("Constrain panel y at these values"), std::vector<float>{}};
    fhicl::Sequence<float> panelvaluesz{Name("PanelValuesZ"), Comment("Constrain panel z at these values"), std::vector<float>{}};
    fhicl::Sequence<float> panelvaluesrx{Name("PanelValuesRX"), Comment("Constrain panel rx at these values"), std::vector<float>{}};
    fhicl::Sequence<float> panelvaluesry{Name("PanelValuesRY"), Comment("Constrain panel ry at these values"), std::vector<float>{}};
    fhicl::Sequence<float> panelvaluesrz{Name("PanelValuesRZ"), Comment("Constrain panel rz at these values"), std::vector<float>{}};
    fhicl::Atom<float> panelsigmax{Name("PanelSigmaX"), Comment("Constrain panel x with this sigma"),0};
    fhicl::Atom<float> panelsigmay{Name("PanelSigmaY"), Comment("Constrain panel y with this sigma"),0};
    fhicl::Atom<float> panelsigmaz{Name("PanelSigmaZ"), Comment("Constrain panel z with this sigma"),0};
    fhicl::Atom<float> panelsigmarx{Name("PanelSigmaRX"), Comment("Constrain panel rx with this sigma"),0};
    fhicl::Atom<float> panelsigmary{Name("PanelSigmaRY"), Comment("Constrain panel ry with this sigma"),0};
    fhicl::Atom<float> panelsigmarz{Name("PanelSigmaRZ"), Comment("Constrain panel rz with this sigma"),0};

    fhicl::Atom<float> maxresid{ Name("MaxResid"), Comment("Max resid")};
    fhicl::Atom<int> minhits{ Name("MinHits"), Comment("Min hits")};
    fhicl::Atom<int> maxinactive{ Name("MaxInactive"), Comment("Max inactive")};
    fhicl::Atom<int> minplanes{ Name("MinPlanes"), Comment("Min planes")};
    fhicl::Atom<int> minpanels{ Name("MinPanels"), Comment("Min panels")};

  };

  typedef art::EDAnalyzer::Table<Config> Parameters;

  explicit AlignKKTrackCollector(const Parameters& settings);
  virtual ~AlignKKTrackCollector() {}

  void beginJob();
  void endJob();
  void beginRun(art::Run const&);
  void analyze(art::Event const&);

private:

  int _diag;
  art::ProductToken<KalSeedCollection> _kktag;
  bool _gzipCompress;
  std::string _milleFilename;
  std::string _steerFilename;
  std::string _paramFilename;
  std::string _extraFilename;
  std::string _constrFilename;
  std::vector<std::string> _steerLines;

  std::vector<bool> _fixedPlanes;
  std::vector<bool> _fixedPanels;
  bool _enablePlaneTranslationDOF, _enablePlaneRotationDOF, _enablePanelTranslationDOF, _enablePanelRotationDOF;
  std::string _weakConstraints;
  std::vector<float> _weakValues, _weakSigmas;
  bool _fixPanelPerPlane;
  bool _panelConstraints;
  std::vector<std::vector<float>> _panelValues;
  std::vector<float> _panelSigmas;

  float _maxresid;
  size_t _minhits, _maxinactive, _minplanes, _minpanels;


  bool _wroteMillepedeParams;
  MilleDataWriter<double> _milleFile;
  std::map<int, int> _DOFCounts;
  std::vector<int> _planeCounts, _panelCounts;
  std::vector<bool> _panelActive;
  int _tracksWritten;
  std::vector<float> _startingAlignPlanes; 
  std::vector<float> _startingAlignPanels; 

  const KalSeedCollection* _kscol;
  const Tracker* _nominalTracker;

  ProditionsHandle<TrackerStatus> _trackerStatus_h;

  std::unique_ptr<DbHandle<TrkAlignPlane>> _trkAlignPlane_h;
  std::unique_ptr<DbHandle<TrkAlignPanel>> _trkAlignPanel_h;

  int getLabel(int const&, int const&, int const&);
  std::vector<int> generateDOFLabels(uint16_t plane, uint16_t panel);
  std::vector<int> generateDOFLabels(StrawId const& strw);
  std::vector<double> pruneInactiveGlobals(uint16_t plane, uint16_t panel, std::vector<double> &derivativesGlobal);
  bool isDOFenabled(int object_class, int object_id, int dof_n);
  bool isDOFfixed(int object_class, int object_id, int dof_n);
  void cacheStartingParams(TrkAlignPlane const& alignConstPlanes, TrkAlignPanel const& alignConstPanels, TrackerStatus const& trackerStatus);
  void writeMillepedeConstraints(Tracker const& tracker);
  void writeMillepedeSteering();
  void writeMillepedeParams();
  void writeMillepedeExtras();

  bool goodTrack(KalSeed const& sts);
};

AlignKKTrackCollector::AlignKKTrackCollector(const Parameters& conf) :
  art::EDAnalyzer(conf), 
  _diag(conf().diaglvl()),
  _kktag(consumes<KalSeedCollection>(conf().kktag())),
  _gzipCompress(conf().millefilegzip()), 
  _milleFilename(conf().millefile()),
  _steerFilename(conf().steerfile()), 
  _paramFilename(conf().paramfile()),
  _extraFilename(conf().extrafile()),
  _constrFilename(conf().constrfile()),
  _steerLines(conf().mpsteers()),
  _enablePlaneTranslationDOF(conf().enableplanetranslation()),
  _enablePlaneRotationDOF(conf().enableplanerotation()),
  _enablePanelTranslationDOF(conf().enablepaneltranslation()),
  _enablePanelRotationDOF(conf().enablepanelrotation()),
  _weakConstraints(conf().weakconstraints()),
  _weakValues(conf().weakvalues()),
  _weakSigmas(conf().weaksigmas()),
  _fixPanelPerPlane(conf().fixpanelperplane()),
  _panelConstraints(conf().panelconstraints()),
  _maxresid(conf().maxresid()),
  _minhits(conf().minhits()),
  _maxinactive(conf().maxinactive()),
  _minplanes(conf().minplanes()),
  _minpanels(conf().minpanels()),
  _wroteMillepedeParams(false), 
  _milleFile(_milleFilename, _gzipCompress),
  _tracksWritten(0)
{
  _fixedPlanes = std::vector<bool>(StrawId::_nplanes,false);
  for (size_t i=0;i<conf().fixplane().size();i++)
    _fixedPlanes[conf().fixplane()[i]] = true;
  _fixedPanels = std::vector<bool>(StrawId::_nupanels,false);
  for (size_t i=0;i<conf().fixpanel().size();i++)
    _fixedPanels[conf().fixpanel()[i]] = true;

  for (size_t p=0;p<StrawId::_nplanes;p++){
    for (size_t j=0;j<_dof_per_plane;j++)
      _DOFCounts[getLabel(1,p,j)] = 0;
  }
  for (size_t p=0;p<StrawId::_nupanels;p++){
    for (size_t j=0;j<_dof_per_panel;j++)
      _DOFCounts[getLabel(2,p,j)] = 0;
  }
  _planeCounts = std::vector<int>(StrawId::_nplanes,0);
  _panelCounts = std::vector<int>(StrawId::_nupanels,0);

  if (_panelConstraints){
    _panelValues.push_back(conf().panelvaluesx());
    _panelValues.push_back(conf().panelvaluesy());
    _panelValues.push_back(conf().panelvaluesz());
    _panelValues.push_back(conf().panelvaluesrx());
    _panelValues.push_back(conf().panelvaluesry());
    _panelValues.push_back(conf().panelvaluesrz());
    _panelSigmas.push_back(conf().panelsigmax());
    _panelSigmas.push_back(conf().panelsigmay());
    _panelSigmas.push_back(conf().panelsigmaz());
    _panelSigmas.push_back(conf().panelsigmarx());
    _panelSigmas.push_back(conf().panelsigmary());
    _panelSigmas.push_back(conf().panelsigmarz());
    for (size_t i=0;i<6;i++){
      if (_panelValues[i].size() != StrawId::_nupanels){
        if (_panelValues[i].size() == 0){
          std::cout << "PanelValues " << i << " not set! using 0 as default" << std::endl;
        } else{
          std::cout << "Warning! Incorrect number of panelValues " << i << " (" << _panelValues[i].size() << " != " << StrawId::_nupanels << ")" << std::endl;
        }
        _panelValues[i] = std::vector<float>(StrawId::_nupanels,0);
      }
    }
  }
}



void AlignKKTrackCollector::beginJob() {
  _trkAlignPlane_h = std::make_unique<DbHandle<TrkAlignPlane>>();
  _trkAlignPanel_h = std::make_unique<DbHandle<TrkAlignPanel>>();
}

void AlignKKTrackCollector::endJob() {
  if (_diag > 0) {
    std::cout << "AlignKKTrackCollector: wrote " 
              << _tracksWritten << " tracks to " 
              << _milleFilename
              << std::endl;
  }

  writeMillepedeConstraints(*_nominalTracker);
  writeMillepedeParams();
  writeMillepedeSteering();
  writeMillepedeExtras();
}

void AlignKKTrackCollector::beginRun(art::Run const& run) {
  GeomHandle<Tracker> track;
  _nominalTracker = track.get();
}

void AlignKKTrackCollector::analyze(art::Event const& event) {
  TrackerStatus const& trackerStatus = _trackerStatus_h.get(event.id());

  auto alignConsts_planes = _trkAlignPlane_h->get(event.id());
  auto alignConsts_panels = _trkAlignPanel_h->get(event.id());

  if (!_wroteMillepedeParams) {
    cacheStartingParams(alignConsts_planes, alignConsts_panels, trackerStatus);
    _wroteMillepedeParams = true;
  }

  auto kkH = event.getValidHandle<KalSeedCollection>(_kktag);

  if (kkH.product() == 0) {
    return;
  }

  KalSeedCollection const& kkcol = *kkH.product();

  for (KalSeed const& ks : kkcol) {
    if (!goodTrack(ks))
      continue;

    std::vector<double> residuals;
    std::vector<double> residual_errs;
    std::vector<std::vector<double>> global_derivs_temp;
    std::vector<std::vector<double>> local_derivs_temp;
    std::vector<std::vector<int>> labels_temp;

    size_t nhits = ks.hits().size();
    for (size_t i=0;i<nhits;i++){
      auto const& hit = ks.hits()[i];
      auto const& hitcalib = ks.hitCalibInfos()[i];

      if (hit.strawHitState() <=  WireHitState::inactive)
        continue;

      StrawId const& straw_id = hit._sid;
      auto plane_id = straw_id.getPlane();
      auto panel_id = straw_id.uniquePanel();

      std::vector<int> labels = generateDOFLabels(straw_id);
      double residual = hit._rdresid;
      double residual_error = sqrt(hit._rdresidmvar);
      std::vector<double> global_derivs;
      std::vector<double> local_derivs;
      for (size_t j=0;j<5;j++)
        local_derivs.push_back(hitcalib._dDdP[j]);
      for (size_t j=0;j<6;j++)
        global_derivs.push_back(hitcalib._dDdPlaneAlign[j]);
      for (size_t j=0;j<6;j++)
        global_derivs.push_back(hitcalib._dDdPanelAlign[j]);
      std::vector<double> global_derivs_pruned = pruneInactiveGlobals(plane_id, panel_id, global_derivs);
      _milleFile.pushHit(local_derivs, global_derivs_pruned, labels, residual, residual_error);

      for (size_t j=0;j<labels.size();j++)
        _DOFCounts[labels[j]]++;
      _planeCounts[plane_id]++;
      _panelCounts[panel_id]++;
    }
    
    // Write the track buffer to file
    _milleFile.flushTrack();

    _tracksWritten++;

    if (_diag > 1) {
      std::cout << "wrote track " << _tracksWritten << std::endl;
    }
  }
}

bool AlignKKTrackCollector::goodTrack(KalSeed const& ks){
  if (!ks.status().hasAnyProperty(TrkFitFlag::kalmanConverged))
    return false;

  std::set<unsigned> planes;
  std::set<unsigned> panels;
  size_t nhits = 0;
  size_t nactive = 0;

  for (auto ihit = ks.hits().begin(); ihit != ks.hits().end(); ++ihit) {
    ++nhits;
    if (ihit->strawHitState() > WireHitState::inactive)
      ++nactive;
    if (fabs(ihit->_rdresid) > _maxresid)
      return false;

    planes.insert(ihit->strawId().plane());
    panels.insert(ihit->strawId().uniquePanel());
  }
  size_t nplanes = planes.size();
  size_t npanels = panels.size();

  if (nactive < _minhits)
    return false;
  if (nhits-nactive > _maxinactive)
    return false;
  if (nplanes < _minplanes || npanels < _minpanels)
    return false;

  return true;
}


int AlignKKTrackCollector::getLabel(int const& object_cls, int const& obj_uid, int const& dof_id) {
  // object class: 0 - 9 - i.e. 1 for planes, 2 for panels
  // object unique id: 0 - 999 supports up to 999 unique objects which is fine for this level of
  // alignment
  // object dof id: 0 - 9

  // (future) a straw label might be:
  // id  plane panel straw dof
  // 3   00    0     00    0
  // min
  // 3 00 0 00 0
  // max
  // 3,356,959

  // 1 000 0

  return object_cls * 10000 + obj_uid * 10 + dof_id;
}

std::vector<int> AlignKKTrackCollector::generateDOFLabels(StrawId const& strw) {
  return generateDOFLabels(strw.getPlane(), strw.uniquePanel());
}

std::vector<int> AlignKKTrackCollector::generateDOFLabels(uint16_t plane, uint16_t panel) {
  std::vector<int> labels;
  labels.reserve(_dof_per_plane + _dof_per_panel);

  for (size_t dof_n = 0; dof_n < _dof_per_plane; dof_n++) {
    if (!isDOFenabled(1, plane, dof_n)) {
      continue;
    }
    labels.push_back(getLabel(1, plane, dof_n));
  }
  for (size_t dof_n = 0; dof_n < _dof_per_panel; dof_n++) {
    if (!isDOFenabled(2, panel, dof_n)) {
      continue;
    }
    labels.push_back(getLabel(2, panel, dof_n));
  }
  
  return labels;
}

std::vector<double> AlignKKTrackCollector::pruneInactiveGlobals(uint16_t plane, uint16_t panel, std::vector<double> &derivativesGlobal) {
  std::vector<double> prunedDerivativesGlobal;
  size_t index = 0;
  for (size_t dof_n = 0; dof_n < _dof_per_plane; dof_n++) {
    if (isDOFenabled(1, plane, dof_n)) {
      prunedDerivativesGlobal.push_back(derivativesGlobal[index]);
    }
    index++;
  }
  for (size_t dof_n = 0; dof_n < _dof_per_panel; dof_n++) {
    if (isDOFenabled(2, panel, dof_n)) {
      prunedDerivativesGlobal.push_back(derivativesGlobal[index]);
    }
    index++;
  }
  
  return prunedDerivativesGlobal;
}

bool AlignKKTrackCollector::isDOFenabled(int object_class, int object_id, int dof_n) {
  if (object_class == 2 && !_enablePanelTranslationDOF && dof_n <= 2) {
    return false;
  }
  if (object_class == 2 && !_enablePanelRotationDOF && dof_n > 2) {
    return false;
  }
  // we aren't sensitive to panel translations along u
  if (object_class == 2 && dof_n == 0)
    return false;
  if (object_class == 1 && !_enablePlaneTranslationDOF && dof_n <= 2) {
    return false;
  }
  if (object_class == 1 && !_enablePlaneRotationDOF && dof_n > 2) {
    return false;
  }
  return true;
}

bool AlignKKTrackCollector::isDOFfixed(int object_class, int object_id, int dof_n) {
  if (object_class == 1) 
    return _fixedPlanes[object_id];
  if (object_class == 2) 
    return _fixedPanels[object_id];
  return false;
}

void AlignKKTrackCollector::writeMillepedeConstraints(Tracker const& nominalTracker){
  std::ofstream output_file(_constrFilename);
  output_file << "! Generated by AlignKKTrackCollector" << std::endl;

  // first check if any DOF are not supposed to be active but had some hits and fix those
  for (uint16_t p = 0; p < StrawId::_nplanes; ++p) {
    if (_fixedPlanes[p])
      continue;
    bool hasactive = false;
    for (size_t pa=0;pa<6;pa++){
      if (_panelActive[p*6+pa])
        hasactive = true;
    }
    if (!hasactive || _planeCounts[p] == 0)
      _fixedPlanes[p] = true;
    if (!hasactive && _planeCounts[p] > 0){
      for (size_t dof_n=0;dof_n<_dof_per_plane;dof_n++){
        if (isDOFenabled(1, p, dof_n) && !isDOFfixed(1, p, dof_n)){
          output_file << "Constraint   " << 0 << std::endl;
          output_file << getLabel(1, p, dof_n) << "  1" << std::endl;
        }
      }
    }
  }
  for (uint16_t p=0;p<StrawId::_nupanels;++p) {
    if (_fixedPanels[p])
      continue;
    if (_panelCounts[p] == 0 || !_panelActive[p])
      _fixedPanels[p] = true;
    if (_panelCounts[p] > 0 && !_panelActive[p]){
      for (size_t dof_n=0;dof_n<_dof_per_panel;dof_n++){
        if (isDOFenabled(2, p, dof_n) && !isDOFfixed(2, p, dof_n)){
          output_file << "Constraint   " << 0 << std::endl;
          output_file << getLabel(2, p, dof_n) << "  1" << std::endl;
        }
      }
    }
  }

  // plane overall translations and rotations are fixed
  output_file << "! planes" << std::endl;
  for (size_t dof_n = 0; dof_n < _dof_per_plane; dof_n++) {

    // x and y parallelogramming not constrained here
    if (dof_n == 3 || dof_n == 4)
      continue;
    // check if any of these DOF are enabled
    bool has_enabled = false;
    for (uint16_t p=0;p<StrawId::_nplanes;++p){
      if (isDOFenabled(1, p, dof_n) && !isDOFfixed(1, p, dof_n)) {
        has_enabled = true;
        break;
      }
    }
    if (!has_enabled)
      continue;
    // measure the initial overall translation or rotation to check that it is zero
    double current_overall = 0;
    for (uint16_t p=0;p<StrawId::_nplanes;p++){
      if (nominalTracker.getPlane(p).planeToDS().rotation().getTheta() == 0)
        current_overall += _startingAlignPlanes[p*6+dof_n];
      else
        current_overall -= _startingAlignPlanes[p*6+dof_n];
    }
    if (fabs(current_overall) > 1e-5){
      std::cout << "AlignKKTrackCollector: Warning - planes have nonzero total misalign for dof " << dof_n << " : " << current_overall << std::endl;
    }
    output_file << "Constraint   0" << std::endl;
    for (uint16_t p = 0; p < StrawId::_nplanes; ++p) {
      if (!isDOFenabled(1, p, dof_n))
        continue;
      if (dof_n == 0 || dof_n == 2 || dof_n == 3 || dof_n == 5){
        // if this plane is rotated, need to invert constraint
        if (nominalTracker.getPlane(p).planeToDS().rotation().getTheta() == 0)
          output_file << getLabel(1, p, dof_n) << "    1" << std::endl;
        else
          output_file << getLabel(1, p, dof_n) << "    -1" << std::endl;
      }else{
        output_file << getLabel(1, p, dof_n) << "    1" << std::endl;
      }
    }
  }
  
  if (_panelConstraints){
    for (uint16_t p=0; p< StrawId::_nupanels; ++p){
      for (size_t dof_n = 0; dof_n < _dof_per_panel; dof_n++) {
        if (!isDOFenabled(2, p, dof_n) || isDOFfixed(2, p, dof_n)) {
          continue;
        }
        output_file << "Measurement " << _panelValues[dof_n][p]-_startingAlignPanels[p*6+dof_n] << " " << _panelSigmas[dof_n] << std::endl;
        output_file << getLabel(2, p, dof_n) << "  1" << std::endl;
      }
    }
  }

  if (_fixPanelPerPlane){
    if (_enablePanelTranslationDOF || _enablePanelRotationDOF){
      output_file << "! fixing one panel dof per plane" << std::endl;
      for (size_t pl=0;pl<StrawId::_nplanes;pl++){
        bool fixed = false;
        bool active = false;
        for (size_t i=0;i<6;i++){
          if (_panelActive[pl*StrawId::_npanels+i] && _panelCounts[pl*StrawId::_npanels+i])
            active = true;
        }
        if (!active)
          continue;
        // find two active panels at 90 degrees
        for (size_t i=0;i<3;i++){
          if (_panelCounts[pl*StrawId::_npanels+i*2] && _panelCounts[pl*StrawId::_npanels+i*2+1] &&
              _panelActive[pl*StrawId::_npanels+i*2] && _panelActive[pl*StrawId::_npanels+i*2+1] &&
              !_fixedPanels[pl*StrawId::_npanels+i*2] && !_fixedPanels[pl*StrawId::_npanels+i*2+1])
          {
            // fix one full panel and v-translation of one panel at 90 degrees
            for (size_t dof_n=0;dof_n<_dof_per_panel;dof_n++){
              if (isDOFenabled(2, pl*StrawId::_npanels+i*2, dof_n) && !isDOFfixed(2, pl*StrawId::_npanels+i*2, dof_n)){
                output_file << "Constraint   " << -1*_startingAlignPanels[(pl*StrawId::_npanels+i*2)*6 + dof_n] << std::endl;
                output_file << getLabel(2,pl*StrawId::_npanels+i*2,dof_n) << "  1" << std::endl;
              }
            }
            if (isDOFenabled(2, pl*StrawId::_npanels+i*2+1, 1) && !isDOFfixed(2, pl*StrawId::_npanels+i*2+1, 1)){
              output_file << "Constraint   " << -1*_startingAlignPanels[(pl*StrawId::_npanels+i*2+1)*6 + 1] << std::endl;
              output_file << getLabel(2,pl*StrawId::_npanels+i*2+1,1) << "  1" << std::endl;
            }
            fixed = true;
            break;
          }
        }
        if (!fixed){
          std::cout << "AlignTrackCollector: Warning - unable to constrain redundant panel dofs in plane " << pl << std::endl;
        }
      }
    }
  }else{
    for (size_t pl=0;pl<StrawId::_nplanes;pl++){
      bool active = false;
      for (size_t i=0;i<6;i++){
        if (isDOFenabled(2, pl*StrawId::_npanels+i, 2) && !isDOFfixed(2, pl*StrawId::_npanels+i, 2)){
          active = true;
        }
      }
      if (active){
        output_file << "Constraint   0" << std::endl;
        for (size_t i=0;i<6;i++){
          if (isDOFenabled(2, pl*StrawId::_npanels+i, 2) && !isDOFfixed(2, pl*StrawId::_npanels+i, 2)){
            if ((i%2) == 0)
              output_file << getLabel(2,pl*StrawId::_npanels+i,2) << "  1" << std::endl;
            else
              output_file << getLabel(2,pl*StrawId::_npanels+i,2) << "  -1" << std::endl;
          }
        }
      }
      active = false;
      for (size_t i=0;i<6;i++){
        if (isDOFenabled(2, pl*StrawId::_npanels+i, 1) && !isDOFfixed(2, pl*StrawId::_npanels+i, 1)){
          active = true;
        }
      }
      if (active){
        output_file << "Constraint   0" << std::endl;
        for (size_t i=0;i<6;i++){
          if (isDOFenabled(2, pl*StrawId::_npanels+i, 1) && !isDOFfixed(2, pl*StrawId::_npanels+i, 1)){
            output_file << getLabel(2,pl*StrawId::_npanels+i,1) << "  " << nominalTracker.getPlane(pl).getPanel(i).vDirection().dot(CLHEP::Hep3Vector(1,0,0)) << std::endl;
          }
        }
        output_file << "Constraint   0" << std::endl;
        for (size_t i=0;i<6;i++){
          if (isDOFenabled(2, pl*StrawId::_npanels+i, 1) && !isDOFfixed(2, pl*StrawId::_npanels+i, 1)){
            output_file << getLabel(2,pl*StrawId::_npanels+i,1) << "  " << nominalTracker.getPlane(pl).getPanel(i).vDirection().dot(CLHEP::Hep3Vector(0,1,0)) << std::endl;
          }
        }
      }
      active = false;
      for (size_t i=0;i<6;i++){
        if (isDOFenabled(2, pl*StrawId::_npanels+i, 1) && !isDOFfixed(2, pl*StrawId::_npanels+i, 1)){
          active = true;
        }
      }
      if (active){
        output_file << "Constraint   0" << std::endl;
        for (size_t i=0;i<6;i++){
          if (isDOFenabled(2, pl*StrawId::_npanels+i, 3) && !isDOFfixed(2, pl*StrawId::_npanels+i, 3)){
            double dpu = nominalTracker.getPlane(pl).getPanel(i).uDirection().dot(nominalTracker.getPlane(pl).uDirection());
            double dpv = nominalTracker.getPlane(pl).getPanel(i).vDirection().dot(nominalTracker.getPlane(pl).uDirection());
            double dpw = nominalTracker.getPlane(pl).getPanel(i).wDirection().dot(nominalTracker.getPlane(pl).uDirection());
            output_file << getLabel(2,pl*StrawId::_npanels+i,3) << "  " << dpu << std::endl;
            output_file << getLabel(2,pl*StrawId::_npanels+i,4) << "  " << dpv << std::endl;
            output_file << getLabel(2,pl*StrawId::_npanels+i,5) << "  " << dpw << std::endl;
          }
        }
        output_file << "Constraint   0" << std::endl;
        for (size_t i=0;i<6;i++){
          if (isDOFenabled(2, pl*StrawId::_npanels+i, 3) && !isDOFfixed(2, pl*StrawId::_npanels+i, 3)){
            double dpu = nominalTracker.getPlane(pl).getPanel(i).uDirection().dot(nominalTracker.getPlane(pl).vDirection());
            double dpv = nominalTracker.getPlane(pl).getPanel(i).vDirection().dot(nominalTracker.getPlane(pl).vDirection());
            double dpw = nominalTracker.getPlane(pl).getPanel(i).wDirection().dot(nominalTracker.getPlane(pl).vDirection());
            output_file << getLabel(2,pl*StrawId::_npanels+i,3) << "  " << dpu << std::endl;
            output_file << getLabel(2,pl*StrawId::_npanels+i,4) << "  " << dpv << std::endl;
            output_file << getLabel(2,pl*StrawId::_npanels+i,5) << "  " << dpw << std::endl;
          }
        }
        output_file << "Constraint   0" << std::endl;
        for (size_t i=0;i<6;i++){
          if (isDOFenabled(2, pl*StrawId::_npanels+i, 3) && !isDOFfixed(2, pl*StrawId::_npanels+i, 3)){
            double dpu = nominalTracker.getPlane(pl).getPanel(i).uDirection().dot(nominalTracker.getPlane(pl).wDirection());
            double dpv = nominalTracker.getPlane(pl).getPanel(i).vDirection().dot(nominalTracker.getPlane(pl).wDirection());
            double dpw = nominalTracker.getPlane(pl).getPanel(i).wDirection().dot(nominalTracker.getPlane(pl).wDirection());
            output_file << getLabel(2,pl*StrawId::_npanels+i,3) << "  " << dpu << std::endl;
            output_file << getLabel(2,pl*StrawId::_npanels+i,4) << "  " << dpv << std::endl;
            output_file << getLabel(2,pl*StrawId::_npanels+i,5) << "  " << dpw << std::endl;
          }
        }

      }
    }
  }

  output_file << "! weak modes follow" << std::endl;
  if (_weakConstraints == "None"){
    return;
    // if fixing specific panels then no further overall constraints
  }else{

    // calculate current values
    double xskew = 0;
    double yskew = 0;
    double xparallel = 0;
    double yparallel = 0;
    double zsqueeze = 0;
    double ztwist = 0;
    double xskew_constraint = 0;
    double yskew_constraint = 0;
    double xparallel_constraint = 0;
    double yparallel_constraint = 0;
    double zsqueeze_constraint = 0;
    double ztwist_constraint = 0;
    for (uint16_t p = 0; p < StrawId::_nplanes; ++p) {
      if (nominalTracker.getPlane(p).planeToDS().rotation().getTheta() == 0){
        if (isDOFenabled(1,p,0) && !isDOFfixed(1,p,0)){
          xskew += _startingAlignPlanes[p*6+0] * (p-17.5);
          xskew_constraint += _weakValues[0] * (p/35. - 0.5) * (p-17.5);
        }
        if (isDOFenabled(1,p,1) && !isDOFfixed(1,p,1)){
          yskew += _startingAlignPlanes[p*6+1] * (p-17.5);
          yskew_constraint += _weakValues[1] * (p/35. - 0.5) * (p-17.5);
        }
        if (isDOFenabled(1,p,2) && !isDOFfixed(1,p,2)){
          zsqueeze += _startingAlignPlanes[p*6+2] * (p-17.5);
          zsqueeze_constraint += _weakValues[2] * (p/35. - 0.5) * (p-17.5);
        }
        if (isDOFenabled(1,p,5) && !isDOFfixed(1,p,5)){
          ztwist += _startingAlignPlanes[p*6+5] * (p-17.5);
          ztwist_constraint += _weakValues[5] * (p/35. - 0.5) * (p-17.5);
        }
        if (isDOFenabled(1,p,3) && !isDOFfixed(1,p,3)){
          xparallel += _startingAlignPlanes[p*6+3];
          xparallel_constraint += _weakValues[3];
        }
        if (isDOFenabled(1,p,4) && !isDOFfixed(1,p,4)){
          yparallel += _startingAlignPlanes[p*6+4];
          yparallel_constraint += _weakValues[4];
        }
      }else{
        if (isDOFenabled(1,p,0) && !isDOFfixed(1,p,0)){
          xskew -= _startingAlignPlanes[p*6+0] * (p-17.5);
          xskew_constraint += _weakValues[0] * (p/35. - 0.5) * (p-17.5);
        }
        if (isDOFenabled(1,p,1) && !isDOFfixed(1,p,1)){
          yskew += _startingAlignPlanes[p*6+1] * (p-17.5);
          yskew_constraint += _weakValues[1] * (p/35. - 0.5) * (p-17.5);
        }
        if (isDOFenabled(1,p,2) && !isDOFfixed(1,p,2)){
          zsqueeze -= _startingAlignPlanes[p*6+2] * (p-17.5);
          zsqueeze_constraint += _weakValues[2] * (p/35. - 0.5) * (p-17.5);
        }
        if (isDOFenabled(1,p,5) && !isDOFfixed(1,p,5)){
          ztwist -= _startingAlignPlanes[p*6+5] * (p-17.5);
          ztwist_constraint += _weakValues[5] * (p/35. - 0.5) * (p-17.5);
        }
        if (isDOFenabled(1,p,3) && !isDOFfixed(1,p,3)){
          xparallel -= _startingAlignPlanes[p*6+3];
          xparallel_constraint += _weakValues[3];
        }
        if (isDOFenabled(1,p,4) && !isDOFfixed(1,p,4)){
          yparallel += _startingAlignPlanes[p*6+4];
          yparallel_constraint += _weakValues[4];
        }
      }
    }

    // Fix all weak modes to 0
    // or Constrain weak modes by measurements

    // X skew
    bool has_enabled = false;
    for (uint16_t p=0;p<StrawId::_nplanes;++p){
      if (isDOFenabled(1, p, 0) && !isDOFfixed(1, p, 0)) {
        has_enabled = true;
        break;
      }
    }
    if (has_enabled && _weakSigmas[0] >= 0){
      if (_weakConstraints == "Fix"){
        output_file << "Constraint " << xskew_constraint-xskew << std::endl;
      }else{
        output_file << "Measurement " << xskew_constraint-xskew << " " << _weakSigmas[0] << std::endl;
      }
      for (uint16_t p = 0; p < StrawId::_nplanes; ++p) {
        if (!isDOFenabled(1, p, 0) || isDOFfixed(1, p, 0))
          continue;
        // if this plane is rotated, need to invert constraint
        if (nominalTracker.getPlane(p).planeToDS().rotation().getTheta() == 0)
          output_file << getLabel(1, p, 0) << "    " << (p-17.5) << std::endl;
        else
          output_file << getLabel(1, p, 0) << "    " << -1*(p-17.5) << std::endl;
      }
    }
    
    // Y skew
    has_enabled = false;
    for (uint16_t p=0;p<StrawId::_nplanes;++p){
      if (isDOFenabled(1, p, 1) && !isDOFfixed(1, p, 1)) {
        has_enabled = true;
        break;
      }
    }
    if (has_enabled && _weakSigmas[1] >= 0){
      if (_weakConstraints == "Fix"){
        output_file << "Constraint " << yskew_constraint-yskew << std::endl;
      }else{
        output_file << "Measurement " << yskew_constraint-yskew << " " << _weakSigmas[1] << std::endl;
      }
      for (uint16_t p = 0; p < StrawId::_nplanes; ++p) {
        if (!isDOFenabled(1, p, 1) || isDOFfixed(1, p, 1))
          continue;
        output_file << getLabel(1, p, 1) << "    " << (p-17.5) << std::endl;
      }
    }

    // Z squeeze
    has_enabled = false;
    for (uint16_t p=0;p<StrawId::_nplanes;++p){
      if (isDOFenabled(1, p, 2) && !isDOFfixed(1, p, 2)) {
        has_enabled = true;
        break;
      }
    }
    if (has_enabled && _weakSigmas[2] >= 0){
      if (_weakConstraints == "Fix"){
        output_file << "Constraint " << zsqueeze_constraint-zsqueeze << std::endl;
      }else{
        output_file << "Measurement " << zsqueeze_constraint-zsqueeze << " " << _weakSigmas[2] << std::endl;
      }
      for (uint16_t p = 0; p < StrawId::_nplanes; ++p) {
        if (!isDOFenabled(1, p, 2) || isDOFfixed(1, p, 2))
          continue;
        // if this plane is rotated, need to invert constraint
        if (nominalTracker.getPlane(p).planeToDS().rotation().getTheta() == 0)
          output_file << getLabel(1, p, 2) << "    " << (p-17.5) << std::endl;
        else
          output_file << getLabel(1, p, 2) << "    " << -1*(p-17.5) << std::endl;
      }
    }
    // X parallel
    has_enabled = false;
    for (uint16_t p=0;p<StrawId::_nplanes;++p){
      if (isDOFenabled(1, p, 3) && !isDOFfixed(1, p, 3)) {
        has_enabled = true;
        break;
      }
    }
    if (has_enabled && _weakSigmas[3] >= 0){
      if (_weakConstraints == "Fix"){
        output_file << "Constraint " << xparallel_constraint-xparallel << std::endl;
      }else{
        output_file << "Measurement " << xparallel_constraint-xparallel << " " << _weakSigmas[3] << std::endl;
      }
      for (uint16_t p = 0; p < StrawId::_nplanes; ++p) {
        if (!isDOFenabled(1, p, 3) || isDOFfixed(1, p, 3))
          continue;
        // if this plane is rotated, need to invert constraint
        if (nominalTracker.getPlane(p).planeToDS().rotation().getTheta() == 0)
          output_file << getLabel(1, p, 3) << "    " << 1 << std::endl;
        else
          output_file << getLabel(1, p, 3) << "    " << -1 << std::endl;
      }
    }
    // Y parallel
    has_enabled = false;
    for (uint16_t p=0;p<StrawId::_nplanes;++p){
      if (isDOFenabled(1, p, 4) && !isDOFfixed(1, p, 4)) {
        has_enabled = true;
        break;
      }
    }
    if (has_enabled && _weakSigmas[4] >= 0){
      if (_weakConstraints == "Fix"){
        output_file << "Constraint " << yparallel_constraint-yparallel << std::endl;
      }else{
        output_file << "Measurement " << yparallel_constraint-yparallel << " " << _weakSigmas[4] << std::endl;
      }
      for (uint16_t p = 0; p < StrawId::_nplanes; ++p) {
        if (!isDOFenabled(1, p, 4) || isDOFfixed(1, p, 4))
          continue;
        output_file << getLabel(1, p, 4) << "    " << 1 << std::endl;
      }
    }
    // z twist
    has_enabled = false;
    for (uint16_t p=0;p<StrawId::_nplanes;++p){
      if (isDOFenabled(1, p, 5) && !isDOFfixed(1, p, 5)) {
        has_enabled = true;
        break;
      }
    }
    if (has_enabled && _weakSigmas[5] >= 0){
      if (_weakConstraints == "Fix"){
        output_file << "Constraint " << ztwist_constraint-ztwist << std::endl;
      }else{
        output_file << "Measurement " << ztwist_constraint-ztwist << " " << _weakSigmas[5] << std::endl;
      }
      for (uint16_t p = 0; p < StrawId::_nplanes; ++p) {
        if (!isDOFenabled(1, p, 5) || isDOFfixed(1, p, 5))
          continue;
        // if this plane is rotated, need to invert constraint
        if (nominalTracker.getPlane(p).planeToDS().rotation().getTheta() == 0)
          output_file << getLabel(1, p, 5) << "    " << (p-17.5) << std::endl;
        else
          output_file << getLabel(1, p, 5) << "    " << -1*(p-17.5) << std::endl;
      }
    }
  }

}

void AlignKKTrackCollector::writeMillepedeSteering() {
  std::ofstream output_file(_steerFilename);

  output_file << "! Steering file generated by AlignKKTrackCollector" << std::endl
              << "Cfiles" << std::endl
              << _paramFilename << std::endl
              << _constrFilename << std::endl
              << _milleFilename << std::endl
              << std::endl;

  for (std::string const& line : _steerLines) {
    output_file << line << std::endl;
  }
  
  //FIXME
  output_file << std::endl
              << "method inversion 10 0.001" << std::endl
              << "end" << std::endl;
}

void AlignKKTrackCollector::writeMillepedeParams() {
  // write a params.txt telling millepede what the alignment constants
  // were for this track collection iteration
  std::ofstream output_file(_paramFilename);
  output_file << "! Generated by AlignKKTrackCollector" << std::endl;
  output_file << "! columns: label, alignment parameter start value, "
                 "presigma (-ve: fixed, 0: variable)" << std::endl;

  output_file << "Parameter" << std::endl;
  output_file << "! Plane DOFs:" << std::endl;

  for (uint16_t p = 0; p < StrawId::_nplanes; ++p) {
    for (size_t dof_n = 0; dof_n < _dof_per_plane; dof_n++) {
      if (!isDOFenabled(1, p, dof_n)) {
        continue;
      }

      // label(int)   initial value (float)   presigma (-ve: fixed, 0: variable)
      output_file << getLabel(1, p, dof_n) << "  " << 0 << "  " 
                  << (isDOFfixed(1, p, dof_n) ? "-1" : "0") << std::endl;
    }
  }
  output_file << std::endl << "! end plane DOFs" << std::endl;

  output_file << std::endl << "! Panel DOFs:" << std::endl;

  for (uint16_t p = 0; p < StrawId::_nupanels; ++p) {
    for (size_t dof_n = 0; dof_n < _dof_per_panel; dof_n++) {
      if (!isDOFenabled(2, p, dof_n)) {
        continue;
      }

      // label(int)   initial value (float)   presigma (-ve: fixed, 0: variable)
      output_file << getLabel(2, p, dof_n) << "  " << 0 << "  " << (isDOFfixed(2, p, dof_n) ? "-1" : "0") << std::endl;
    }
  }
  output_file << std::endl << "! end panel DOFs" << std::endl;

  output_file.close();
}

void AlignKKTrackCollector::writeMillepedeExtras() {
  std::ofstream output_file2(_extraFilename);
  for (uint16_t p = 0; p < StrawId::_nplanes; ++p) {
    for (size_t dof_n = 0; dof_n < _dof_per_plane; dof_n++) {
      if (!isDOFenabled(1, p, dof_n)) {
        continue;
      }

      // label(int)   initial value (float)   presigma (-ve: fixed, 0: variable)
      output_file2 << getLabel(1, p, dof_n) << "  " << _startingAlignPlanes[p*6+dof_n] << "  " << _DOFCounts[getLabel(1,p,dof_n)] << std::endl;
    }
  }
  for (uint16_t p = 0; p < StrawId::_nupanels; ++p) {
    for (size_t dof_n = 0; dof_n < _dof_per_panel; dof_n++) {
      if (!isDOFenabled(2, p, dof_n)) {
        continue;
      }

      // label(int)   initial value (float)   presigma (-ve: fixed, 0: variable)
      output_file2 << getLabel(2, p, dof_n) << "  " << _startingAlignPanels[p*6+dof_n] << "  " << _DOFCounts[getLabel(2,p,dof_n)] << std::endl;
    }
  }

  output_file2.close();
}

void AlignKKTrackCollector::cacheStartingParams(TrkAlignPlane const& alignConstPlanes,
                                               TrkAlignPanel const& alignConstPanels,
                                               TrackerStatus const& trackerStatus) {
  _startingAlignPlanes = std::vector<float>(StrawId::_nplanes*6,0);
  _startingAlignPanels = std::vector<float>(StrawId::_nupanels*6,0);
  for (uint16_t p = 0; p < StrawId::_nplanes; ++p) {
    auto const& rowpl = alignConstPlanes.rowAt(p);

    std::vector<float> pl_consts{rowpl.dx(), rowpl.dy(), rowpl.dz(),
                                 rowpl.rx(), rowpl.ry(), rowpl.rz()};

    for (size_t i=0;i<6;i++) 
      _startingAlignPlanes[p*6+i] = pl_consts[i];
  }
  for (uint16_t p = 0; p < StrawId::_nupanels; ++p) {
    auto const& rowpa = alignConstPanels.rowAt(p);

    std::vector<float> pa_consts{rowpa.dx(), rowpa.dy(), rowpa.dz(),
                                 rowpa.rx(), rowpa.ry(), rowpa.rz()};

    for (size_t i=0;i<6;i++) 
      _startingAlignPanels[p*6+i] = pa_consts[i];
  }

  _panelActive = std::vector<bool>(StrawId::_nupanels,true);
  static StrawStatus mask = StrawStatus(StrawStatus::absent) |
    StrawStatus(StrawStatus::nowire) |
    StrawStatus(StrawStatus::noHV) |
    StrawStatus(StrawStatus::noLV) |
    StrawStatus(StrawStatus::nogas) |
    StrawStatus(StrawStatus::lowgasgain) |
    StrawStatus(StrawStatus::noHVPreamp) |
    StrawStatus(StrawStatus::noCalPreamp) |
    StrawStatus(StrawStatus::disabled);
  for (size_t pl=0;pl<StrawId::_nplanes;pl++){
    for (size_t pa=0;pa<StrawId::_npanels;pa++){
      if (trackerStatus.panelStatus(StrawId(pl,pa,0)).hasAnyProperty(mask)){
        _panelActive[pl*StrawId::_npanels+pa] = false;
      }
    }
  }
}

} // namespace mu2e

using mu2e::AlignKKTrackCollector;
DEFINE_ART_MODULE(AlignKKTrackCollector)
