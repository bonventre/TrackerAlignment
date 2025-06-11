cet_make_library(
    SOURCE
      SConscript
      AlignmentUtilities.cc
    LIBRARIES PUBLIC
      Offline::CosmicReco
      Offline::MCDataProducts
      Offline::RecoDataProducts
      Offline::DataProducts
      Offline::GeometryService
      Offline::ProditionsService
      Offline::TrackerConditions
      Offline::TrackerGeom
      Offline::Mu2eUtilities
      Offline::GeneralUtilities
      Offline::RecoDataProducts

      LIBRARY_NAME AlignmentUtilities
)



cet_build_plugin(AlignTrackCollector art::module
  REG_SOURCE AlignTrackCollector_module.cc
    LIBRARIES REG
      tracker-alignment::AlignmentUtilities
      Offline::CosmicReco
      Offline::MCDataProducts
      Offline::RecoDataProducts
      Offline::DataProducts
      Offline::GeometryService
      Offline::ProditionsService
      Offline::TrackerConditions
      Offline::TrackerGeom
      Offline::Mu2eUtilities
      Offline::GeneralUtilities
      Offline::RecoDataProducts
      art_root_io::TFileService_service
      ROOT::Postscript
)

install_source(SUBDIRS src)
