SET (LIB_${PROJECT_NAME}_PYSOURCES 
        I3MuonGun.cxx
        CanCan.cxx
        Generator.cxx
        Surface.cxx
        Flux.cxx
        RadialDistribution.cxx
        EnergyDistribution.cxx
        WeightCalculator.cxx
        CORSIKAGenerationProbability.cxx
        module.cxx
)
IF (IS_DIRECTORY ${CMAKE_SOURCE_DIR}/PROPOSAL)
    LIST (APPEND LIB_${PROJECT_NAME}_PYSOURCES
        CompactTrack.cxx
        TrackBinner.cxx
        MuonPropagator.cxx
        histogram.cxx
    )
    ADD_DEFINITIONS(-DUSE_PROPOSAL)
ENDIF (IS_DIRECTORY ${CMAKE_SOURCE_DIR}/PROPOSAL)

i3_add_pybindings(MuonGun
  ${LIB_${PROJECT_NAME}_PYSOURCES}
  USE_TOOLS boost python log4cplus
  USE_PROJECTS MuonGun
  )

