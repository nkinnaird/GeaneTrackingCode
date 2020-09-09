# build libraries
art_make( MODULE_LIBRARIES
          artg4_util
          artg4_geantInit
          gm2geom_coordSystems
          gm2geom_coordSystems_CoordSystemsStore_service
	  artg4_pluginActions_physicalVolumeStore_physicalVolumeStore_service
	  SERVICE_LIBRARIES
          gm2geom_common_Gm2Constants_service 
          LIB_LIBRARIES
          ${GSL}
          ${GSLCBLAS}
          artg4_util
          artg4_geantInit
          artg4_gm2Geane
          gm2tracker_fitter
          gm2geom_strawtracker
          gm2geom_coordSystems
          gm2geom_coordSystems_CoordSystemsStore_service
	  artg4_pluginActions_physicalVolumeStore_physicalVolumeStore_service
        )

# copy headers
install_headers() 
