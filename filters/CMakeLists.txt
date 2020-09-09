art_make( MODULE_LIBRARIES
            artg4_util
            artg4_geantInit
            gm2tracker_utils
	          gm2util_coordSystems
            artg4_gm2Geane
	    gm2geom_fields
	    gm2geom_fields_gm2FieldManager_service
        artg4_pluginActions_physicalVolumeStore_physicalVolumeStore_service
          SERVICE_LIBRARIES
            gm2geom_coordSystems
            gm2geom_coordSystems_CoordSystemsStore_service
	          gm2util_common
	          gm2geom_fields
        )

install_headers()

#copied from reco initially - should remove unessecary libraries from here
