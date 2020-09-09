art_make( SERVICE_LIBRARIES 
          gm2geom_common_Gm2Constants_service 
	  gm2geom_fields
	  gm2geom_fields_gm2FieldManager_service
          LIB_LIBRARIES
          ${GSL}
          ${GSLCBLAS} 
	  gm2geom_fields
	  gm2geom_fields_gm2FieldManager_service
        )
        
install_headers() 
