########################################################
# cmake file for building k4TrkAlgorithms Package
CMAKE_MINIMUM_REQUIRED(VERSION 3.11 FATAL_ERROR)
########################################################

# project name
PROJECT( k4TrkAlgorithms )


# project version
SET( ${PROJECT_NAME}_VERSION_MAJOR 0 )
SET( ${PROJECT_NAME}_VERSION_MINOR 3 )
SET( ${PROJECT_NAME}_VERSION_PATCH 0 )

### LIBRARY #################################################################

FILE ( GLOB _plugin_sources src/*.c* )

gaudi_add_module(k4TrkAlgorithmsPlugins
                 SOURCES	${_plugin_sources}
                 LINK   	Gaudi::GaudiKernel
                        	EDM4HEP::edm4hep
                        	ROOT::Hist
                        	${ROOT_LIBRARIES}
                        	DD4hep::DDCore
                        	k4FWCore::k4FWCore)

INSTALL( TARGETS k4TrkAlgorithmsPlugins
         EXPORT ${CMAKE_PROJECT_NAME}TARGETS
         RUNTIME DESTINATION "${CMAKE_INSTALL_BINDIR}" COMPONENT bin
         LIBRARY DESTINATION "${CMAKE_INSTALL_LIBDIR}" COMPONENT shlib
         PUBLIC_HEADER DESTINATION "${CMAKE_INSTALL_INCLUDEDIR}/@{CMAKE_PROJECT_NAME}"
         COMPONENT dev )


# include directories
TARGET_INCLUDE_DIRECTORIES(k4TrkAlgorithmsPlugins 
                   PRIVATE ./include 
                           ${CMAKE_CURRENT_BINARY_DIR}
                           ${Gaudi_INCLUDE_DIRS}
                           ${ROOT_INCLUDE_DIRS})

SET( GAUDI_GENCONF_DIR "genConfDir" )
if (${Gaudi_VERSION} VERSION_LESS 35.1)
	SET( GAUDI_GENCONF_DIR "genConf" )
endif()

