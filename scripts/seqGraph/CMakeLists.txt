cmake_minimum_required(VERSION 3.10)
project(seqGraph VERSION 1.0.1 DESCRIPTION "seq graph")

set(CMAKE_CXX_STANDARD 14)
set(DEFAULT_BUILD_TYPE RELEASE)
if(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
    message(STATUS "Setting build type to '${DEFAULT_BUILD_TYPE}' as none was specified.")
    set(CMAKE_BUILD_TYPE "${DEFAULT_BUILD_TYPE}" CACHE STRING "Choose the type of build." FORCE)
    # Set the possible values of build type for cmake-gui
    set_property(CACHE CMAKE_BUILD_TYPE PROPERTY STRINGS "Debug" "Release" "MinSizeRel" "RelWithDebInfo")
endif()
include(GNUInstallDirs)
set(CMAKE_CXX_FLAGS_DEBUG "$ENV{CXX_FLAGS} -O0 -Wall -ggdb -fkeep-inline-functions")
set(CMAKE_CXX_FLAGS_RELEASE "$ENV{CXXFLAGS} -O3 ")

set(SOURCE_FILES src/Edge.cpp
        src/EndPoint.cpp
        src/Exception.cpp
        src/Graph.cpp
        src/Junction.cpp
        src/Vertex.cpp
        src/Weight.cpp
        src/matching.cpp include/cxxopts.hpp)

add_library(${PROJECT_NAME} SHARED ${SOURCE_FILES} include/cxxopts.hpp)

add_executable(matching main.cpp ${SOURCE_FILES} include/util.h)
add_executable(rDistance utils/reads_overlap.cpp utils/reads_overlap.h include/util.h)
add_executable(eref extract_ref.cpp )
add_executable(eResult utils/e_result_reads.cpp utils/e_result_reads.h)
add_executable(bDistance utils/bDistance.cpp utils/bDistance.h include/util.h)
target_link_libraries(rDistance -lhts)
target_link_libraries(eResult -lhts)
target_link_libraries(bDistance -lhts)
target_link_libraries(matching)
target_link_libraries(eref -lpthread)
target_include_directories(${PROJECT_NAME} PUBLIC
        $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
        $<INSTALL_INTERFACE:include>)

set_target_properties(${PROJECT_NAME} PROPERTIES
        VERSION ${PROJECT_VERSION}
        SOVERSION 1)

install(TARGETS ${PROJECT_NAME} EXPORT seqGraphConfig
        ARCHIVE  DESTINATION ${CMAKE_INSTALL_LIBDIR}
        LIBRARY  DESTINATION ${CMAKE_INSTALL_LIBDIR}
        RUNTIME  DESTINATION ${CMAKE_INSTALL_BINDIR})
install(DIRECTORY include/ DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/${PROJECT_NAME})

install(EXPORT seqGraphConfig DESTINATION share/seqGraph/cmake)

export(TARGETS ${PROJECT_NAME} FILE seqGraphConfig.cmake)

