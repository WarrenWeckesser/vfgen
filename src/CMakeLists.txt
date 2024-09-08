#
#  CMakeLists.txt
#
#  This file is to be used with CMake to build VFGEN.
#
#  Copyright (C) 2008 Warren Weckesser
#
#  This program is free software; you can redistribute it and/or modify
#  it under the terms of the GNU General Public License, Version 2, as
#  published by the Free Software Foundation.
#
#  This program is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#
#  You should have received a copy of the GNU General Public License along
#  with this program.  The file LICENSE-GPL2 in the VFGEN source directory
#  contains a copy. You may also obtain a copy by writing to the Free Software
#  Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.
#

cmake_minimum_required(VERSION 2.8.12...3.18.4)

project(vfgen CXX)
#
# Reminders:
#
# :: Define the environment variables CMAKE_INCLUDE_PATH and CMAKE_LIBRARY_PATH if
#    ginac or mxml are installed in locations other than /usr/lib or /usr/local/lib.
# :: Define the environment variable CMAKE_INSTALL_PREFIX if you want 'make install'
#    to install the vfgen executable in a directory other than /usr/local/bin.
#    The default prefix is /usr/local; vfgen will be installed in
#        ${CMAKE_INSTALL_PREFIX}/bin
#

#
# Set the version number here, in the form "MAJOR.MINOR.PATCH".
# Then set up a definition of the macro VERSION that will be available
# in the source code.
#
#SET(VFGEN_VERSION_MAJOR "2")
#SET(VFGEN_VERSION_MINOR "5")
#SET(VFGEN_VERSION_PATCH "0-dev")
#SET(VFGEN_VERSION ${VFGEN_VERSION_MAJOR}.${VFGEN_VERSION_MINOR}.${VFGEN_VERSION_PATCH})

file(STRINGS "../version.txt" VFGEN_VERSION)
add_definitions(-DVERSION="${VFGEN_VERSION}")

set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++17")

if(CMAKE_COMPILER_IS_GNUCXX)
    add_definitions(-Wall)
    add_definitions(-Wextra)
    add_definitions(-Wpedantic)
    add_definitions(-Werror)
endif(CMAKE_COMPILER_IS_GNUCXX)

#
# For convenience, define the list of C++ files in a macro.
#
set(VFGEN_CPP_SRC_FILES codegen_utils.cpp ginac_aux_functions.cpp strutils.cpp MyVec.cpp vf_adolc.cpp vf_auto.cpp vf_boostodeint.cpp vf.cpp vf_cvode.cpp vf_cvode7.cpp vf_dde23.cpp vf_ddebiftool.cpp vf_dde_solver.cpp vf_delay2ode.cpp vf_differentials.cpp vf_dstool.cpp vf_evf.cpp vfgen.cpp vf_gsl.cpp vf_javamath.cpp vf_javascript.cpp vf_julia.cpp vf_latex.cpp vf_lsoda.cpp vf_matcont.cpp vf_matlab.cpp vf_octave.cpp vf_pddecont.cpp vf_pydstool.cpp vf_pygsl.cpp vf_r.cpp vf_radau5.cpp vf_scilab.cpp vf_scipy.cpp vf_taylor.cpp vf_xml.cpp vf_xpp.cpp)

add_executable(vfgen ${VFGEN_CPP_SRC_FILES})

#
# Look for the required libraries.
# VFGEN requires Mini-XML (www.easysw.com/~mike/mxml) and GiNaC (www.ginac.de).
#
find_path(GINAC_INCLUDE_DIR ginac/ginac.h /usr/local/include /usr/include)
find_path(MXML_INCLUDE_DIR mxml.h /usr/local/include /usr/include)

find_library(GINAC_LIBRARY ginac /usr/local/lib /usr/lib)
find_library(MXML_LIBRARY mxml /usr/local/lib /usr/lib)

if(GINAC_INCLUDE_DIR AND GINAC_LIBRARY)
    message("Found ginac.")
else(GINAC_INCLUDE_DIR AND GINAC_LIBRARY)
    message(SEND_ERROR "missing ginac library")
endif(GINAC_INCLUDE_DIR AND GINAC_LIBRARY)

if(MXML_INCLUDE_DIR AND MXML_LIBRARY)
    message("Found mxml.")
else(MXML_INCLUDE_DIR AND MXML_LIBRARY)
    message(SEND_ERROR "missing mxml library")
endif(MXML_INCLUDE_DIR AND MXML_LIBRARY)

include_directories(${GINAC_INCLUDE_DIR} ${MXML_INCLUDE_DIR})

target_link_libraries(vfgen ${GINAC_LIBRARY} ${MXML_LIBRARY})

set(CMAKE_INSTALL_PREFIX $ENV{CMAKE_INSTALL_PREFIX})
install(TARGETS vfgen RUNTIME DESTINATION bin)


#SET(CPACK_PACKAGE_DESCRIPTION_SUMMARY "VFGEN is a vector field source code generator.")
#SET(CPACK_PACKAGE_VENDOR "Warren Weckesser")
#SET(CPACK_PACKAGE_VERSION_MAJOR ${VFGEN_VERSION_MAJOR})
#SET(CPACK_PACKAGE_VERSION_MINOR ${VFGEN_VERSION_MINOR})
#SET(CPACK_PACKAGE_VERSION_PATCH ${VFGEN_VERSION_PATCH})
#SET(CPACK_SOURCE_GENERATOR "TGZ")
#INCLUDE(CPack)
