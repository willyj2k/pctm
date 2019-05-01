# Install script for directory: /Users/pkmnfreak/pctm/CGL/src

# Set the install prefix
if(NOT DEFINED CMAKE_INSTALL_PREFIX)
  set(CMAKE_INSTALL_PREFIX "/usr/local")
endif()
string(REGEX REPLACE "/$" "" CMAKE_INSTALL_PREFIX "${CMAKE_INSTALL_PREFIX}")

# Set the install configuration name.
if(NOT DEFINED CMAKE_INSTALL_CONFIG_NAME)
  if(BUILD_TYPE)
    string(REGEX REPLACE "^[^A-Za-z0-9_]+" ""
           CMAKE_INSTALL_CONFIG_NAME "${BUILD_TYPE}")
  else()
    set(CMAKE_INSTALL_CONFIG_NAME "")
  endif()
  message(STATUS "Install configuration: \"${CMAKE_INSTALL_CONFIG_NAME}\"")
endif()

# Set the component getting installed.
if(NOT CMAKE_INSTALL_COMPONENT)
  if(COMPONENT)
    message(STATUS "Install component: \"${COMPONENT}\"")
    set(CMAKE_INSTALL_COMPONENT "${COMPONENT}")
  else()
    set(CMAKE_INSTALL_COMPONENT)
  endif()
endif()

# Is this installation the result of a crosscompile?
if(NOT DEFINED CMAKE_CROSSCOMPILING)
  set(CMAKE_CROSSCOMPILING "FALSE")
endif()

if("x${CMAKE_INSTALL_COMPONENT}x" STREQUAL "xUnspecifiedx" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/lib" TYPE STATIC_LIBRARY FILES "/Users/pkmnfreak/pctm/build/CGL/src/libCGL.a")
  if(EXISTS "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib/libCGL.a" AND
     NOT IS_SYMLINK "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib/libCGL.a")
    execute_process(COMMAND "/Library/Developer/CommandLineTools/usr/bin/ranlib" "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib/libCGL.a")
  endif()
endif()

if("x${CMAKE_INSTALL_COMPONENT}x" STREQUAL "xUnspecifiedx" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/include/CGL" TYPE FILE FILES
    "/Users/pkmnfreak/pctm/CGL/src/CGL.h"
    "/Users/pkmnfreak/pctm/CGL/src/vector2D.h"
    "/Users/pkmnfreak/pctm/CGL/src/vector3D.h"
    "/Users/pkmnfreak/pctm/CGL/src/vector4D.h"
    "/Users/pkmnfreak/pctm/CGL/src/matrix3x3.h"
    "/Users/pkmnfreak/pctm/CGL/src/matrix4x4.h"
    "/Users/pkmnfreak/pctm/CGL/src/quaternion.h"
    "/Users/pkmnfreak/pctm/CGL/src/complex.h"
    "/Users/pkmnfreak/pctm/CGL/src/color.h"
    "/Users/pkmnfreak/pctm/CGL/src/osdtext.h"
    "/Users/pkmnfreak/pctm/CGL/src/viewer.h"
    "/Users/pkmnfreak/pctm/CGL/src/base64.h"
    "/Users/pkmnfreak/pctm/CGL/src/tinyxml2.h"
    "/Users/pkmnfreak/pctm/CGL/src/renderer.h"
    )
endif()

