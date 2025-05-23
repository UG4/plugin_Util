# Copyright (c) 2011-2015:  G-CSC, Goethe University Frankfurt
# Author: Sebastian Reiter
# 
# This file is part of UG4.
# 
# UG4 is free software: you can redistribute it and/or modify it under the
# terms of the GNU Lesser General Public License version 3 (as published by the
# Free Software Foundation) with the following additional attribution
# requirements (according to LGPL/GPL v3 §7):
# 
# (1) The following notice must be displayed in the Appropriate Legal Notices
# of covered and combined works: "Based on UG4 (www.ug4.org/license)".
# 
# (2) The following notice must be displayed at a prominent place in the
# terminal output of covered works: "Based on UG4 (www.ug4.org/license)".
# 
# (3) The following bibliography is recommended for citation and must be
# preserved in all covered files:
# "Reiter, S., Vogel, A., Heppner, I., Rupp, M., and Wittum, G. A massively
#   parallel geometric multigrid solver on hierarchically distributed grids.
#   Computing and visualization in science 16, 4 (2013), 151-164"
# "Vogel, A., Reiter, S., Rupp, M., Nägel, A., and Wittum, G. UG4 -- a novel
#   flexible software system for simulating pde based models on high performance
#   computers. Computing and visualization in science 16, 4 (2013), 165-179"
# 
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
# GNU Lesser General Public License for more details.

################################################################################
# SAMPLE
#
# This script shows how easy it is to add a build script for your plugins.
# Simply specify the plugins name and a list of your source files.
################################################################################

# Name of your plugin and sources.
set(pluginName	Util)
set(SOURCES		domain_util.cpp util_plugin.cpp)

set (CMAKE_CXX_STANDARD 14)

# Note: 
# In C++ 17, std::unary_function has been removed. Unfortunately, UG4's boost is outdated... 
# FYI Starting with 1.76.0 Configure boost:
# https://www.boost.org/doc/libs/1_76_0/libs/config/doc/html/boost_config/boost_macro_reference.html 
# set (CMAKE_CXX_STANDARD 17)
# add_compile_definitions(BOOST_NO_CXX98_FUNCTION_BASE)

################################################################################
# The code below doesn't have to be changed (usually)
################################################################################
# cmake_minimum_required(VERSION 2.6)

project(UG_PLUGIN_${pluginName})

# include the definitions and dependencies for ug-plugins.
include(${UG_ROOT_CMAKE_PATH}/ug_plugin_includes.cmake)

################################################################################
# Classic binding (static or dynamic plugin).
################################################################################
if(NOT USE_PYBIND11)
if(buildEmbeddedPlugins)
	# add the sources to ug4's sources
	EXPORTSOURCES(${CMAKE_CURRENT_SOURCE_DIR} ${SOURCES})
else(buildEmbeddedPlugins)
	# create a shared library from the sources and link it against ug4.
	add_library(${pluginName} SHARED ${SOURCES})
	target_link_libraries (${pluginName} ug4)
endif(buildEmbeddedPlugins)
endif(NOT USE_PYBIND11)

################################################################################
# Python binding (static plugin, dynamic python interface).
################################################################################
if(USE_PYBIND11)

	# Static plugin.
	# include_directories(${SUPERLU_INCLUDE_DIRS})
	add_library(${pluginName} STATIC ${SOURCES})

	SET(myPluginSources util_plugin_pybind.cpp)
	# SET(myLibraries ug4 SuperLU6)
	
	# First argument must match module name in PYBIND11_MODULE call
	python_add_library(pyutil MODULE ${myPluginSources} ${SOURCES} WITH_SOABI)
	# ug4pybind_add_module(pyutil ${myPluginSources} ${myLibraries})
	target_link_libraries (pyutil PRIVATE ${pluginName} ug4_s)
	set_target_properties(pyutil PROPERTIES INSTALL_RPATH "$ORIGIN/..:$ORIGIN/../../../lib")
	install(TARGETS pyutil LIBRARY DESTINATION ug4py COMPONENT pymodules)
endif(USE_PYBIND11)
