/*
 * Copyright (c) 2011-2015:  G-CSC, Goethe University Frankfurt
 * Author: Sebastian Reiter
 * 
 * This file is part of UG4.
 * 
 * UG4 is free software: you can redistribute it and/or modify it under the
 * terms of the GNU Lesser General Public License version 3 (as published by the
 * Free Software Foundation) with the following additional attribution
 * requirements (according to LGPL/GPL v3 §7):
 * 
 * (1) The following notice must be displayed in the Appropriate Legal Notices
 * of covered and combined works: "Based on UG4 (www.ug4.org/license)".
 * 
 * (2) The following notice must be displayed at a prominent place in the
 * terminal output of covered works: "Based on UG4 (www.ug4.org/license)".
 * 
 * (3) The following bibliography is recommended for citation and must be
 * preserved in all covered files:
 * "Reiter, S., Vogel, A., Heppner, I., Rupp, M., and Wittum, G. A massively
 *   parallel geometric multigrid solver on hierarchically distributed grids.
 *   Computing and visualization in science 16, 4 (2013), 151-164"
 * "Vogel, A., Reiter, S., Rupp, M., Nägel, A., and Wittum, G. UG4 -- a novel
 *   flexible software system for simulating pde based models on high performance
 *   computers. Computing and visualization in science 16, 4 (2013), 165-179"
 * 
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
 * GNU Lesser General Public License for more details.
 */


// Own headers.
#include "domain_util.h"

// replace this with util_domain_dependent.h or util_algebra_dependent.h
// to speed up compilation time
#include "bridge/bridge.h"
#include "bridge/util.h"
#include "bridge/util_domain_dependent.h"
#include "bridge/util_algebra_dependent.h"
// #include "bridge/util_domain_algebra_dependent.h"

#include "common/ug_config.h"
#include "common/error.h"


#include "lib_disc/domain.h"
#include "lib_disc/domain_util.h"

#include "ug.h"

#include <string>

using namespace std;
using namespace ug::bridge;

#ifdef UG_PARALLEL
#include "pcl/pcl_util.h"
#endif

namespace ug{
namespace Util{

/** 
 *  \defgroup util_plugin Util
 *  \ingroup plugins
 *  This is a small sample plugin.
 *  \{
 */


/**
 * Class exporting the functionality of the plugin. All functionality that is to
 * be used in scripts or visualization must be registered here.
 */
struct Functionality
{

/**
 * Function called for the registration of Domain and Algebra dependent parts
 * of the plugin. All Functions and Classes depending on both Domain and Algebra
 * are to be placed here when registering. The method is called for all
 * available Domain and Algebra types, based on the current build options.
 *
 * @param reg				registry
 * @param parentGroup		group for sorting of functionality
 */
template <typename TDomain, typename TAlgebra, typename TRegistry=ug::bridge::Registry>
static void DomainAlgebra(TRegistry& reg, string grp)
{
//	useful defines
	string suffix = GetDomainAlgebraSuffix<TDomain,TAlgebra>();
	string tag = GetDomainAlgebraTag<TDomain,TAlgebra>();

}

/**
 * Function called for the registration of Domain dependent parts
 * of the plugin. All Functions and Classes depending on the Domain
 * are to be placed here when registering. The method is called for all
 * available Domain types, based on the current build options.
 *
 * @param reg				registry
 * @param parentGroup		group for sorting of functionality
 */
template <typename TDomain, typename TRegistry=ug::bridge::Registry>
static void Domain(TRegistry& reg, string grp)
{
//	useful defines
	string suffix = GetDomainSuffix<TDomain>();
	string tag = GetDomainTag<TDomain>();

	reg.add_function("CheckSubsets", &CheckSubsets<TDomain>);

	{
		// This is more complicated, as the return type should be included in name.
		std::string grpname = std::string("CreateDomain");
		reg.add_function(grpname.append(suffix),
				static_cast<SmartPtr<TDomain>(*)(const std::string&, int, const std::vector<std::string> &)>
			(&CreateDomain<TDomain>), grpname);
	}
}

/**
 * Function called for the registration of Dimension dependent parts
 * of the plugin. All Functions and Classes depending on the Dimension
 * are to be placed here when registering. The method is called for all
 * available Dimension types, based on the current build options.
 *
 * @param reg				registry
 * @param parentGroup		group for sorting of functionality
 */
template <int dim, typename TRegistry=ug::bridge::Registry>
static void Dimension(TRegistry& reg, string grp)
{
//	useful defines
	string suffix = GetDimensionSuffix<dim>();
	string tag = GetDimensionTag<dim>();

}

/**
 * Function called for the registration of Algebra dependent parts
 * of the plugin. All Functions and Classes depending on Algebra
 * are to be placed here when registering. The method is called for all
 * available Algebra types, based on the current build options.
 *
 * @param reg				registry
 * @param parentGroup		group for sorting of functionality
 */
template <typename TAlgebra, typename TRegistry=ug::bridge::Registry>
static void Algebra(TRegistry& reg, string grp)
{
//	useful defines
	string suffix = GetAlgebraSuffix<TAlgebra>();
	string tag = GetAlgebraTag<TAlgebra>();

}

/**
 * Function called for the registration of Domain and Algebra independent parts
 * of the plugin. All Functions and Classes not depending on Domain and Algebra
 * are to be placed here when registering.
 *
 * @param reg				registry
 * @param parentGroup		group for sorting of functionality
 */
template <typename TRegistry=ug::bridge::Registry>
static void Common(TRegistry& reg, string grp)
{
	//reg.add_function("PluginSaysHello", &PluginSaysHello, grp)
	//	.add_function("PluginCrashes", &PluginCrashes, grp)
	//	.add_function("PluginCrashesFatal", &PluginCrashesFatal, grp);
}

}; // end Functionality

// end group sample_plugin
/// \}

}// end of namespace Util


template <typename TRegistry=ug::bridge::Registry>
void RegisterBridge_Util(TRegistry& reg, string grp)
{
	grp.append("/Util");
	typedef Util::Functionality Functionality;

	try{
#ifndef UG_USE_PYBIND11
		//RegisterCommon<Functionality>(*reg,grp);
		//RegisterDimensionDependent<Functionality>(*reg,grp);
		RegisterDomainDependent<Functionality>(reg,grp);
		//RegisterAlgebraDependent<Functionality>(*reg,grp);
		//RegisterDomainAlgebraDependent<Functionality>(*reg,grp);
#else
		//RegisterCommon<Functionality, TRegistry>(*reg,grp);
		//RegisterDimensionDependent<Functionality, TRegistry>(*reg,grp);
		RegisterDomainDependent<Functionality, TRegistry>(reg,grp);
		//RegisterAlgebraDependent<Functionality>(*reg,grp);
		//RegisterDomainAlgebraDependent<Functionality>(*reg,grp);
#endif

	}
	UG_REGISTRY_CATCH_THROW(grp);
}



#ifdef UG_USE_PYBIND11 // Expose for pybind11.
namespace Util{
	void InitUGPlugin(ug::pybind::Registry* reg, string grp)
	{
		ug::RegisterBridge_Util<ug::pybind::Registry>(*reg, grp);
	}
}
#endif


}//	end of namespace ug


/**
 * This function is called when the plugin is loaded.
 */
extern "C" void InitUGPlugin_Util(Registry* reg, string grp)
{ ug::RegisterBridge_Util(*reg, grp); }

extern "C" UG_API void FinalizeUGPlugin_Util()
{}

