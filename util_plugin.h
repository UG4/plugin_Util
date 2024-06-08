# ifndef __util__util_bridge_h__
# define __util__util_bridge_h__

#pragma once

#include "bridge/util.h"

extern "C" void InitUGPlugin_Util(ug::bridge::Registry* reg, std::string grp);
extern "C" void FinalizeUGPlugin_Util();

#ifdef UG_USE_PYBIND11

#include "bindings/pybind/ug_pybind.h"

namespace ug {
namespace bridge {
namespace UtilBridge{
	void InitUGPlugin(ug::pybind::Registry* reg, std::string grp);
}
}
}
#endif


# endif //__util__util_bridge_h__
