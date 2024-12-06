

#include "util_plugin.h"


#ifdef UG_USE_PYBIND11

PYBIND11_MODULE(pyutil, m)
{
	m.doc() = "Util module";
	m.attr("__name__") = "ug4py.pyutil";

	ug::pybind::Registry registry(m);
	std::string name("Utilities");

	ug::Util::InitUGPlugin(&registry, name);
}

#endif
