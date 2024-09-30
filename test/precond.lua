print("load script")
ug_load_script("plugins/plugin_Util/lua/lua-include.lua")
InitUG(2, AlgebraType("CPU", 1))
local iluDesc ={
	type = "ilu",
	damping = 0.1,
	overlap = true
}
local gsDesc ={
	type = "gs",
	consistentInterfaces = false,
	overlap = false
}
local sgsDesc = {
	type = "sgs",
	consistentInterfaces = true,
	overlap = false
}
local cgsDesc = {
	type = "cgs"
}
local solverutil = SolverUtil()
print("calling CreatePreconditioner (ilu)")
local precond =util.test.CreatePreconditioner(iluDesc, SolverUtil)
print("ilu"..precond:config_string())
local precond =util.test.CreatePreconditioner(gsDesc, SolverUtil)
print("gs:"..precond:config_string())
local precond =util.test.CreatePreconditioner(sgsDesc, SolverUtil)
print("sgs:"..precond.config_string())
local precond =util.test.CreatePreconditioner(cgsDesc, SolverUtil)
print("cgs:"..precond.config_string())
local lineSearchDesc = {
                       	type			= "standard",
                       	maxSteps		= 5,
                       	lambdaStart		= 1,
                       	lambdaReduce		= 0.4,
                       	acceptBest 		= true,
                       	checkAll		= false,
                       	verbose         = true,
                       }
local ls = util.test.CreateLineSearch(lineSearchDesc)