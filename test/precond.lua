print("load script")
ug_load_script("ug_util.lua")
ug_load_script("../lua/lua-include.lua")
InitUG(2, AlgebraType("CPU", 1))
local iluDesc ={
	type = "ilu",
	damping = 0.1,
	overlap = true
}
local jacDesc ={
	type = "jac",
	damping = 0.66
}
local gsDesc ={
	type = "gs",
	consistentInterfaces = false,
	overlap = false
}
local sgsDesc = {
	type = "sgs",
	consistentInterfaces = false,
	overlap = false
}
local egsDesc = {
	type = "egs"
}

dom = util.CreateDomain("grid.ugx", 1, {})

local gmgDesc = {
	type = "gmg",
	approxSpace = approxSpace,
	smoother = {type = "sgs"},
	preSmooth = 1,
	postSmooth = 1, 
	baseLevel = 0,
	baseSolver = "lu",
	rap = false
}

--local cgsDesc = {
--	type = "cgs"
--}
local solverutil = SolverUtil()
print("calling CreatePreconditioner (ilu)")
local precond = util.test.CreatePreconditioner(iluDesc, SolverUtil)
print("ilu"..precond:config_string())
local precond = util.test.CreatePreconditioner(jacDesc, SolverUtil)
print("jac:"..precond:config_string())
local precond = util.test.CreatePreconditioner(gsDesc, SolverUtil)
print("gs:"..precond:config_string())
local precond = util.test.CreatePreconditioner(sgsDesc, SolverUtil)
print("sgs:"..precond:config_string())
local precond = util.test.CreatePreconditioner(egsDesc, SolverUtil)
print("egs:"..precond:config_string())
local precond = util.test.CreatePreconditioner(gmgDesc, SolverUtil)
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