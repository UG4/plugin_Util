print("load script")
ug_load_script("plugins/plugin_Util/lua/lua-include.lua")
InitUG(2, AlgebraType("CPU", 1))
local precondDesc ={
	type = "ilu",
	damping = 0.1,
	overlap = true
}
local solverutil = SolverUtil()
print("calling CreatePreconditioner")
local precond =util.test.CreatePreconditioner(precondDesc, SolverUtil)
print(precond:config_string())

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