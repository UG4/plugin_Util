print("[precond_test] load script")
ug_load_script("ug_util.lua")

ug_load_script("../lua/lua-include.lua")
InitUG(2, AlgebraType("CPU", 1))
local iluDesc ={
	type = "ilu",
	damping = 0.3,
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


dom = util.CreateDomain("grids/grid.ugx", 1, {})
approxSpace = ApproximationSpace(dom)
local gmgDesc = {
	type = "gmg",
	-- approxSpace = approxSpace,
	smoother = {type = "sgs"},
	baseLevel = 0,
	baseSolver = "lu",
	rap = false,
        cycle = "V",
        --"discretization" = null,
        gatheredBaseSolverIfAmbiguous = false,
        preSmooth = 3,
        postSmooth = 3,
        rap = false,
        rim = false,
        emulateFullRefined = false,
        transfer = "std",
        debug = false,
        mgStats = null
}

local cgsDesc = {
	type = "cgs",
	alpha = 1.0,
	beta = 1.0,
	weights = false,
	relax = 1.0	
}

approxSpace:print_statistic()

print("\n[>>>precond_test<<<] type(gmgDesc['approxSpace']): ", type(gmgDesc["approxSpace"]), "\n")

print("[>>> precond_test <<<] \ngmgDesc: ")

for key, value in pairs(gmgDesc) do
    if type(value)  == "table" then
        print("     " .. key .. " = {")
        for subKey, subValue in pairs(value) do
            print("        " .. subKey .. " = " .. tostring(subValue))
        end
        print("     }")
    else
        print("     " .. key .. " = " .. tostring(value))
    end
end


local solverutil = SolverUtil()
print("\n[>>>precond_test<<<] initialize solverutil ..")
solverutil:setApproximationSpace("approxSpace",approxSpace)

print("\n[>>>precond_test<<<] calling CreatePreconditioner")
local precond = util.test.CreatePreconditioner(iluDesc, solverutil)
print("\n[>>>precond_test<<<] ilu(config_string):\n"..precond:config_string())

local precond = util.test.CreatePreconditioner(jacDesc, solverutil)
print("\n[>>>precond_test<<<] jac(config_string):\n"..precond:config_string())

local precond = util.test.CreatePreconditioner(gsDesc, solverutil)
print("\n[>>>precond_test<<<] gs(config_string):\n"..precond:config_string())

local precond = util.test.CreatePreconditioner(sgsDesc, solverutil)
print("\n[>>>precond_test<<<] sgs(config_string):\n"..precond:config_string())

local precond = util.test.CreatePreconditioner(egsDesc, solverutil)
print("\n[>>>precond_test<<<] egs(config_string):\n"..precond:config_string())

local precond = util.test.CreatePreconditioner(cgsDesc, solverutil)
print("\n[>>>precond_test<<<] cgs(config_string):\n"..precond:config_string())

local precond = util.test.CreatePreconditioner(gmgDesc,solverutil)
print("\n[>>>precond_test<<<] gmg:(config_string):\n"..precond:config_string())

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
