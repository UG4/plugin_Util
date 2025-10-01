print("load script")
ug_load_script("ug_util.lua")
ug_load_script("../lua/lua-include.lua")
InitUG(2, AlgebraType("CPU", 1))
dom = util.CreateDomain("grid.ugx", 1, {})
approxSpace = ApproximationSpace(dom)
approxSpace:print_statistic()
--local iluDesc = {
--	type = "ilu",
--	damping = 0.1,
--	overlap = true
--}
--local jacDesc = {
--	type = "jac",
--	damping = 0.66
--}
--local gsDesc = {
--	type = "gs",
--	consistentInterfaces = false,
--	overlap = false
--}
--local sgsDesc = {
--	type = "sgs",
--	consistentInterfaces = false,
--	overlap = false
--}
--local egsDesc = {
--	type = "egs"
--}
--local gmgDesc = {
--	type = "gmg",
--	-- approxSpace = approxSpace,
--	smoother = {type = "sgs"},
--	preSmooth = 1,
--	postSmooth = 1, 
--	baseLevel = 0,
--	baseSolver = "lu",
--	rap = false,
--    cycle = "V",
--    "discretization" = null,
--    gatheredBaseSolverIfAmbiguous = false,
--    preSmooth = 3,
--    postSmooth = 3,
--    rap = false,
--    rim = false,
--    emulateFullRefined = false,
--    transfer = "std",
--    debug = false,
--    mgStats = null
--}
--local cgsDesc = {
--	type = "cgs",
--	alpha = 1.0,
--	beta = 1.0,
--	weights = false,
--	relax = 1.0	
--}
--local sscDesc = {
--        type = "ssc", 
--        damping = 1.0, 
--        vertex = {{"p"}, {"ux", "uy"}}
--}
--local schurDesc = {
--	type = "schur",
--	dirichletSolver	= "lu",
--	skeletonSolver	= "lu"	
--}

--print(type(gmgDesc["approxSpace"]))
print("----- Create SolverUtil -----")
local solverutil = SolverUtil()
solverutil:setApproximationSpace("approxSpace", approxSpace)
--local solverutil = SolverUtil()
--solverutil:setDebug("debug", myDebugWriter)
print("----- Testing NewtonSolver with DebugWriter -----")

local newtonDesc = {
    type = "newton",
    debug = { 
        debug = true, 
        vtk = true, 
        conn_viewer = false 
    }
}
local newtonSolver = util.test.CreateNewtonSolver(newtonDesc, solverutil)
print("NewtonSolver created!")

--local precond = util.test.CreatePreconditioner(iluDesc, solverutil)
--print("ilu: "..precond:config_string())
--local precond = util.test.CreatePreconditioner(jacDesc, solverutil)
--print("jac: "..precond:config_string())
--local precond = util.test.CreatePreconditioner(gsDesc, solverutil)
--print("gs: "..precond:config_string())
--local precond = util.test.CreatePreconditioner(sgsDesc, solverutil)
--print("sgs: "..precond:config_string())
--local precond = util.test.CreatePreconditioner(egsDesc, solverutil)
--print("egs: "..precond:config_string())
--local precond = util.test.CreatePreconditioner(cgsDesc, solverutil)
--print("cgs: "..precond:config_string())
--local precond = util.test.CreatePreconditioner(gmgDesc,solverutil)
--print("gmg: "..precond:config_string())
--local precond = util.test.CreatePreconditioner(sscDesc,solverutil)
--print("ssc: "..precond:config_string())
--local precond = util.test.CreatePreconditioner(schurDesc,solverutil)
--print("schur: "..precond:config_string())
--local lineSearchDesc = {
--                       	type			= "standard",
--                       	maxSteps		= 5,
--                       	lambdaStart		= 1,
--                       	lambdaReduce		= 0.4,
--                       	acceptBest 		= true,
--                       	checkAll		= false,
--                       	verbose         = true,
--                       }
--local ls = util.test.CreateLineSearch(lineSearchDesc)

