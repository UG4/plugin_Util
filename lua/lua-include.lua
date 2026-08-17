ug_load_script("ug_util.lua")
ug_load_script("externals/JSONForUG4/json-lua/json.lua")
util = util or {}
util.test = util.test or {}


function util.test.CreateConvCheck(convCheckDesc, solverutil)
    print("\n[lua-include.lua] ### create ConvCheck ### \n")
    -- create JSON object
    local jschema=JSON();
    -- create JSON validator
    validator=JSONSchemaValidator();
    -- load schema
    local path = ug_get_root_path().."/plugins/plugin_Util/json_schemata/ugcore/solver/convcheck.json"
    print("[lua-include.lua] schema path: "..path)
    JSON_load_from_file(jschema, path);
    validator:init(jschema)

    print(JSON_dump(jschema))

    local descriptor=JSON();
    local json_util = JSON();
    JSON_parse(descriptor, util.json.encode(convCheckDesc))
    print("|"..JSON_dump(descriptor).."|\n")
    validator:validate_document(descriptor)


    local convCheck = CreateConvCheckCPU1(descriptor)
    
    print(convCheck:config_string())
    return convCheck
end

function util.test.CreateLinearSolver(descriptor, solverutil)
    print("\n[lua-include.lua] ### create LinearSolver ### \n")
    local jsondesc = JSON()
    local json_util = JSON()
    JSON_parse(jsondesc, util.json.encode(descriptor))
    
    print("[lua-include.lua] function provider_CreateLinearSolver")
    local functionProvider = SolverUtilFunctionProvider()
    local solver = functionProvider:CreateLinearSolver(jsondesc, solverutil)
    print(solver:config_string())
    return solver
end

function util.test.CreateNewtonSolver(descriptor, solverutil)
    print("\n[lua-include.lua] ### create NewtonSolver ### \n")
    local jsondesc = JSON()
    local json_util = JSON()
    JSON_parse(jsondesc, util.json.encode(descriptor))
    print("[lua-include.lua] function provider_CreateNewtonSolver")
    local functionProvider = SolverUtilFunctionProvider()
    local solver = functionProvider:CreateNewtonSolver(jsondesc, solverutil)
    print(solver:config_string())
    return solver
end

function util.test.CreateSolver(descriptor, solverutil)
    print("\n[lua-include.lua] ### create Solver ### \n")
    if solverDesc.type == "newton" then
        solver = util.test.CreateNewtonSolver(solverDesc, solverutil)
    else
        solver = util.test.CreateLinearSolver(solverDesc, solverutil)
    end
    return solver
end

function util.test.CreatePreconditioner(desc, solverutil)
    print("\n[lua-include.lua] ### create Preconditioner ###")
    print("[lua-include.lua] Lua descriptor:")

    -- Iterate over all entries of the Lua descriptor.
    for key, value in pairs(desc) do
        if type(value) == "table" then
            print("                 " .. key .. " = {")
            for subKey, subValue in pairs(value) do
                print("                     " .. subKey .. " = " .. tostring(subValue))
            end
            print("                 }")
        else
            print("                 " .. key .. " = " .. tostring(value))
        end
    end

    -- Create the UG4 JSON object.
    local jsondesc = JSON()

    -- Convert the Lua table to JSON and parse it into the UG4 JSON object.
    JSON_parse(jsondesc, util.json.encode(desc))

    -- Print the JSON representation that will be passed to the C++ side.
    print("\n[lua-include.lua] JSON passed to C++:")
    print("                 " .. JSON_dump(jsondesc))

    -- Print a message before creating the function provider.
    print("\n[lua-include.lua] create function provider")

    -- Create the SolverUtil function provider.
    local functionProvider = SolverUtilFunctionProvider()

    -- Create the preconditioner from the JSON descriptor.
    local precond = functionProvider:CreatePreconditioner(jsondesc, solverutil)

    -- Print the resulting UG4 configuration.
    print("\n[lua-include.lua] resulting configuration:")
    print("  " .. precond:config_string())

    -- Return the created preconditioner.
    return precond
end

function util.test.CreateLineSearch(desc)
    print("\n[lua-include.lua] ### create LineSearch ### \n")
    local jsondesc = JSON()
    --parse descriptor as JSON
    JSON_parse(jsondesc, util.json.encode(desc))
    print("                    |"..JSON_dump(jsondesc).."|\n")
    --call our registered c++ function
    local functionProvider = SolverUtilFunctionProvider()
    ls = functionProvider:CreateLineSearch(jsondesc)
    --print(ls:config_string())
    --return linsearch
    return ls
end
