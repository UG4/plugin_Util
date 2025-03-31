
ug_load_script("ug_util.lua")
ug_load_script("externals/JSONForUG4/json-lua/json.lua")
util = util or {}
util.test = util.test or {}

function util.test.CreateConvCheck(convCheckDesc, solverutil)

    -- create JSON object
    local jschema=JSON();
    -- create JSON validator
    validator=JSONSchemaValidator();
    -- load schema
    local path = ug_get_root_path().."/plugins/plugin_Util/json_schemata/ugcore/solver/convcheck.json"
    print("schema path: "..path)
    JSON_load_from_file(jschema, path);
    validator:init(jschema)

    print(JSON_dump(jschema))

    local descriptor=JSON();
    local json_util = JSON();
    JSON_parse(descriptor, util.json.encode(convCheckDesc))
    JSON_parse(json_util, util.json.encode(solverutil))
    print("|"..JSON_dump(descriptor).."|")
    validator:validate_document(descriptor)


    local convCheck = CreateConvCheckCPU1(descriptor)
    print(convCheck:config_string())
    return convCheck
end

function util.test.CreateSolver(descriptor, solverutil)


end

function util.test.CreatePreconditioner(desc, solverutil)
    
    local jsondesc = JSON()
    JSON_parse(jsondesc, util.json.encode(desc))
    --print("|"..JSON_dump(jsondesc).."|")
    -- local solverut = SolverUtil()
    print("create function provider")
    local functionProvider = SolverUtilFunctionProvider()
    local precond = functionProvider:CreatePreconditioner(jsondesc, solverutil)
    print(precond:config_string())
    return precond
end

function util.test.CreateLineSearch(desc)

    local jsondesc = JSON()
    --parse descriptor as JSON
    JSON_parse(jsondesc, util.json.encode(desc))
    print("|"..JSON_dump(jsondesc).."|")
    --call our registered c++ function
    local functionProvider = SolverUtilFunctionProvider()
    ls = functionProvider:CreateLineSearch(jsondesc)
    print(ls:config_string())
    --return linsearch
    return ls
end