**port CreateLineSearch LUA -> C++**

```plaintext
CreateSolver
├── PrepareSolverUtil
├── util.tableDesc.IsPreset
├── util.tableDesc.ToNameAndDesc
├── util.solver.defaults
├── NewtonSolver (through the creation of a "newton" solver)
│   ├── set_linear_solver
│   │   └── CreateLinearSolver
│   │       ├── PrepareSolverUtil
│   │       ├── util.tableDesc.IsPreset
│   │       ├── util.tableDesc.ToNameAndDesc
│   │       ├── util.solver.defaults
│   │       ├── LinearSolver, CG, BiCGStab, GMRES, LU, SuperLU, AgglomeratingSolver (conditional on the type of linear solver)
│   │       ├── set_preconditioner
│   │       │   └── CreatePreconditioner
│   │       │       ├── PrepareSolverUtil
│   │       │       ├── util.tableDesc.IsPreset
│   │       │       ├── util.tableDesc.ToNameAndDesc
│   │       │       ├── util.solver.defaults
│   │       │       ├── various preconditioners (e.g., ILU, ILUT, Jacobi, BlockGaussSeidel, etc.)
│   │       │       └── CreateOrdering
│   │       │           └── NativeCuthillMcKeeOrdering, BoostCuthillMcKeeOrdering, LexOrdering, etc. (depending on the specific ordering)
│   │       ├── set_convergence_check
│   │       │   └── CreateConvCheck
│   │       │       ├── PrepareSolverUtil
│   │       │       ├── util.tableDesc.IsPreset
│   │       │       ├── util.tableDesc.ToNameAndDesc
│   │       │       ├── util.solver.defaults
│   │       │       └── ConvCheck, CompositeConvCheck (depending on the convergence check descriptor)
│   │       └── SetDebugWriter
│   ├── set_convergence_check
│   │   └── CreateConvCheck
│   └── set_line_search
│       └── CreateLineSearch
│           ├── util.tableDesc.IsPreset
│           ├── util.tableDesc.ToNameAndDesc
│           ├── util.solver.defaults
│           └── StandardLineSearch (depending on the line search descriptor)
└── CreateLinearSolver (fallback if name is not "newton")
```plaintext
```

In diesem Abhängigkeitsbaum wollen wir das Blatt "CreateLineSeach" in C++ übersetzen. Die Lua-Funktion finden wir in der Datei `solver_util.lua`, standardmaessig unter dem Pfad `ug4/ugcore/scripts/util/solver_util.lua`
```markdown
```lua
function util.solver.CreateLineSearch(lineSearchDesc)
	if util.tableDesc.IsPreset(lineSearchDesc) then return lineSearchDesc end

	local name, desc = util.tableDesc.ToNameAndDesc(lineSearchDesc)
	local defaults	 = util.solver.defaults.lineSearch[name]
	if desc == nil then desc = defaults end
	
	local ls = nil
	if name == "standard" then
	--	battle booleans
		if desc.acceptBest == nil then desc.acceptBest = defaults.acceptBest end
		if desc.checkAll == nil then desc.checkAll = defaults.checkAll end
		
		ls = StandardLineSearch(desc.maxSteps or defaults.maxSteps,
								desc.lambdaStart or defaults.lambdaStart,
								desc.lambdaReduce or defaults.lambdaReduce,
								desc.acceptBest,
								desc.checkAll)
		if desc.verbose ~= nil then ls:set_verbose(desc.verbose) end
		if desc.suffDesc ~= nil then ls:set_suff_descent_factor(desc.suffDesc) end
		if desc.maxDefect ~= nil then ls:set_maximum_defect(desc.maxDefect) end
	end
	
	util.solver.CondAbort(ls == nil, "Invalid line-search specified: " .. name)
	if desc then
		desc.instance = ls
	end

	return ls
end
```

als Eingabe in die lua-Funktion bekommen wir den lineSearch descriptor (lua - Tabelle), wo Parameter der lineSearch definiert werden. Dabei ist darauf zu achten, dass default Werte fuer nicht definierte Parameter verwendet werden sollen. Default Werte fuer Loeserkonfigurationen koennen der Datei  `ug4/ugcore/scripts/util/solver.defaults.json` entnommen werden. Bei genauer Betrachtung der Funktion CreateLineSearch sehen wir, dass ein Object ls zurückgegeben wird. Dies wird durch den Aufruf ls = StandardLineSearch(...) erzeugt. Hierbei handelt es sich um einen Konstruktor eines C++ Objects, welches lua zur Verfuegung gestellt wurde. Als erstes sollten wir also diese Klasse im C++ - Code finden. Dafuer empfiehlt sich eine angemessene Entwicklerumgebung mit denen man das UG4 Projekt durchsuchen kann. Genauer suchen wir nach der *.cpp* Datei (meistens mit `bridge` im Namen...), welche die StandardLineSearch Klasse fuer lua registriert. Wenn wir UG4 nach `StandardLineSearch` durchsuchen, finden wir die Datei `algebra_bridge.cpp`
```cpp
//	StandardLineSearch
	{
		std::string grp = parentGroup; grp.append("/Discretization/Nonlinear");
		typedef StandardLineSearch<vector_type> T;
		typedef ILineSearch<vector_type> TBase;
		string name = string("StandardLineSearch").append(suffix);
		reg.add_class_<T, TBase>(name, grp)
			.add_constructor()
			.template add_constructor<void (*)(int, number, number, bool)>()
			.template add_constructor<void (*)(int, number, number, bool, bool)>()
			.add_method("set_maximum_steps", &T::set_maximum_steps, "Max. steps of the line search", "steps")
			.add_method("set_lambda_start", &T::set_lambda_start, "Initial step scaling factor", "start")
			.add_method("set_reduce_factor", &T::set_reduce_factor, "Step reduction factor", "factor")
			.add_method("set_suff_descent_factor", &T::set_suff_descent_factor, "Sufficient descent", "factor")
			.add_method("set_accept_best", &T::set_accept_best, "If to accept the best", "bAcceptBest")
			.add_method("set_check_all", &T::set_check_all, "If to test all the steps", "bCheckAll")
			.add_method("set_maximum_defect", &T::set_maximum_defect, "Max. norm of the defect", "maxDef")
			.add_method("set_verbose", &T::set_verbose, "Verbosity level", "verboseLevel")
			.add_method("set_offset", &T::set_offset, "Indentation", "strOffset")
			.set_construct_as_smart_pointer(true);
		reg.add_class_to_group(name, "StandardLineSearch", tag);
	}
```
Wir sehen, dass durch den lua - Aufruf StandardLineSearch die templated C++ - Klasse `StandardLineSearch<vector_type>` aufgerufen wird. Das Templateargument beschreibt dabei die verwendete Algebra und hat uns eigentlich wenig zu interessieren. Die Zeile `set_construct_as_smartpointer(true)` sagt uns, 
dass ein SmartPointer auf das isntanziierte Objekt zurueckgegeben wird. Mit der Eingabe (JSON Deskriptor)  und Ausgabe (`SmartPtr<StandardLineSearch<vector_type>`) koennen wir nun das geruest unserer C++ Funktion erstellen:
```cpp

template<typename TAlgebra>
    SmartPtr<StandardLineSearch<typename TAlgebra::vector_type>> CreateLineSearch(nlohmann::json& desc){
        // typedef for convenience
        typedef StandardLineSearch<typename TAlgebra::vector_type>      line_search_type;
        SmartPtr<line_search_type> ls;
        
        // TODO 
        
        return ls;
    }
```
Nun sollten wir uns um die interne Logik der Funktion kuemmern. Wir gehen dafuer Schritt fuer Schritt durch den Deskriptor, setzen default - Werte und gegebenenfalls nicht-default Werte. Ein moeglicher Deskriptor sieht wie folgt aus:
```
{
	type			= "standard",
	maxSteps		= 10,
	lambdaStart		= 1,
	lambdaReduce		= 0.5,
	acceptBest 		= true,
	checkAll		= false
}
```
Als erstes empfiehlt es sich die default Werte zu laden: 
```cpp

nlohmann::json json_default_preconds = json_predefined_defaults::solvers.at("lineSearch");

```
Nun Fragen wir sukzessive alle Eintraege des Deskriptors ab, wobei wir zunaechst den default Wert laden und anschliessend gegebenenfalls ueberschreiben. :
```cpp

// handle type of line search
        // default
        std::string type = "standard";
        // input type
        if(desc.contains("type")){
            type = desc["type"];
        }

        if(type == "standard"){
            // handle parameters of standard line search
            int maxSteps = json_default_lineSearch[type]["maxSteps"];
            if(desc.contains("maxSteps")){
                maxSteps = desc["maxSteps"];
            }
            number lambdaStart = json_default_lineSearch[type]["lambdaStart"];
            if(desc.contains("lambdaStart")){
                lambdaStart = desc["lambdaStart"];
            }
            number lambdaReduce = json_default_lineSearch[type]["lambdaReduce"];
            if(desc.contains("lambdaReduce")){
                lambdaReduce = desc["lambdaReduce"];
            }
            bool acceptBest = json_default_lineSearch[type]["acceptBest"];
            if(desc.contains("acceptBest")){
                acceptBest = desc["acceptBest"];
            }
            bool checkAll = json_default_lineSearch[type]["checkAll"];
            if(desc.contains("checkAll")){
                checkAll = desc["checkAll"];
            }

```
Anschliessend instanziieren wir die `StandardLineSearch` mit den gefundenen Parametern:
```cpp
// create line search with chosen parameters
            ls = make_sp(new line_search_type(maxSteps,
                                              lambdaStart,
                                              lambdaReduce,
                                              acceptBest,
                                              checkAll));
```
und fuegen optionale Parameter hinzu:
```cpp
// set conditional parameters
            if(desc.contains("verbose")){
                bool verbose = desc["verbose"];
                ls->set_verbose(verbose);
            }
            if(desc.contains("suffDesc")){
                number suffDesc = desc["suffDesc"];
                ls->set_suff_descent_factor(suffDesc);
            }
            if(desc.contains("maxDefect")){
                number maxDefect = desc["maxDefect"];
                ls->set_maximum_defect(maxDefect);
            }
```
unsere komplette Funktion sieht dann so aus:
```cpp
template<typename TAlgebra>
    SmartPtr<StandardLineSearch<typename TAlgebra::vector_type>> CreateLineSearch(nlohmann::json& desc){
        // typedef for convenience
        typedef StandardLineSearch<typename TAlgebra::vector_type>      line_search_type;
        SmartPtr<line_search_type> ls;

        // load defaults
        nlohmann::json json_default_lineSearch = json_predefined_defaults::solvers["lineSearch"];

        // handle type of line search
        // default
        std::string type = "standard";
        // input type
        if(desc.contains("type")){
            type = desc["type"];
        }

        if(type == "standard"){
            // handle parameters of standard line search
            int maxSteps = json_default_lineSearch[type]["maxSteps"];
            if(desc.contains("maxSteps")){
                maxSteps = desc["maxSteps"];
            }
            number lambdaStart = json_default_lineSearch[type]["lambdaStart"];
            if(desc.contains("lambdaStart")){
                lambdaStart = desc["lambdaStart"];
            }
            number lambdaReduce = json_default_lineSearch[type]["lambdaReduce"];
            if(desc.contains("lambdaReduce")){
                lambdaReduce = desc["lambdaReduce"];
            }
            bool acceptBest = json_default_lineSearch[type]["acceptBest"];
            if(desc.contains("acceptBest")){
                acceptBest = desc["acceptBest"];
            }
            bool checkAll = json_default_lineSearch[type]["checkAll"];
            if(desc.contains("checkAll")){
                checkAll = desc["checkAll"];
            }
            // create line search with chosen parameters
            ls = make_sp(new line_search_type(maxSteps,
                                              lambdaStart,
                                              lambdaReduce,
                                              acceptBest,
                                              checkAll));
            // set conditional parameters
            if(desc.contains("verbose")){
                bool verbose = desc["verbose"];
                ls->set_verbose(verbose);
            }
            if(desc.contains("suffDesc")){
                number suffDesc = desc["suffDesc"];
                ls->set_suff_descent_factor(suffDesc);
            }
            if(desc.contains("maxDefect")){
                number maxDefect = desc["maxDefect"];
                ls->set_maximum_defect(maxDefect);
            }

        }
        // force exit if line search is invalid
        CondAbort(ls.invalid(), "Invalid line-search specified: " + type);
        desc["instance"] = ls;
        return ls;
    }
```
**Rueckwaertskompatibilitaet**
Wenn moeglich soll unsere neue Implementation kompatibel mit alten lua-Skripten sein. ein Aufruf der alten lua Funktion in einem lua skript soll sich also genauso verhalten wie zuvor. Dies koennen wir realisieren, indem wir in der lua Funktion einfach unsere neue C++ Funktion aufrufen. zunaechst legen wir dafuer lua Funktionen im Verzeichnis `plugin_Util/lua/lua-include` ab. Wenn alles fertig portiert ist, sollten diese in den Kern von UG4 gezogen werden. Damit unsere neue C++ - Funktion von lua aus aufgerufen werden kann, muessen wir diese fuer lua registrieren. Das geschieht in der Datei `plugin_util.cpp`.
Eine Hilfsklasse `SolverUtilFunctionProvider` wurde zur Verfuegung gestellt und bereits registriert. Wenn eine neue C++ Funktion implementiert wurde, fuegen wir diese dem function provider hinzu:
```cpp
template<typename TDomain, typename TAlgebra>
class SolverUtilFunctionProvider{
public:
    typedef typename TAlgebra::vector_type              vector_type;
    SolverUtilFunctionProvider(){};
	// INSERT CALL TO NEW FUNCTION HERE
	SmartPtr<StandardLineSearch<vector_type>> GetCreateLineSearch(nlohmann::json &desc){
        	return CreateLineSearch<TAlgebra>(desc);
    	}
    	// END INSERT
};
```
anschliessend muessen wir die hinzugefuegte methode im registry bekannt machen (datei `util_plugin.cpp`):
```cpp
{
		//FUNCTION PROVIDER GETS REGISTERED HERE
            typedef SolverUtilFunctionProvider<TDomain, TAlgebra> T;
            typedef typename TAlgebra::vector_type vector_type;
            string name = string("SolverUtilFunctionProvider").append(suffix);
            reg.template add_class_<T>(name, grp).add_constructor()
                    // INSERT NEW METHOD HERE 
                    .add_method("CreateLineSearch", &T::GetCreateLineSearch)
                    // END INSERT
                    .set_construct_as_smart_pointer(true);
            reg.add_class_to_group(name, "SolverUtilFunctionProvider", tag);
        }
```
Nun koennen wir `CreateLineSearch` in lua wie folgt aufrufen:  
(1) kovertiere deskriptor (lua tabelle) in json format:
```lua
-- we assume our lua table descriptor is desc
local jsondesc = JSON()
JSON_parse(jsondesc, util.json.encode(desc)
-- now we have a json descriptor jsondesc
```
(2) instanziiere function provider
```lua
-- create function provider
local funcProv = SolverUtilFunctionProvider()
```
(3) rufe usnere regsitrierte funktion `CreateLineSearch auf`
```lua
-- call CreateSolver
local ls = funcProv:CreateLineSearch(jsondesc)

``` 
insgesamt sieht unsere neue, rueckwaertskompatible lua Funktion also wie folgt aus:
```lua

function util.test.CreateLineSearch(desc)

    local jsondesc = JSON()
    
    --parse descriptor as JSON
    JSON_parse(jsondesc, util.json.encode(desc))
    
    -- create function provider
    local functionProvider = SolverUtilFunctionProvider()
    
    --call our registered c++ function
    local ls = functionProvider:CreateLineSearch(jsondesc)
    print(ls:config_string())
    --return linsearch
    return ls
end

```