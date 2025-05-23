//
// Created by julian on 6/4/24.
//

#ifndef UG4_SOLVER_UTIL_H
#define UG4_SOLVER_UTIL_H

#include <nlohmann/json.hpp>

#include "common/common.h"
#include "common/log.h"
#include "common/util/smart_pointer.h"
#include "common/util/plugin_util.h"
#include <memory>
#include <unordered_map>
#include <variant>
#include <stdexcept>
#include "registry/registry.h"
#include "../../ugcore/ugbase/lib_disc/function_spaces/grid_function_util.h"
#include "lib_algebra/operator/convergence_check.h"
#include "lib_algebra/operator/linear_solver/linear_solver.h"
#include "lib_algebra/ordering_strategies/algorithms/boost_cuthill_mckee_ordering.cpp"
#include "lib_algebra/ordering_strategies/algorithms/boost_minimum_degree_ordering.cpp"
#include "lib_algebra/ordering_strategies/algorithms/SCC_ordering.cpp"
#include "lib_algebra/ordering_strategies/algorithms/topological_ordering.cpp"
#include "lib_disc/ordering_strategies/algorithms/lexorder.h"
#include "lib_disc/ordering_strategies/algorithms/riverorder.h"
// include solver components
#include "../../externals/JSONForUG4/json-cxx/include/nlohmann/adl_serializer.hpp"
#include "../../externals/JSONForUG4/json-cxx/include/nlohmann/json.hpp"
#include "../../externals/JSONForUG4/json-cxx/single_include/nlohmann/json.hpp"
#include "../../externals/JSONForUG4/json-cxx/single_include/nlohmann/json.hpp"
#include "lib_disc/function_spaces/grid_function.h"
#include "lib_algebra/operator/preconditioner/ilu.h"
#include "lib_algebra/operator/preconditioner/ilut.h"
#include "lib_algebra/operator/preconditioner/jacobi.h"
#include "lib_algebra/operator/preconditioner/block_gauss_seidel.h"
#include "lib_algebra/operator/preconditioner/gauss_seidel.h"
#include "lib_disc/operator/linear_operator/element_gauss_seidel/element_gauss_seidel.h"
#include "lib_disc/operator/linear_operator/element_gauss_seidel/component_gauss_seidel.h"
#include "lib_disc/operator/linear_operator/multi_grid_solver/mg_solver.h"
#include "lib_disc/operator/non_linear_operator/newton_solver/newton.h"
#include "lib_disc/operator/non_linear_operator/line_search.h"
#include "lib_algebra/operator/linear_solver/cg.h"
#include "lib_algebra/operator/linear_solver/bicgstab.h"
#include "lib_algebra/operator/linear_solver/gmres.h"
#include "lib_algebra/operator/linear_solver/agglomerating_solver.h"
#include "lib_algebra/operator/linear_solver/lu.h"
#include "../SuperLU6/super_lu.h"
#include "lib_disc/spatial_disc/domain_disc.h"

namespace ug
{
    namespace Util
    {

        void CondAbort(bool condition, std::string message)
        {
            UG_ASSERT(!condition, "ERROR in util.solver: " << message);
        }

        template <typename TDomain, typename TAlgebra>
        SmartPtr<IOrderingAlgorithm<TAlgebra, std::vector<size_t>>>
        CreateOrdering(nlohmann::json &desc);
        template <typename TAlgebra>
        SmartPtr<StdConvCheck<typename TAlgebra::vector_type>> CreateConvCheck(nlohmann::json &descriptor)
        {
            typedef typename TAlgebra::vector_type vector_type;
            // set parameters from descriptor attributes
            bool verbose = false;
            if (descriptor.contains("verbose")){
                verbose = descriptor["verbose"];
            }
            std::string name = "standard";
            if (descriptor.contains("type")){
                name = descriptor["type"];
            }
            number iterations = 100;
            if (descriptor.contains("iterations")){
                iterations = descriptor["iterations"];
            }
            number reduction = 1e-6;
            if (descriptor.contains("reduction")){
                reduction = descriptor["reduction"];
            }
            number absolute = 1e-12;
            if (descriptor.contains("absolute")){
                absolute = descriptor["absolute"];
            }
            bool suppress_unsuccessful = false;
            if (descriptor.contains("suppress_unsuccessful")){
                suppress_unsuccessful = descriptor["suppress_unsuccessful"];
            }

            // create convergence check
            SmartPtr<StdConvCheck<vector_type>> convCheck = make_sp<StdConvCheck<vector_type>>(
                new StdConvCheck<vector_type>(iterations,
                                              absolute,
                                              reduction,
                                              verbose,
                                              suppress_unsuccessful));
            return convCheck;
        }

        // Variant type for storing different solver components
        template<typename TDomain, typename TAlgebra>
        using SolverComponent = std::variant<
            SmartPtr<ApproximationSpace<TDomain>>,
            SmartPtr<ILU<TAlgebra>>,
            SmartPtr<ILUTPreconditioner<TAlgebra>>,
            SmartPtr<Jacobi<TAlgebra>>,
            SmartPtr<BlockGaussSeidel<TAlgebra, false, false>>,
            SmartPtr<BlockGaussSeidel<TAlgebra, true, false>>,
            SmartPtr<BlockGaussSeidel<TAlgebra, false, true>>,
            SmartPtr<BlockGaussSeidel<TAlgebra, true, true>>,
            SmartPtr<SymmetricGaussSeidel<TAlgebra>>,
            SmartPtr<GaussSeidel<TAlgebra>>,
            SmartPtr<ElementGaussSeidel<TDomain, TAlgebra>>,
            SmartPtr<ComponentGaussSeidel<TDomain, TAlgebra>>,
            SmartPtr<NewtonSolver<TAlgebra>>,
            SmartPtr<IAssemble<TAlgebra>>,
            SmartPtr<AssembledMultiGridCycle<TDomain, TAlgebra>>,
            SmartPtr<ug::GridFunctionDebugWriter<TDomain, TAlgebra>>,
            SmartPtr<LinearSolver<typename TAlgebra::vector_type>>,
            SmartPtr<CG<typename TAlgebra::vector_type>>,
            SmartPtr<BiCGStab<typename TAlgebra::vector_type>>,
            SmartPtr<GMRES<typename TAlgebra::vector_type>>,
            SmartPtr<SuperLUSolver<TAlgebra>>,
            SmartPtr<LU<TAlgebra>>,
            SmartPtr<AgglomeratingSolver<TAlgebra>>>;

        template<typename TDomain, typename TAlgebra>
        class SolverUtil{
            /// Container class for objects used in the CreateSolver routine.
            /*
             * components are stored as std::variant and populate an unordered
             * map for access. Calls to various subroutines are made
             * */
        public:
            void setDiscretization(const std::string key, SmartPtr<IAssemble<TAlgebra>> component){
                setComponent(key, component);
            }

            void setApproximationSpace(const std::string key, SmartPtr<ApproximationSpace<TDomain>> component){
                setComponent(key, component);
            }

            void setDebug(const std::string key, SmartPtr<ug::GridFunctionDebugWriter<TDomain, TAlgebra>> component){
                setComponent(key, component);
            }

            void setComponent(const std::string &key, SolverComponent<TDomain, TAlgebra> component){
                components[key] = component;
            }

            SolverComponent<TDomain, TAlgebra> getComponent(const std::string &key) const{
                auto it = components.find(key);
                if (it != components.end()){
                    return it->second;
                }
                UG_THROW("SolverUtil.getComponent(" << key << "): no component named " << key << " found");
            }

            bool hasComponent(const std::string &key) const{
                auto it = components.find(key);
                if (it != components.end()){
                    return true;
                }
                return false;
            }

            template<typename T>
            SmartPtr<T> getComponentAs(const std::string &key) const{
                return std::get<SmartPtr<T>>(getComponent(key));
            }

        private:
            std::unordered_map<std::string, SolverComponent<TDomain, TAlgebra>> components;
        };

        template<typename TAlgebra>
        SmartPtr<StandardLineSearch<typename TAlgebra::vector_type>> CreateLineSearch(nlohmann::json &desc){
            // typedef for convenience
            typedef StandardLineSearch<typename TAlgebra::vector_type> line_search_type;
            SmartPtr<line_search_type> ls;

            // load defaults
            nlohmann::json json_default_lineSearch = json_predefined_defaults::solvers["lineSearch"];

            // handle type of line search
            // default
            std::string type = "standard";
            // input type
            if (desc.contains("type")){
                type = desc["type"];
            }

            if (type == "standard"){
                // handle parameters of standard line search
                int maxSteps = json_default_lineSearch[type]["maxSteps"];
                if (desc.contains("maxSteps")){
                    maxSteps = desc["maxSteps"];
                }
                number lambdaStart = json_default_lineSearch[type]["lambdaStart"];
                if (desc.contains("lambdaStart")){
                    lambdaStart = desc["lambdaStart"];
                }
                number lambdaReduce = json_default_lineSearch[type]["lambdaReduce"];
                if (desc.contains("lambdaReduce")){
                    lambdaReduce = desc["lambdaReduce"];
                }
                bool acceptBest = json_default_lineSearch[type]["acceptBest"];
                if (desc.contains("acceptBest")){
                    acceptBest = desc["acceptBest"];
                }
                bool checkAll = json_default_lineSearch[type]["checkAll"];
                if (desc.contains("checkAll")){
                    checkAll = desc["checkAll"];
                }
                // create line search with chosen parameters
                ls = make_sp(new line_search_type(maxSteps,
                                                  lambdaStart,
                                                  lambdaReduce,
                                                  acceptBest,
                                                  checkAll));
                // set conditional parameters
                if (desc.contains("verbose")){
                    bool verbose = desc["verbose"];
                    ls->set_verbose(verbose);
                }
                if (desc.contains("suffDesc")){
                    number suffDesc = desc["suffDesc"];
                    ls->set_suff_descent_factor(suffDesc);
                }
                if (desc.contains("maxDefect")){
                    number maxDefect = desc["maxDefect"];
                    ls->set_maximum_defect(maxDefect);
                }
            }
            // force exit if line search is invalid
            CondAbort(ls.invalid(), "Invalid line-search specified: " + type);

            return ls;
        }

        template<typename TDomain, typename TAlgebra>
        SmartPtr<ILinearOperatorInverse<typename TAlgebra::vector_type>>
        //SmartPtr<ILinearIterator<typename TAlgebra::vector_type>>
        CreateLinearSolver(nlohmann::json &solverDesc, SolverUtil<TDomain, TAlgebra> &solverutil){
            typedef typename TAlgebra::vector_type TVector;
            typedef LinearSolver<TVector> TLinSolv;
            //typedef ILinearOperatorInverse<TVector> TLinSolv;
            typedef ILinearOperatorInverse<TVector> TReturn;
            //SmartPtr<TLinSolv> linSolver;
            SmartPtr<TReturn> linSolver;
            // load defaults
            nlohmann::json json_default_linearSolver = json_predefined_defaults::solvers["linearSolver"];

            std::string precondtype = "ilu";
            std::string convChecktype = "standard";

            bool createPrecond = false;
            bool createConvCheck = false;

            // if no descriptor given, create default linear solver
            // if (!solverDesc.contains("linearSolver")){
            // if (!solverDesc.contains("linearSolver")){
            //    UG_LOG("default LinearSolver\n")
            //    SmartPtr<ILU<TAlgebra>> ilu = make_sp<ILU<TAlgebra>>(new ILU<TAlgebra>());
            //    SmartPtr<StdConvCheck<TVector>> convCheck = make_sp<StdConvCheck<TVector>>(
            //        new StdConvCheck<TVector>(100, 1e-9, 1e-12));
            //    SmartPtr<TLinSolv> default_linear_solver = make_sp(new TLinSolv());
            //    default_linear_solver->set_convergence_check(convCheck);
            //    default_linear_solver->set_preconditioner(ilu);
            //    return default_linear_solver;
            //}

            std::string type = "linear";
            if (solverDesc.contains("type")){
                type = solverDesc["type"];
            }
            if (type == "linear"){
                UG_LOG("type solver is linear\n")
                // create linear
                auto linear = make_sp(new TLinSolv());
                linSolver = linear.template cast_dynamic<TReturn>();

                // configure linear
                precondtype = json_default_linearSolver["linear"]["precond"];
                convChecktype = json_default_linearSolver["linear"]["convCheck"];
                createPrecond = true;
                createConvCheck = true;
            }
            else if (type == "cg"){
                UG_LOG("type solver is cg\n")
                // create CGSolver
                typedef CG<TVector> TCGSolver;
                SmartPtr<TCGSolver> CG = make_sp(new TCGSolver());
                //linSolver = CG.template cast_dynamic<TLinSolv>();
                linSolver = CG.template cast_dynamic<TReturn>();
                // configure cg
                precondtype = json_default_linearSolver["cg"]["precond"];
                convChecktype = json_default_linearSolver["cg"]["convCheck"];
                createPrecond = true;
                createConvCheck = true;
            }
            else if (type == "bicgstab"){
                UG_LOG("type solver is bicgstab\n")
                // create BiCGStabSolver
                typedef BiCGStab<TVector> TBiCGStabSolver;
                SmartPtr<TBiCGStabSolver> BICGSTAB = make_sp(new TBiCGStabSolver());
                // linSolver = BICGSTAB.template cast_dynamic<TLinSolv>();
                linSolver = BICGSTAB.template cast_dynamic<TReturn>();
                // configure bicgstab
                precondtype = json_default_linearSolver["bicgstab"]["precond"];
                convChecktype = json_default_linearSolver["bicgstab"]["convCheck"];
                createPrecond = true;
                createConvCheck = true;
            }
            else if (type == "gmres"){
                UG_LOG("type solver is gmres\n")
                // create gmres
                typedef GMRES<TVector> TGMRESSolver;
                SmartPtr<TGMRESSolver> GMRES;
                int restart = 5;
                if (json_default_linearSolver.contains("gmres") &&
                    json_default_linearSolver["gmres"].contains("restart")){
                    restart = json_default_linearSolver["gmres"]["restart"];
                }
                if (solverDesc.contains("restart")){
                    restart = solverDesc["restart"];
                }
                GMRES = make_sp(new TGMRESSolver(restart));
                //linSolver = GMRES.template cast_dynamic<TLinSolv>();
                linSolver = GMRES.template cast_dynamic<TReturn>();
                createPrecond = true;
                createConvCheck = true;
            }
            else if (type == "lu"){
                UG_LOG("CreateLinearSolver, type solver is lu \n")
                // AgglomeratingSolver

            #ifdef UG_USE_SUPERLU
                UG_LOG("CreateLinearSolver,name_lu, AgglomeratingSolver(SuperLU()) \n");
                auto SuperLU = make_sp(new AgglomeratingSolver<TAlgebra>(make_sp(new SuperLUSolver<TAlgebra>())));
                // linSolver = make_sp(new AgglomeratingSolver<TAlgebra>(make_sp(new SuperLUSolver<TAlgebra>())))
                //            .template cast_dynamic<TLinSolv>();
                auto test = SuperLU.template cast_dynamic<ILinearOperatorInverse<typename TAlgebra::vector_type>>();
                linSolver = SuperLU.template cast_dynamic<TReturn>();

            #else
                UG_LOG("CreateLinearSolver,name_lu, AgglomeratingSolver(LU_solver) \n")
                typedef LU<TAlgebra> TLU;
                SmartPtr<TLU> LU_solver = make_sp(new TLU());

                bool showProgress = json_default_linearSolver["lu"]["showProgress"];
                if (solverDesc.contains("showProgress") && !solverDesc["showProgress"].is_null()){
                    showProgress = solverDesc["showProgress"];
                }
                LU_solver->set_show_progress(showProgress);

                bool info = json_default_linearSolver["lu"]["info"];
                if (solverDesc.contains("info")){
                    info = solverDesc["info"];
                }
                LU_solver->set_info(info);
                linSolver = make_sp(new AgglomeratingSolver<TAlgebra>(LU_solver))
                            .template cast_dynamic<TLinSolv>();

            #endif
            }
            else if (type == "uglu"){
                UG_LOG("CreateLinearSolver,name_uglu, AgglomeratingSolver(LU_solver) \n")
                // create LU Solver
                typedef LU<TAlgebra> TLU;
                SmartPtr<TLU> LU_solver = make_sp(new TLU());

                bool showProgress = json_default_linearSolver["uglu"]["showProgress"];
                if (solverDesc.contains("showProgress") && !solverDesc["showProgress"].is_null()){
                    showProgress = solverDesc["showProgress"];
                }
                LU_solver->set_show_progress(showProgress);

                bool info = json_default_linearSolver["uglu"]["info"];
                if (solverDesc.contains("info")){
                    info = solverDesc["info"];
                }
                LU_solver->set_info(info);

                // create AgglomeratingSolver(LU_Solver))
                linSolver = make_sp(new AgglomeratingSolver<TAlgebra>(LU_solver)).template cast_dynamic<TReturn>();
            }
            else if (type == "superlu"){
                UG_LOG("CreateLinearSolver,name_superlu, AgglomeratingSolver(SuperLU()) \n")
                // create SuperLU Solver
                typedef SuperLUSolver<TAlgebra> TSupLUSolv;
                SmartPtr<TSupLUSolv> superlu = make_sp(new TSupLUSolv());

                // typedef for convenience
                typedef typename TAlgebra::vector_type Tvector;
                typedef ILinearOperatorInverse<Tvector, Tvector> TSolverBase;

                // cast to ILinearOperatorInverse<Tvector, Tvector>
                SmartPtr<TSolverBase> linOpInverse;
                linOpInverse = superlu.template cast_static<TSolverBase>();

                // create AgglomeratingSolver(SuperLU_Solver))
                typedef AgglomeratingSolver<TAlgebra> TAggSolver;
                SmartPtr<TAggSolver> superLU = make_sp(new TAggSolver(linOpInverse));
                linSolver = superLU.template cast_dynamic<TReturn>();
            }
            
            // Checks for a valid solver
            CondAbort(linSolver.invalid(), "Invalid linear solver specified: " + type);

            // Checks if solverDesc has precond and overwrite default
            if (solverDesc.contains("precond") && solverDesc["precond"].is_string()){
                precondtype = solverDesc["precond"];
            }
            // Checks if solverDesc has convCheck and overwrite default
            if (solverDesc.contains("convCheck") && solverDesc["convCheck"].is_string()){
                convChecktype = solverDesc["convCheck"];
            }

            // Preconditioner set on linear, cd, bicgstab, gmres
            if (createPrecond){
                nlohmann::json precondDesc;

                if (solverDesc.contains("precond") && solverDesc["precond"].is_object()){
                    UG_LOG("CreateLinearSolver, solverDesc is object \n")
                    precondDesc = solverDesc["precond"];
                }
                else{
                    UG_LOG("CreateLinearSolver, solverDesc is string \n")
                    precondDesc["type"] = precondtype;
                }
                UG_LOG("CreateLinearSolver, precondDesc: \n")
                UG_LOG(precondDesc.dump(4) << "\n")
                SmartPtr<ILinearIterator<typename TAlgebra::vector_type>> preconditioner = CreatePreconditioner(
                    precondDesc, solverutil);
                (linSolver.template cast_dynamic<IPreconditionedLinearOperatorInverse<typename TAlgebra::vector_type>>())->set_preconditioner(preconditioner);
            }

            // Convergence Check set
            if (createConvCheck){
                nlohmann::json convCheckDesc;

                if (solverDesc.contains("convCheck") && solverDesc["convCheck"].is_object()){
                    convCheckDesc = solverDesc["convCheck"];
                }
                else{
                    convCheckDesc["type"] = convChecktype;
                }
                UG_LOG("CreateLinearSolver, convCheckDesc: \n")
                UG_LOG(convCheckDesc.dump(4) << "\n")
                SmartPtr<StdConvCheck<TVector>> convCheck = CreateConvCheck<TAlgebra>(convCheckDesc);
                linSolver->set_convergence_check(convCheck);
            }

            // TODO: SetDebugWriter(linSolver, solverDesc, defaults, solverutil)

            return linSolver; //.template cast_dynamic<TReturn>();
        }


        template<typename TDomain, typename TAlgebra>
        SmartPtr<NewtonSolver<TAlgebra>>
        CreateNewtonSolver(nlohmann::json &solverDesc, SolverUtil<TDomain, TAlgebra> &solverutil){
            UG_LOG("Creating NewtonSolver\n")
            std::cout << "[DEBUG] JSON input to CreateNewtonSolver:\n" << solverDesc.dump(4) << std::endl;
            nlohmann::json json_default_nonlinearSolver = json_predefined_defaults::solvers["nonlinearSolver"];

            // PrepareSolverUtil<TDomain, TAlgebra>();
            std::string type = "linear";
            if (solverDesc.contains("type")){
                type = solverDesc["type"];
            }
            if (type == "newton"){
                UG_LOG("CreateNewtonSolver, type solver is newton \n")

                auto newtonSolver = make_sp(new NewtonSolver<TAlgebra>());

                // linear solver
                nlohmann::json linSolverDesc;

                // get descriptor for linear solver
                if (solverDesc.contains("linSolver")){
                    if (solverDesc["linSolver"].is_string()){
                        UG_LOG("CreateNewtonSolver, linSolver is a string: \n");
                        // linSolverType = solverDesc["linSolver"];
                        linSolverDesc["type"] = solverDesc["linSolver"];
                    }
                    else if (solverDesc["linSolver"].is_object()){
                        UG_LOG("CreateNewtonSolver, linSolver is an object \n");
                        linSolverDesc = solverDesc["linSolver"];
                    }
                    else{
                        UG_LOG(">> CreateNewtonSolver, linSolver must be a string or an object \n")
                    }
                }
                else{
                    // use default
                    if (json_default_nonlinearSolver["newton"]["linSolver"].is_string()){
                        linSolverDesc["type"] = json_default_nonlinearSolver["newton"]["linSolver"];
                    }
                    else if (json_default_nonlinearSolver["newton"]["linSolver"].is_object()){
                        linSolverDesc = json_default_nonlinearSolver["newton"]["linSolver"];
                    }
                }
                UG_LOG("CreateNewtonSolver, linSolverDesc:\n");
                UG_LOG(linSolverDesc.dump(4) << "\n");

                SmartPtr<ILinearOperatorInverse<typename TAlgebra::vector_type>> linSolver;
                linSolver = CreateLinearSolver(linSolverDesc, solverutil);
                newtonSolver->set_linear_solver(linSolver);

                // convergence check
                std::string convCheckType = json_default_nonlinearSolver["newton"]["convCheck"];
                // get descriptor for convergence check
                if (solverDesc.contains("convCheck")){
                    convCheckType = solverDesc["convCheck"];
                }
                nlohmann::json convCheckDesc;
                convCheckDesc["type"] = convCheckType;
                SmartPtr<StdConvCheck<typename TAlgebra::vector_type>> convCheck;
                convCheck = CreateConvCheck<TAlgebra>(convCheckDesc);
                newtonSolver->set_convergence_check(convCheck);

                // line search
                std::string lineSearchType = "standard";
                if (json_default_nonlinearSolver["newton"]["lineSearch"].is_string()){
                    lineSearchType = json_default_nonlinearSolver["newton"]["lineSearch"];
                }
                // get descriptor for line search
                if (solverDesc.contains("lineSearch")){
                    lineSearchType = solverDesc["lineSearch"];
                }
                nlohmann::json lineSearchDesc;
                lineSearchDesc["type"] = lineSearchType;
                // if line_Search not null
                if (!lineSearchType.empty() && lineSearchType != "none"){
                    SmartPtr<StandardLineSearch<typename TAlgebra::vector_type>> lineSearch;
                    lineSearch = CreateLineSearch<TAlgebra>(lineSearchDesc);
                    newtonSolver->set_line_search(lineSearch);
                }
                // set reassemble_J_freq if exists
                if (solverDesc.contains("reassemble_J_freq") && solverDesc["reassemble_J_freq"].is_number()){
                    newtonSolver->set_reassemble_J_freq(solverDesc["reassemble_J_freq"]);
                }

                // TODO: SetDebugWriter(newtonSolver, solverDesc, defaults, solverutil);

                return newtonSolver;
            }
            UG_THROW("CreateNewtonSolver: Only 'newton' solver type is supported in this function.")
        }

        template <typename TDomain, typename TAlgebra>
        SmartPtr<ILinearIterator<typename TAlgebra::vector_type>>
        CreatePreconditioner(nlohmann::json &desc, SolverUtil<TDomain, TAlgebra> &solverutil){
            typedef typename TAlgebra::vector_type TVector;
            typedef std::vector<size_t> ordering_container_type;
            typedef ILinearIterator<typename TAlgebra::vector_type> TPrecond;

            // TODO: Implement preconditioner creation behavior based on 'desc' and 'solverutil'
            // TODO: Check if preset

            nlohmann::json json_default_preconds = json_predefined_defaults::solvers.at("preconditioner");

            nlohmann::json json_default_transfer = json_predefined_defaults::solvers["transfer"]["std"];

            nlohmann::json json_default_mgStats = json_predefined_defaults::solvers["mgStats"]["standard"];

            SmartPtr<TPrecond> preconditioner;

            SmartPtr<ApproximationSpace<TDomain>> approxSpace = NullSmartPtr();

            if (solverutil.hasComponent("approxSpace")){
                // approxSpace = solverutil.getComponent("approxSpace");
            }

            std::string type = desc["type"];
            if (type == "ilu"){
                UG_LOG("CreatePreconditioner ilu \n")
                // create ilu
                typedef ILU<TAlgebra> TILU;
                SmartPtr<TILU> ILU = make_sp(new TILU());

                // configure ilu
                number beta = json_default_preconds["ilu"]["beta"];
                if (desc.contains("beta")){
                    beta = desc["beta"];
                }
                ILU->set_beta(beta);

                number damping = json_default_preconds["ilu"]["damping"];
                if (desc.contains("damping")){
                    damping = desc["damping"];
                }
                ILU->set_damp(damping);

                bool sort = json_default_preconds["ilu"]["sort"];
                if (desc.contains("sort")){
                    sort = desc["sort"];
                }
                ILU->set_sort(sort);

                number sortEps = json_default_preconds["ilu"]["sortEps"];
                if (desc.contains("sortEps")){
                    sortEps = desc["sortEps"];
                }
                ILU->set_sort_eps(sortEps);

                number inversionEps = json_default_preconds["ilu"]["inversionEps"];
                if (desc.contains("inversionEps")){
                    inversionEps = desc["inversionEps"];
                }
                ILU->set_inversion_eps(inversionEps);

                bool consistentInterfaces = json_default_preconds["ilu"]["consistentInterfaces"];
                if (desc.contains("consistentInterfaces")){
                    consistentInterfaces = desc["consistentInterfaces"];
                }
                ILU->enable_consistent_interfaces(consistentInterfaces);

                bool overlap = json_default_preconds["ilu"]["overlap"];
                if (desc.contains("overlap")){
                    overlap = desc["overlap"];
                }
                ILU->enable_overlap(overlap);

                if (desc.contains("ordering"))
                {
                    SmartPtr<IOrderingAlgorithm<TAlgebra, ordering_container_type>> ordering;
                    ordering = CreateOrdering<TDomain, TAlgebra>(desc);

                    ILU->set_ordering_algorithm(ordering);
                }
                // set ordering (no default value)
                // TODO: ordering = CreateOrdering
                // TODO: precond.set_ordering_algorithm(ordering)

                // Cast ILU to IPreconditioner for return
                preconditioner = ILU.template cast_static<TPrecond>();
            }
            else if (type == "ilut"){
                UG_LOG("CreatePreconditioner, creating ilut \n")
                // create ilut
                typedef ILUTPreconditioner<TAlgebra> TILUT;

                number threshold = json_default_preconds["ilut"]["threshold"];
                if (desc.contains("threshold")){
                    threshold = desc["threshold"];
                }
                SmartPtr<TILUT> ILUT = make_sp(new TILUT(threshold));
                std::string name = json_default_preconds["ilut"]["ordering"];

                if (desc.contains("ordering"))
                {

                    SmartPtr<IOrderingAlgorithm<TAlgebra, ordering_container_type>> ordering;
                    ordering = CreateOrdering<TDomain, TAlgebra>(desc);
                    ILUT->set_ordering_algorithm(ordering);
                }
                else
                {
                    typedef NativeCuthillMcKeeOrdering<TAlgebra, ordering_container_type> TNCMK;
                    SmartPtr<TNCMK> ordering = make_sp(new TNCMK());
                    ILUT->set_ordering_algorithm(ordering);
                }
            }
            // TODO: ordering = CreateOrdering
            // TODO: precond.set_ordering_algorithm(ordering)

            else if (type == "jac")
            {
                UG_LOG("CreatePreconditioner jacobi \n")
                typedef Jacobi<TAlgebra> TJAC;
                SmartPtr<TJAC> JAC = make_sp(new TJAC());
                UG_LOG("jacobi, damping default\n")
                number damping = json_default_preconds["jac"]["damping"];
                UG_LOG("jacobi, damping desc\n")
                if (desc.contains("damping")){
                    damping = desc["damping"];
                }
                UG_LOG("jacobi, set damping\n")
                JAC->set_damp(damping);
                preconditioner = JAC.template cast_static<TPrecond>();
            }
            else if (type == "gs"){
                UG_LOG("CreatePreconditioner gauss seidel\n")
                typedef GaussSeidel<TAlgebra> TGS;
                SmartPtr<TGS> GS = make_sp(new TGS());
                UG_LOG("consistentInterfaces default\n")

                bool consistentInterfaces = json_default_preconds["gs"]["consistentInterfaces"];
                UG_LOG("consistentInterfaces desc\n")
                if (desc.contains("consistentInterfaces")){
                    consistentInterfaces = desc["consistentInterfaces"];
                }
                UG_LOG("enable consistentInterfaces\n")
                GS->enable_consistent_interfaces(consistentInterfaces);

                bool overlap = json_default_preconds["gs"]["overlap"];
                if (desc.contains("overlap")){
                    overlap = desc["overlap"];
                }
                UG_LOG("enable overlap\n")
                GS->enable_overlap(overlap);
                preconditioner = GS.template cast_static<TPrecond>();
            }
            else if (type == "sgs"){
                UG_LOG("CreatePreconditioner symmetric gauss seidel\n")
                typedef SymmetricGaussSeidel<TAlgebra> TSGS;
                SmartPtr<TSGS> SGS = make_sp(new TSGS());
                bool consistentInterfaces = json_default_preconds["sgs"]["consistentInterfaces"];
                if (desc["sgs"].contains("consistentInterfaces")){
                    consistentInterfaces = desc["sgs"]["consistentInterfaces"];
                }
                SGS->enable_consistent_interfaces(consistentInterfaces);

                bool overlap = json_default_preconds["sgs"]["overlap"];
                if (desc["sgs"].contains("overlap")){
                    overlap = desc["sgs"]["overlap"];
                }
                SGS->enable_overlap(overlap);
                preconditioner = SGS.template cast_static<TPrecond>();
            }
            else if (type == "egs"){
                UG_LOG("CreatePreconditioner element gauss seidel\n")
                typedef ElementGaussSeidel<TDomain, TAlgebra> TEGS;
                SmartPtr<TEGS> EGS = make_sp(new TEGS());
                preconditioner = EGS.template cast_static<TPrecond>();
            }
            else if (type == "cgs"){
                UG_LOG("CreatePreconditioner component gauss seidel\n")
                std::vector<std::string> vFullRowCmp;
                typedef ComponentGaussSeidel<TDomain, TAlgebra> TCGS;

                number relax = 1.0;
                if (json_default_preconds.contains("cgs") && json_default_preconds["cgs"].contains("relax")){
                    relax = json_default_preconds["cgs"]["relax"];
                }
                if (desc.contains("relax")){
                    relax = desc["relax"];
                }

                // UG_LOG("json_default_preconds: " << json_default_preconds.dump(4) << "\n");
                // UG_LOG("desc: " << desc.dump(4) << "\n");

                SmartPtr<TCGS> CGS = make_sp(new TCGS(relax, vFullRowCmp));

                number alpha = 1.0;
                if (json_default_preconds.contains("cgs") && json_default_preconds["cgs"].contains("alpha")){
                    alpha = json_default_preconds["cgs"]["alpha"];
                }
                if (desc.contains("alpha")){
                    alpha = desc["alpha"];
                }
                UG_LOG("cgs, set alpha\n")
                CGS->set_alpha(alpha);

                number beta = 1.0;
                if (json_default_preconds.contains("cgs") && json_default_preconds["cgs"].contains("beta")){
                    beta = json_default_preconds["cgs"]["beta"];
                }
                if (desc.contains("beta")){
                    beta = desc["beta"];
                }
                UG_LOG("cgs, set beta\n")
                CGS->set_beta(beta);

                bool weights = false;
                if (json_default_preconds.contains("cgs") && json_default_preconds["cgs"].contains("weights")){
                    weights = json_default_preconds["cgs"]["weights"];
                }
                if (desc.contains("weights")){
                    weights = desc["weights"];
                }
                UG_LOG("cgs, enable weights\n")
                CGS->set_weights(weights);

                preconditioner = CGS.template cast_static<TPrecond>();
            }
            else if (type == "ssc"){
                // TODO:Duy create ssc
            }
            else if (type == "gmg"){
                UG_LOG("CreatePreconditioner Geometric MultiGrid (GMG)\n");
                // idea: dont convert the objects of the lua table into json give them directly to the CreatePreconditioner Function
                typedef AssembledMultiGridCycle<TDomain, TAlgebra> TGMG;
                SmartPtr<ApproximationSpace<TDomain>> approxSpace;
                SmartPtr<TGMG> GMG = make_sp(new TGMG());
                if (solverutil.hasComponent("approxSpace")){
                    UG_LOG("ApproximationSpace found!\n");
                    approxSpace = std::get<SmartPtr<ApproximationSpace<TDomain>>>(
                        solverutil.getComponent("approxSpace"));

                    SmartPtr<TGMG> GMG = make_sp(new TGMG(approxSpace));
                }
                else{
                    UG_LOG("An ApproximationSpace is required to create a 'gmg' solver.\n");
                    exit(0);
                }
                std::string baseSolverType = json_default_preconds["gmg"]["baseSolver"];
                if (desc["gmg"].contains("baseSolver")){
                    baseSolverType = desc["gmg"]["baseSolver"];
                }
                nlohmann::json baseSolverDesc;
                baseSolverDesc["type"] = baseSolverType;
                // SmartPtr<LinearSolver<typename TAlgebra::vector_type>> baseSolver = CreateLinearSolver(baseSolverDesc, solverutil);
                UG_LOG("Solver found!\n");
                GMG->set_base_solver(CreateLinearSolver(baseSolverDesc, solverutil));

                std::string smootherType = json_default_preconds["gmg"]["smoother"];
                if (desc["gmg"].contains("smoother")){
                    smootherType = desc["gmg"]["smoother"];
                }
                nlohmann::json smootherTypeDesc;
                smootherTypeDesc["type"] = smootherType;
                SmartPtr<ILinearIterator<TVector>> smoother = CreatePreconditioner(smootherTypeDesc, solverutil);
                UG_LOG("Smoother found!\n");
                GMG->set_smoother(smoother);

                number baseLevel = json_default_preconds["gmg"]["baseLevel"];
                if (desc["gmg"].contains("baseLevel")){
                    baseLevel = desc["gmg"]["baseLevel"];
                }
                UG_LOG("base_level found!\n");
                GMG->set_base_level(baseLevel);

                std::string cycleType = json_default_preconds["gmg"]["cycle"];
                if (desc["gmg"].contains("cycle")){
                    cycleType = desc["gmg"]["cycle"];
                }
                UG_LOG("CycleType found!\n");
                GMG->set_cycle_type(cycleType);

                // Not so sure if it should contain the actual parameter
                if (solverutil.hasComponent("discretization")){
                    SmartPtr<IAssemble<TAlgebra>> discretization = std::get<SmartPtr<IAssemble<TAlgebra>>>(
                        solverutil.getComponent("discretization"));
                    UG_LOG("Discretization found!\n");
                    GMG->set_discretization(discretization);
                }
                number preSmooth = json_default_preconds["gmg"]["preSmooth"];
                if (desc["gmg"].contains("preSmooth")){
                    preSmooth = desc["gmg"]["preSmooth"];
                }
                UG_LOG("presmooth found!\n");
                GMG->set_num_presmooth(preSmooth);

                number postSmooth = json_default_preconds["gmg"]["postSmooth"];
                if (desc["gmg"].contains("postSmooth")){
                    postSmooth = desc["gmg"]["postSmooth"];
                }
                UG_LOG("postsmooth found!\n");
                GMG->set_num_postsmooth(postSmooth);

                bool rap = json_default_preconds["gmg"]["rap"];
                if (desc["gmg"].contains("rap")){
                    rap = desc["gmg"]["rap"];
                }
                UG_LOG("rap found!\n");
                GMG->set_rap(rap);

                bool rim = json_default_preconds["gmg"]["rim"];
                if (desc["gmg"].contains("rim")){
                    rim = desc["gmg"]["rim"];
                }
                UG_LOG("smooth_on_surface_rim found!\n");
                GMG->set_smooth_on_surface_rim(rim);

                bool emulateFullRefined = json_default_preconds["gmg"]["emulateFullRefined"];
                if (desc["gmg"].contains("emulateFullRefined")){
                    emulateFullRefined = desc["gmg"]["emulateFullRefined"];
                }
                UG_LOG("emulate_full_refined_grid found!\n");
                GMG->set_emulate_full_refined_grid(emulateFullRefined);

                UG_LOG("transfer started!\n")

                typedef StdTransfer<TDomain, TAlgebra> TFD;
                SmartPtr<TFD> TF = make_sp(new TFD());
                UG_LOG("SmartPtr created\n");
                std::cout << json_default_transfer.dump() << std::endl;
                number restrictionDamp = json_default_transfer["restrictionDamp"];
                UG_LOG("Defaultrestiction!\n");
                if (desc["gmg"].contains("transfer") && desc["gmg"]["transfer"].contains("restrictionDamp")){
                    restrictionDamp = desc["gmg"]["transfer"]["restrictionDamp"];
                    UG_LOG("restiction!\n");
                }
                UG_LOG("restictionDamp!");
                TF->set_restriction_damping(restrictionDamp);

                number prolongationDamp = json_default_transfer["prolongationDamp"];
                if (desc["gmg"].contains("transfer") && desc["gmg"]["transfer"].contains("prolongationDamp")){
                    prolongationDamp = desc["gmg"]["transfer"]["prolongationDamp"];
                }
                UG_LOG("prolongation!");
                TF->set_prolongation_damping(prolongationDamp);

                bool optimization = json_default_transfer["enableP1LagrangeOptimization"];
                if (desc["gmg"].contains("transfer") && desc["gmg"]["transfer"].
                    contains("enableP1LagrangeOptimization")){
                    optimization = desc["gmg"]["transfer"]["enableP1LagrangeOptimization"];
                }
                UG_LOG("lagrange_opt!");
                TF->enable_p1_lagrange_optimization(optimization);

                // TF->set_debug(TF, transferDesc, transferDefault, solverutil);
                UG_LOG("transfer found!\n");
                GMG->set_transfer(TF);

                // entweder mit boolean Objekt bauen oder direkt Objekt abfragen

                SmartPtr<GridFunctionDebugWriter<TDomain, TAlgebra>> debugDesc;
                if (desc["gmg"].contains("debug")){
                    bool debug = desc["gmg"]["debug"];
                    if (debug){
                        bool vtk = true;
                        bool conn_viewer = false;
                        typedef GridFunctionDebugWriter<TDomain, TAlgebra> TDW;
                        SmartPtr<TDW> DWW = make_sp(new TDW(approxSpace));
                        DWW->set_conn_viewer_output(conn_viewer);
                        DWW->set_vtk_output(vtk);
                        GMG->set_debug(DWW);
                    }
                }

                if (solverutil.hasComponent("Debugger") and !desc["gmg"].contains("debug")){
                    debugDesc = SetDebugger(desc, solverutil);
                    GMG->set_debug(debugDesc);
                }

                bool gatheredBaseSolverIfAmbiguous = json_default_preconds["gmg"]["gatheredBaseSolverIfAmbiguous"];
                UG_LOG("gatheredBaseSolverIfAMb found!\n")
                if (desc["gmg"].contains("gatheredBaseSolverIfAmbiguous")){
                    gatheredBaseSolverIfAmbiguous = desc["gmg"]["gatheredBaseSolverIfAmbiguous"];
                }
                GMG->set_gathered_base_solver_if_ambiguous(gatheredBaseSolverIfAmbiguous);
                UG_LOG("beginn mgStats!\n")
                UG_LOG("mgStats found!\n")
                std::cout << json_default_mgStats.dump() << std::endl;
                typedef MGStats<TDomain, TAlgebra> MGS;
                SmartPtr<MGS> MGSD = make_sp(new MGS());

                std::string prefix = json_default_mgStats["filenamePrefix"];
                if (desc.contains("mgStats") && desc.contains("standard") && desc.contains("filenamePrefix")){
                    prefix = desc["gmg"]["mgStats"]["standard"]["filenamePrefix"];
                }
                MGSD->set_filename_prefix(prefix.c_str());

                bool exitError = json_default_mgStats["exitOnError"];
                if (desc.contains("mgStats") && desc.contains("standard") && desc.contains("exitOnError")){
                    exitError = desc["gmg"]["mgStats"]["standard"]["exitOnError"];
                }
                MGSD->set_exit_on_error(exitError);

                bool errorVec = json_default_mgStats["writeErrVecs"];
                if (desc.contains("mgStats") && desc.contains("standard") && desc.contains("writeErrVecs")){
                    errorVec = desc["gmg"]["mgStats"]["standard"]["writeErrVecs"];
                }
                MGSD->set_write_err_vecs(errorVec);

                bool errorDiff = json_default_mgStats["writeErrDiffs"];
                if (desc.contains("mgStats") && desc.contains("standard") && desc.contains("writeErrDiffs")){
                    errorDiff = desc["gmg"]["mgStats"]["standard"]["writeErrDiffs"];
                }
                MGSD->set_write_err_diffs(errorDiff);
                if (!json_default_mgStats["activeStages"].is_null()){
                    std::vector<int> activeStage = json_default_mgStats["activeStages"].get<std::vector<int>>();
                    if (desc.contains("mgStats") && desc.contains("standard") && desc.contains("activeStages")){
                        activeStage = desc["gmg"]["mgStats"]["standard"]["activeStages"].get<std::vector<int>>();
                    }
                    MGSD->set_active_stages(activeStage);
                }
                GMG->set_mg_stats(MGSD);

                preconditioner = GMG.template cast_static<TPrecond>();
            }
            else if (type == "schur"){
                // TODO:Duy create schur
            }

            // return Preconditioner
            return preconditioner;
        }

        template<typename TDomain, typename TAlgebra>
        void PrepareSolverUtil(nlohmann::json &desc, nlohmann::json &solverutil){
            typedef SolverUtil<TDomain, TAlgebra> TSolverUtil;
            // Create SolverUtil container class
            SmartPtr<TSolverUtil> solv_util = make_sp(new TSolverUtil());

            if (solverutil.contains("ApproxSpace")){
                solv_util->setComponent("ApproxSpace", solverutil["ApproxSpace"]);
            }
        }
        typedef std::vector<size_t> ordering_container_type;
        template <typename TDomain, typename TAlgebra>
        SmartPtr<IOrderingAlgorithm<TAlgebra, ordering_container_type>>
        CreateOrdering(nlohmann::json &desc)
        {
            typedef std::vector<size_t> ordering_container_type;

            /*  std::variant<
                  SmartPtr<NativeCuthillMcKeeOrdering<TAlgebra, ordering_container_type>>,
                  SmartPtr<BoostCuthillMcKeeOrdering<TAlgebra, ordering_container_type>>,
                  SmartPtr<BoostMinimumDegreeOrdering<TAlgebra, ordering_container_type>>,
                  SmartPtr<TopologicalOrdering<TAlgebra, ordering_container_type>>,
                  SmartPtr<SCCOrdering<TAlgebra, ordering_container_type>>>
                  */

            typedef IOrderingAlgorithm<TAlgebra, ordering_container_type> TBase;
            SmartPtr<TBase> ordering;
            if (desc["ordering"].is_string() && !(desc["ordering"].is_null()))
            {
                std::string name = desc["ordering"];
                if (name == "NativeCuthillMcKee")
                {
                    typedef NativeCuthillMcKeeOrdering<TAlgebra, ordering_container_type> TNCMK;
                    SmartPtr<TNCMK> order = make_sp(new TNCMK());
                    ordering = order.template cast_dynamic<TBase>();
                }
                if (name == "ReverseNativeCuthillMcKee")
                {
                    typedef NativeCuthillMcKeeOrdering<TAlgebra, ordering_container_type> TNCMK;
                    SmartPtr<TNCMK> order = make_sp(new TNCMK());
                    order->set_reverse(true);
                    ordering = order.template cast_dynamic<TBase>();
                }
                if (name == "BoostCuthillMcKee")
                {
                    typedef BoostCuthillMcKeeOrdering<TAlgebra, ordering_container_type> TBCMK;
                    SmartPtr<TBCMK> order = make_sp(new TBCMK());
                    ordering = order.template cast_dynamic<TBase>();
                }
                if (name == "ReverseBoostCuthillMcKee")
                {
                    typedef BoostCuthillMcKeeOrdering<TAlgebra, ordering_container_type> TBCMK;
                    SmartPtr<TBCMK> order = make_sp(new TBCMK());
                    order->set_reverse(true);
                    ordering = order.template cast_dynamic<TBase>();
                }
                if (name == "BoostMinimumDegree")
                {
                    typedef BoostMinimumDegreeOrdering<TAlgebra, ordering_container_type> TBMD;
                    SmartPtr<TBMD> order = make_sp(new TBMD());
                    ordering = order.template cast_dynamic<TBase>();
                }
                if (name == "Topological")
                {
                    typedef TopologicalOrdering<TAlgebra, ordering_container_type> TT;
                    SmartPtr<TT> order = make_sp(new TT());
                    ordering = order.template cast_dynamic<TBase>();
                }
            }
            if (desc["ordering"].is_object())
            {
                nlohmann::json orderingTypeDesc;
                if (desc["ordering"]["name"] == "NativeCuthillMcKee")
                {
                    typedef NativeCuthillMcKeeOrdering<TAlgebra, ordering_container_type> TNCMK;
                    SmartPtr<TNCMK> order = make_sp(new TNCMK());
                    ordering = order.template cast_dynamic<TBase>();
                }
                if (desc["ordering"]["name"] == "ReverseNativeCuthillMcKee")
                {
                    typedef NativeCuthillMcKeeOrdering<TAlgebra, ordering_container_type> TNCMK;
                    SmartPtr<TNCMK> order = make_sp(new TNCMK());
                    order->set_reverse(true);
                    ordering = order.template cast_dynamic<TBase>();
                }
                if (desc["ordering"]["name"] == "BoostCuthillMcKee")
                {
                    typedef BoostCuthillMcKeeOrdering<TAlgebra, ordering_container_type> TBCMK;
                    SmartPtr<TBCMK> order = make_sp(new TBCMK());
                    ordering = order.template cast_dynamic<TBase>();
                }
                if (desc["ordering"]["name"] == "ReverseBoostCuthillMcKee")
                {
                    typedef BoostCuthillMcKeeOrdering<TAlgebra, ordering_container_type> TBCMK;
                    SmartPtr<TBCMK> order = make_sp(new TBCMK());
                    order->set_reverse(true);
                    ordering = order.template cast_dynamic<TBase>();
                }
                if (desc["ordering"]["name"] == "BoostMinimumDegree")
                {
                    typedef BoostMinimumDegreeOrdering<TAlgebra, ordering_container_type> TBMD;
                    SmartPtr<TBMD> order = make_sp(new TBMD());
                    ordering = order.template cast_dynamic<TBase>();
                }
                if (desc["ordering"]["name"] == "Topological")
                {
                    typedef TopologicalOrdering<TAlgebra, ordering_container_type> TT;
                    SmartPtr<TT> order = make_sp(new TT());
                    ordering = order.template cast_dynamic<TBase>();
                }
                if (desc["ordering"]["name"] == "SSC")
                {
                    typedef SCCOrdering<TAlgebra, ordering_container_type> TSSC;
                    SmartPtr<TSSC> order = make_sp(new TSSC());
                    std::string sub = desc["ordering"]["subalgo"];
                    nlohmann::json subDesc;
                    subDesc["name"] = sub;
                    using TBase = IOrderingAlgorithm<TAlgebra, ordering_container_type>;
                    SmartPtr<TBase> subAlgo = CreateOrdering<TDomain, TAlgebra>(subDesc);
                    order->set_ordering_subalgorithm(subAlgo);
                    ordering = order.template cast_dynamic<TBase>();
                }
                if (desc["ordering"]["name"] == "Lex")
                {
                    typedef LexOrdering<TAlgebra, TDomain, ordering_container_type> TL;
                    SmartPtr<TL> order = make_sp(new TL());
                    std::string dir = desc["ordering"]["dir"];
                    order->set_direction(dir.c_str());
                    ordering = order.template cast_dynamic<TBase>();
                }
                if (desc["ordering"]["name"] == "River")
                {
                    typedef RiverOrdering<TAlgebra, TDomain, ordering_container_type> TR;
                    SmartPtr<TR> order = make_sp(new TR());
                    std::string sources = desc["ordering"]["sources"];
                    order->select_sources(sources.c_str());
                    ordering = order.template cast_dynamic<TBase>();
                }
            }

            /* using TBase = IOrderingAlgorithm<TAlgebra, ordering_container_type>;
             SmartPtr<TBase> order = std::visit([](auto &algo) -> SmartPtr<TBase>
                                                { return SmartPtr<TBase>(algo.get()); }, ordering);
             */
            return ordering;
        }

        /*
        Set Debugger for different preconditioners and solvers

        parameters nlohmann::json &desc
        parameters SolverUtil<TDomain, TAlgebra> &solverutil
        return GridFunctionDebugWriter<TDomain, TAlgebra>
        */
        template<typename TDomain, typename TAlgebra>
        SmartPtr<GridFunctionDebugWriter<TDomain, TAlgebra>> SetDebugger(
            nlohmann::json &desc, SolverUtil<TDomain, TAlgebra> &solverutil){
            SmartPtr<GridFunctionDebugWriter<TDomain, TAlgebra>> debugger;
            debugger = std::get<SmartPtr<GridFunctionDebugWriter<TDomain, TAlgebra>>>(
                solverutil.getComponent("debugger"));
            return debugger;
        }
        // Needs tests independent from other util.lua
        template <typename TDomain, typename TAlgebra>
        std::tuple<MatrixOperator<typename TAlgebra::matrix_type, typename TAlgebra::vector_type>,
                   GridFunction<ug::GridFunction<TDomain, TAlgebra>, typename TAlgebra::vector_type>,
                   GridFunction<ug::GridFunction<TDomain, TAlgebra>, typename TAlgebra::vector_type>>
        SolveLinearProblem(SmartPtr<DomainDiscretization<TDomain, TAlgebra>> &domainDisc,
                           nlohmann::json &solverDesc, SolverUtil<TDomain, TAlgebra> &solverutil,
                           std::string outFilePrefix,
                           SmartPtr<UserData<number, TDomain::dim, number>> startValueCB = nullptr)
        {

            typedef ug::GridFunction<TDomain, TAlgebra> TFct;
            typedef typename TAlgebra::vector_type vector_type;
            typedef MatrixOperator<typename TAlgebra::matrix_type, typename TAlgebra::vector_type> MO;
            SmartPtr<MO> A = make_sp(new MO());
            typedef GridFunction<TFct, vector_type> GF;
            SmartPtr<GF> u = make_sp(new GF());
            SmartPtr<GF> b = make_sp(new GF());

            A = AssembledLinearOperator(domainDisc);
            u = GridFunction(domainDisc->approximation_space());
            b = GridFunction(domainDisc->approximation_space());
            if (!startValueCB.valid()){
                u->set(0.0);
            }
            else if (startValueCB->is_const()) {
                u->set(startValueCB->const_value());
            }
            domainDisc->adjust_solution(u);
            domainDisc->assemble_linear(A, b);

            std::string type;

            if (solverDesc.contains("type")) {
                type = solverDesc["type"];
                if (type == "newton") {
                    auto solver = CreateNewtonSolver(solverDesc, solverutil);
                    solver->init(A, u);
                    solver->apply(u, b);
                }
                else {
                    auto solver = CreateLinearSolver(solverDesc, solverutil);
                    solver->init(A, u);
                    solver->apply(u, b);
                }
            }

            if (outFilePrefix.empty())
            {
                WriteGridFunctionToVTK(u, outFilePrefix);
                    SaveVectorForConnectionViewer(u, (outFilePrefix + ".vec").c_str());
            }

            return A, u, b;
        }
        // Need tests dependent from util.SolveLinearTimeProblem
        template<typename TDomain, typename TAlgebra>
        nlohmann::json SolveLinearTimeProblem(SmartPtr<DomainDiscretization<TDomain, TAlgebra>> &domainDisc,
                nlohmann::json &solverDesc, nlohmann::json &timeDesc, nlohmann::json &outDesc,
                SolverUtil<TDomain, TAlgebra> &solverutil,
                SmartPtr<UserData<number, TDomain::dim, number>> startValueCB = nullptr){
            static const int dim = TDomain::dim;
            nlohmann::json json_default_solver;
            nlohmann::json json_solver;
        if (timeDesc.is_object()) {
            json_solver = timeDesc;
            std::string name;
            if (json_solver.contains("name") || json_solver.contains("type")) {
                name = json_solver["name"];
                name = json_solver["type"];
            }
            if (json_predefined_defaults::solvers.contains(name)) {
                json_default_solver = json_predefined_defaults::solvers[name];
            }
        }

            number start = json_default_solver["start"];
            if (json_solver.contains("start")) {
                start = json_solver["start"];
            }
            number stop = json_default_solver["stop"];
            if (json_solver.contains("stop")) {
                stop = json_solver["stop"];
            }
            std::string scheme = json_default_solver["scheme"];
            if (json_solver.contains("scheme")) {
                scheme = json_solver["scheme"];
            }
            number dt = json_default_solver["dt"];
            if (json_solver.contains("dt")) {
                dt = json_solver["dt"];
            }



        typedef ug::GridFunction<TDomain, TAlgebra> TFct;
        typedef typename TAlgebra::vector_type vector_type;

        typedef GridFunction<TFct, vector_type> GF;
        SmartPtr<GF> u = make_sp(new GF());
        u = GridFunction(domainDisc->approximation_space());
        if (!startValueCB.valid()){
                u->set(0.0);
        }
        else if (startValueCB->is_const()) {
                u->set(startValueCB->const_value());
        }


            std::string type;
            auto solver = CreateSolver(solverDesc, solverutil);
            if (solverDesc.contains("type")) {
                type = solverDesc["type"];
                if (type == "newton") {
                    auto solver = CreateNewtonSolver(solverDesc, solverutil);
                }
                else {
                    auto solver = CreateLinearSolver(solverDesc, solverutil);
                }
            }



        typedef VTKOutput<dim> VTKO;
            SmartPtr<VTKO> vtkOut = make_sp(new VTKO());

        std::string outFilePrefix;
        if (outDesc.is_string()) {
            outFilePrefix = outDesc.get<std::string>();
        }
        if (outDesc.is_object()) {
            if(outDesc.contains("prefix")){
            outFilePrefix += outDesc["prefix"];
            }
            if (outDesc.contains("vtkOut")) {
                vtkOut = outDesc["vtkOut"];
            }
        }
            nlohmann::json SolverLinearTimeProblemDesc;
            SolverLinearTimeProblemDesc["u"] = u;
            SolverLinearTimeProblemDesc["DomainDesc"] = domainDisc;
            SolverLinearTimeProblemDesc["Solver"] = solver;
            SolverLinearTimeProblemDesc["vtkOut"] = vtkOut;
            SolverLinearTimeProblemDesc["scheme"] = scheme;
            SolverLinearTimeProblemDesc["number"] = 1;
            SolverLinearTimeProblemDesc["start"] = start;
            SolverLinearTimeProblemDesc["stop"] = stop;
            SolverLinearTimeProblemDesc["dt"] = dt;
            return SolverLinearTimeProblemDesc;
        }

        // Need tests dependent from util.SolveNonLinearTimeProblem
        template<typename TDomain, typename TAlgebra>
        nlohmann::json SolveNonLinearTimeProblem(SmartPtr<DomainDiscretization<TDomain, TAlgebra>> &domainDisc,
                nlohmann::json &solverDesc, nlohmann::json &timeDesc, nlohmann::json &outDesc,
                SolverUtil<TDomain, TAlgebra> &solverutil,
                SmartPtr<UserData<number, TDomain::dim, number>> startValueCB = nullptr){
            static const int dim = TDomain::dim;
            nlohmann::json json_default_solver;
            nlohmann::json json_solver;
            if (timeDesc.is_object()) {
                json_solver = timeDesc;
                std::string name;
                if (json_solver.contains("name") || json_solver.contains("type")) {
                    name = json_solver["name"];
                    name = json_solver["type"];
                }
                if (json_predefined_defaults::solvers.contains(name)) {
                    json_default_solver = json_predefined_defaults::solvers[name];
                }
            }

            number start = json_default_solver["start"];
            if (json_solver.contains("start")) {
                start = json_solver["start"];
            }
            number stop = json_default_solver["stop"];
            if (json_solver.contains("stop")) {
                stop = json_solver["stop"];
            }
            std::string scheme = json_default_solver["scheme"];
            if (json_solver.contains("scheme")) {
                scheme = json_solver["scheme"];
            }
            number dt = json_default_solver["dt"];
            if (json_solver.contains("dt")) {
                dt = json_solver["dt"];
            }
            typedef ug::GridFunction<TDomain, TAlgebra> TFct;
            typedef typename TAlgebra::vector_type vector_type;

            typedef GridFunction<TFct, vector_type> GF;
            SmartPtr<GF> u = make_sp(new GF());
            u = GridFunction(domainDisc->approximation_space());
            if (!startValueCB.valid()){
                u->set(0.0);
            }
            else if (startValueCB->is_const()) {
                u->set(startValueCB->const_value());
            }
            std::string type;
            auto solver = CreateSolver(solverDesc, solverutil);
            if (solverDesc.contains("type")) {
                type = solverDesc["type"];
                if (type == "newton") {
                    auto solver = CreateNewtonSolver(solverDesc, solverutil);
                }
                else {
                    auto solver = CreateLinearSolver(solverDesc, solverutil);
                }
            }





            typedef VTKOutput<dim> VTKO;
            SmartPtr<VTKO> vtkOut = make_sp(new VTKO());
            std::string outFilePrefix;
            if (outDesc.is_string()) {
                outFilePrefix = outDesc.get<std::string>();
            }
            if (outDesc.is_object()) {
                if(outDesc.contains("prefix")){
                    outFilePrefix += outDesc["prefix"];
                }
                if (outDesc.contains("vtkOut")) {
                    vtkOut = outDesc["vtkOut"];
                }
            }
            nlohmann::json SolverLinearTimeProblemDesc;
            SolverLinearTimeProblemDesc["u"] = u;
            SolverLinearTimeProblemDesc["DomainDesc"] = domainDisc;
            SolverLinearTimeProblemDesc["Solver"] = solver;
            SolverLinearTimeProblemDesc["vtkOut"] = vtkOut;
            SolverLinearTimeProblemDesc["scheme"] = scheme;
            SolverLinearTimeProblemDesc["number"] = 1;
            SolverLinearTimeProblemDesc["start"] = start;
            SolverLinearTimeProblemDesc["stop"] = stop;
            SolverLinearTimeProblemDesc["dt"] = dt;
            return SolverLinearTimeProblemDesc;
        }
        /*
         * Helper class to provide c++ util functions in lua.
         * We can instaciate this object in lua via
         * local functionProvider = SolverUtilFunctionProvider()
         * and automatically get the correct templated class,
         * e.g SolverUtilFunctionProvider2dCPU1
         * this means functionprovider automatically chooses
         * the correct templated util functions!
         */

        template<typename TDomain, typename TAlgebra>
        class SolverUtilFunctionProvider{
        public:
            typedef std::vector<size_t> ordering_container_type;
            typedef typename TAlgebra::vector_type vector_type;
            typedef typename TAlgebra::matrix_type matrix_type;
            const static int dim = TDomain::dim;

            SolverUtilFunctionProvider(){
            };

            SmartPtr<ILinearIterator<typename TAlgebra::vector_type>>
            GetCreatePreconditioner(nlohmann::json &desc, SolverUtil<TDomain, TAlgebra> &solverutil){
                return CreatePreconditioner<TDomain, TAlgebra>(desc, solverutil);
            }

            SmartPtr<ILinearIterator<typename TAlgebra::vector_type>>
            GetCreateLinearSolver(nlohmann::json &desc, SolverUtil<TDomain, TAlgebra> &solverutil){
                return CreateLinearSolver<TDomain, TAlgebra>(desc, solverutil);
            }

            SmartPtr<NewtonSolver<TAlgebra>>
            GetCreateNewtonSolver(nlohmann::json &desc, SolverUtil<TDomain, TAlgebra> &solverutil){
                return CreateNewtonSolver<TDomain, TAlgebra>(desc, solverutil);
            }

            SmartPtr<StandardLineSearch<vector_type>>
            GetCreateLineSearch(nlohmann::json &desc){
                return CreateLineSearch<TAlgebra>(desc);
            }

            SmartPtr<StdConvCheck<vector_type>>
            GetCreateConvCheck(nlohmann::json &desc){
                return CreateConvCheck<TAlgebra>(desc);
            }
            SmartPtr<GridFunctionDebugWriter<TDomain, TAlgebra>> GetSetDebugger(
                nlohmann::json &desc, SolverUtil<TDomain, TAlgebra> &solverutil)
            {
                return SetDebugger<TDomain, TAlgebra>(desc, solverutil);
            }
            SmartPtr<IOrderingAlgorithm<TAlgebra, ordering_container_type>>
            GetCreateOrdering(nlohmann::json &desc)
            {
                return CreateOrdering<TAlgebra, ordering_container_type>(desc);
            }
            std::tuple<MatrixOperator<typename TAlgebra::matrix_type, typename TAlgebra::vector_type>,
                   GridFunction<ug::GridFunction<TDomain, TAlgebra>, typename TAlgebra::vector_type>,
                   GridFunction<ug::GridFunction<TDomain, TAlgebra>, typename TAlgebra::vector_type>>

            GetSolveLinearProblem(SmartPtr<DomainDiscretization<TDomain, TAlgebra>> &domainDisc,
                           nlohmann::json &solverDesc, SolverUtil<TDomain, TAlgebra> &solverutil,
                           std::string outFilePrefix,
                           SmartPtr<UserData<number, TDomain::dim, number>> startValueCB = nullptr) {
                return SolveLinearProblem<TDomain,TAlgebra>(solverDesc, solverutil, startValueCB, outFilePrefix);
            }
            nlohmann::json
            GetSolveLinearTimeProblem(SmartPtr<DomainDiscretization<TDomain, TAlgebra>> &domainDisc,
                nlohmann::json &solverDesc, nlohmann::json &timeDesc, nlohmann::json &outDesc,
                SolverUtil<TDomain, TAlgebra> &solverutil,
                SmartPtr<UserData<number, TDomain::dim, number>> startValueCB = nullptr) {
                return SolveLinearTimeProblem<TDomain,TAlgebra>(domainDisc, solverDesc, timeDesc, outDesc, solverutil,startValueCB);
            }
            nlohmann::json
            GetSolveNonLinearTimeProblem(SmartPtr<DomainDiscretization<TDomain, TAlgebra>> &domainDisc,
                nlohmann::json &solverDesc, nlohmann::json &timeDesc, nlohmann::json &outDesc,
                SolverUtil<TDomain, TAlgebra> &solverutil,
                SmartPtr<UserData<number, TDomain::dim, number>> startValueCB = nullptr) {
                return SolveNonLinearTimeProblem<TDomain,TAlgebra>(domainDisc, solverDesc, timeDesc, outDesc, solverutil,startValueCB);
            }
        };
    } // namespace util
} // namespace ug
#endif // UG4_SOLVER_UTIL_H
