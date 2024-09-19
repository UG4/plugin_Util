//
// Created by julian on 6/4/24.
//

#ifndef UG4_SOLVER_UTIL_H
#define UG4_SOLVER_UTIL_H

#include <nlohmann/json.hpp>

#include "common/common.h"
#include "common/log.h"
#include "common/util/smart_pointer.h"
#include <memory>
#include <unordered_map>
#include <variant>
#include <stdexcept>
#include "lib_algebra/operator/convergence_check.h"
#include "lib_algebra/operator/linear_solver/linear_solver.h"
// include solver components
#include "lib_disc/function_spaces/grid_function.h"
#include "lib_algebra/operator/preconditioner/ilu.h"
#include "lib_algebra/operator/preconditioner/ilut.h"
#include "lib_algebra/operator/preconditioner/jacobi.h"
#include "lib_algebra/operator/preconditioner/block_gauss_seidel.h"
#include "lib_algebra/operator/preconditioner/gauss_seidel.h"
#include "lib_disc/operator/linear_operator/element_gauss_seidel/element_gauss_seidel.h"
#include "lib_disc/operator/non_linear_operator/newton_solver/newton.h"
namespace ug{
namespace Util {

    void CondAbort(bool condition, std::string message){
        UG_ASSERT(!condition, "ERROR in util.solver: " << message);
    }
    template<typename TAlgebra>
    SmartPtr<StdConvCheck <typename TAlgebra::vector_type>> CreateConvCheck(nlohmann::json& descriptor) {

        typedef typename TAlgebra::vector_type vector_type;
        // set parameters from descriptor attributes
        bool verbose = false;
        if(descriptor.contains("verbose")){
            verbose = descriptor["verbose"];
        }
        std::string name = "standard";
        if(descriptor.contains("type")){
            name = descriptor["type"];
        }
        number iterations = 100;
        if(descriptor.contains("iterations")){
            iterations = descriptor["iterations"];
        }
        number reduction = 1e-6;
        if(descriptor.contains("reduction")){
            reduction = descriptor["reduction"];
        }
        number absolute = 1e-12;
        if(descriptor.contains("absolute")){
            absolute = descriptor["absolute"];
        }
        bool suppress_unsuccessful = false;
        if(descriptor.contains("suppress_unsuccessful")){
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
                SmartPtr<NewtonSolver<TAlgebra>>
        // missing ElementGaussSeidel
        >;

        template<typename TDomain, typename TAlgebra>
        class SolverUtil{
        /// Container class for objects used in the CreateSolver routine.
        /*
         * components are stored as std::variant and populate an unordered
         * map for access. Calls to various subroutines are made
         * */
        public:

            void setComponent(const std::string& key, SolverComponent<TDomain, TAlgebra> component){
                components[key] = component;
            }

            SolverComponent<TDomain, TAlgebra> getComponent(const std::string& key) const {
                auto it = components.find(key);
                if (it != components.end()){
                    return it->second;
                }
                UG_THROW("SolverUtil.getComponent(" << key << "): no component named " << key << " found");
            }

            bool hasComponent(const std::string& key) const{

                auto it = components.find(key);
                if (it != components.end()){
                    return true;
                }
                return false;

            }

            template<typename T>
            SmartPtr<T> getComponentAs(const std::string& key) const {
                return std::get<SmartPtr<T>>(getComponent(key));
            }

        private:
            std::unordered_map<std::string, SolverComponent<TDomain, TAlgebra>> components;
        };




        template<typename TDomain, typename TAlgebra>
        SmartPtr<SolverUtil<TDomain, TAlgebra>> CreateSolver(nlohmann::json& solverDesc,
                                                             SolverUtil<TDomain, TAlgebra> solverutil){

            //PrepareSolverUtil<TDomain, TAlgebra>();
            std::string type = "linear";
            if(solverDesc.contains("type")){
                type = solverDesc["type"];
            }
            if(type == "newton"){
                auto newton_solver= make_sp(new NewtonSolver<TAlgebra>());
               // call createlinearsolver

               // get descriptor for linear solver
                if(solverDesc.contains("linSolver")){
                   //CreateLinearSolver<TDomain, TAlgebra>(solverDesc, solverutil);
                }

               // line search
                if(solverDesc.contains("lineSearch")){
                    //CreateLineSearch(solverDesc["lineSearch"], solverutil);
                }
                newton_solver->set_convergence_check(CreateConvCheck<TAlgebra>(solverDesc["convCheck"]), solverutil);
            }

        }





        template<typename TDomain, typename TAlgebra>
        SmartPtr<LinearSolver<typename TAlgebra::vector_type>> CreateLinearSolver(nlohmann::json& desc,
                                                                                  SolverUtil<TDomain, TAlgebra> solverutil){
            typedef typename TAlgebra::vector_type TVector;
            typedef LinearSolver<TVector> TLinSolv;

            // TODO: is preset

            bool create_precond = false;
            bool create_conv_check = false;

            // if no descriptor given, create default linear solver
            if(!desc.contains("LinearSolver")){

                SmartPtr<ILU<TAlgebra>> ilu = make_sp<ILU<TAlgebra>>(new ILU<TAlgebra>());
                SmartPtr<StdConvCheck<TVector>> convCheck = make_sp<StdConvCheck<TVector>>(
                new StdConvCheck<TVector>(100, 1e-9, 1e-12));
                SmartPtr<TLinSolv> default_linear_solver = make_sp(new TLinSolv());
                default_linear_solver->set_convergence_check(convCheck);
                default_linear_solver->set_preconditioner(ilu);
                return default_linear_solver;
                
            }



        }

    template<typename TDomain, typename TAlgebra>
    SmartPtr <IPreconditioner<TAlgebra>>
    CreatePreconditioner(nlohmann::json &desc, SolverUtil<TDomain, TAlgebra> &solverutil) {
        typedef typename TAlgebra::vector_type TVector;
        typedef IPreconditioner<TAlgebra> TPrecond;

        // TODO: Implement preconditioner creation behavior based on 'desc' and 'solverutil'
        // TODO: Check if preset

        nlohmann::json json_default_preconds = json_predefined_defaults::solvers.at("preconditioner");

        SmartPtr <TPrecond> preconditioner;

        SmartPtr<ApproximationSpace<TDomain>> approxSpace = NullSmartPtr();

        if(solverutil.hasComponent("approxSpace")){
            //approxSpace = solverutil.getComponent("approxSpace");
        }

        std::string type = desc["type"];
        if(type == "ilu"){
            // create ilu
            typedef ILU<TAlgebra> TILU;
            SmartPtr<TILU> ILU = make_sp(new TILU());

            // configure ilu
            number beta = json_default_preconds["ilu"]["beta"];
            if(desc.contains("beta")){
                beta = desc["beta"];
            }
            ILU->set_beta(beta);

            number damping = json_default_preconds["ilu"]["damping"];
            if(desc.contains("damping")){
                damping = desc["damping"];
            }
            ILU->set_damp(damping);

            bool sort = json_default_preconds["ilu"]["sort"];
            if(desc.contains("sort")){
                sort = desc["sort"];
            }
            ILU->set_sort(sort);

            number sortEps = json_default_preconds["ilu"]["sortEps"];
            if(desc.contains("sortEps")){
                sortEps = desc["sortEps"];
            }
            ILU->set_sort_eps(sortEps);

            number inversionEps = json_default_preconds["ilu"]["inversionEps"];
            if(desc.contains("inversionEps")){
                inversionEps = desc["inversionEps"];
            }
            ILU->set_inversion_eps(inversionEps);

            bool consistentInterfaces = json_default_preconds["ilu"]["consistentInterfaces"];
            if(desc.contains("consistentInterfaces")){
                consistentInterfaces = desc["consistentInterfaces"];
            }
            ILU->enable_consistent_interfaces(consistentInterfaces);

            bool overlap = json_default_preconds["ilu"]["overlap"];
            if(desc.contains("overlap")){
                overlap = desc["overlap"];
            }
            ILU->enable_overlap(overlap);
            // set ordering (no default value)
            // TODO: ordering = CreateOrdering
            // TODO: precond.set_ordering_algorithm(ordering)

            // Cast ILU to IPreconditioner for return
            preconditioner = ILU.template cast_static<TPrecond>();

        }
        else if(desc.contains("ilut")){
            // create ilut
            typedef ILUTPreconditioner<TAlgebra> TILUT;

            number threshold = json_default_preconds["ilut"]["threshold"];
            if(desc["ilut"].contains("threshold")){
                threshold = desc["ilut"]["threshold"];
            }
            SmartPtr<TILUT> ILUT = make_sp(new TILUT(threshold));

            // TODO: ordering = CreateOrdering
            // TODO: precond.set_ordering_algorithm(ordering)


}
        else if(desc.contains("jac")){
            //TODO:Duy createjac
        }
        else if(desc.contains("gs")){
            
            UG_LOG("creating gauss seidel\n")
            typedef GaussSeidel<TAlgebra> TGS;
            SmartPtr<TGS> GS = make_sp(new TGS());
            UG_LOG("consistentInterfaces default\n")
            bool consistentInterfaces = json_default_preconds["gs"]["consistentInterfaces"];
            UG_LOG("consistentInterfaces desc\n")
            if(desc["gs"].contains("consistentInterfaces")){
                consistentInterfaces = desc["gs"]["consistentInterfaces"];
            }
            UG_LOG("enable consistentInterfaces\n")
            GS->enable_consistent_interfaces(consistentInterfaces);


            bool overlap = json_default_preconds["gs"]["overlap"];
            if(desc["gs"].contains("overlap")){
                overlap = desc["gs"]["overlap"];
            }
            UG_LOG("enable overlap\n")
            GS->enable_overlap(overlap);
            preconditioner = GS.template cast_static<TPrecond>();
        }
        else if(desc.contains("sgs")){
            //TODO:Tim create sgs
        }
        else if(desc.contains("egs")){
            // create egs
        }
        else if(desc.contains("cgs")){
            // create cgs
        }
        else if(desc.contains("ssc")){
            //TODO:Duy create ssc
        }
        else if(desc.contains("gmg")){
            // create gmg
        }
        else if(desc.contains("schur")){
            //TODO:Duy create schur
        }

// return Preconditioner
        return preconditioner;
    }

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

        return ls;
    }

        template<typename TDomain, typename TAlgebra>
        void PrepateSolverUtil(nlohmann::json& desc, nlohmann::json& solverutil){

            typedef SolverUtil<TDomain, TAlgebra> TSolverUtil;
            // Create SolverUtil container class
            SmartPtr<TSolverUtil> solv_util = make_sp(new TSolverUtil());

            if(solverutil.contains("ApproxSpace")){
                solv_util->setComponent("ApproxSpace",solverutil["ApproxSpace"]);
            }
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
    typedef typename TAlgebra::vector_type              vector_type;
    typedef typename TAlgebra::matrix_type              matrix_type;
    const static int dim = TDomain::dim;

    SolverUtilFunctionProvider(){};

    SmartPtr <IPreconditioner<TAlgebra>> GetCreatePreconditioner(nlohmann::json &desc, SolverUtil<TDomain, TAlgebra> &solverutil){
        return CreatePreconditioner<TDomain, TAlgebra>(desc, solverutil);
    }

    SmartPtr<StandardLineSearch<vector_type>> GetCreateLineSearch(nlohmann::json &desc){
        return CreateLineSearch<TAlgebra>(desc);
    }


};
} //namespace util
} //namespace ug
#endif //UG4_SOLVER_UTIL_H
