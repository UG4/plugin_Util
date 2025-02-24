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
#include "lib_disc/operator/linear_operator/multi_grid_solver/mg_solver.h"
#include "lib_disc/operator/non_linear_operator/newton_solver/newton.h"
namespace ug
{
    namespace Util
    {

        void CondAbort(bool condition, std::string message)
        {
            UG_ASSERT(!condition, "ERROR in util.solver: " << message);
        }
        template <typename TAlgebra>
        SmartPtr<StdConvCheck<typename TAlgebra::vector_type>> CreateConvCheck(nlohmann::json &descriptor)
        {

            typedef typename TAlgebra::vector_type vector_type;
            // set parameters from descriptor attributes
            bool verbose = false;
            if (descriptor.contains("verbose"))
            {
                verbose = descriptor["verbose"];
            }
            std::string name = "standard";
            if (descriptor.contains("type"))
            {
                name = descriptor["type"];
            }
            number iterations = 100;
            if (descriptor.contains("iterations"))
            {
                iterations = descriptor["iterations"];
            }
            number reduction = 1e-6;
            if (descriptor.contains("reduction"))
            {
                reduction = descriptor["reduction"];
            }
            number absolute = 1e-12;
            if (descriptor.contains("absolute"))
            {
                absolute = descriptor["absolute"];
            }
            bool suppress_unsuccessful = false;
            if (descriptor.contains("suppress_unsuccessful"))
            {
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
        template <typename TDomain, typename TAlgebra>
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
            SmartPtr<NewtonSolver<TAlgebra>>,
            SmartPtr<AssembledMultiGridCycle<TDomain, TAlgebra>>
            // missing ElementGaussSeidel
            >;

        template <typename TDomain, typename TAlgebra>
        class SolverUtil
        {
            /// Container class for objects used in the CreateSolver routine.
            /*
             * components are stored as std::variant and populate an unordered
             * map for access. Calls to various subroutines are made
             * */
        public:
            void setComponent(const std::string &key, SolverComponent<TDomain, TAlgebra> component)
            {
                components[key] = component;
            }

            SolverComponent<TDomain, TAlgebra> getComponent(const std::string &key) const
            {
                auto it = components.find(key);
                if (it != components.end())
                {
                    return it->second;
                }
                UG_THROW("SolverUtil.getComponent(" << key << "): no component named " << key << " found");
            }

            bool hasComponent(const std::string &key) const
            {

                auto it = components.find(key);
                if (it != components.end())
                {
                    return true;
                }
                return false;
            }

            template <typename T>
            SmartPtr<T> getComponentAs(const std::string &key) const
            {
                return std::get<SmartPtr<T>>(getComponent(key));
            }

        private:
            std::unordered_map<std::string, SolverComponent<TDomain, TAlgebra>> components;
        };

        template <typename TDomain, typename TAlgebra>
        SmartPtr<SolverUtil<TDomain, TAlgebra>> CreateSolver(nlohmann::json &solverDesc,
                                                             SolverUtil<TDomain, TAlgebra> solverutil)
        {

            // PrepareSolverUtil<TDomain, TAlgebra>();
            std::string type = "linear";
            if (solverDesc.contains("type"))
            {
                type = solverDesc["type"];
            }
            if (type == "newton")
            {
                auto newton_solver = make_sp(new NewtonSolver<TAlgebra>());
                // call createlinearsolver

                // get descriptor for linear solver
                if (solverDesc.contains("linSolver"))
                {
                    // CreateLinearSolver<TDomain, TAlgebra>(solverDesc, solverutil);
                }

                // line search
                if (solverDesc.contains("lineSearch"))
                {
                    // CreateLineSearch(solverDesc["lineSearch"], solverutil);
                }
                newton_solver->set_convergence_check(CreateConvCheck<TAlgebra>(solverDesc["convCheck"]), solverutil);
            }
        }

        template <typename TDomain, typename TAlgebra>
        SmartPtr<LinearSolver<typename TAlgebra::vector_type>> CreateLinearSolver(nlohmann::json &desc,
                                                                                  SolverUtil<TDomain, TAlgebra> solverutil)
        {
            typedef typename TAlgebra::vector_type TVector;
            typedef LinearSolver<TVector> TLinSolv;

            // TODO: is preset

            bool create_precond = false;
            bool create_conv_check = false;

            // if no descriptor given, create default linear solver
            if (!desc.contains("LinearSolver"))
            {

                SmartPtr<ILU<TAlgebra>> ilu = make_sp<ILU<TAlgebra>>(new ILU<TAlgebra>());
                SmartPtr<StdConvCheck<TVector>> convCheck = make_sp<StdConvCheck<TVector>>(
                    new StdConvCheck<TVector>(100, 1e-9, 1e-12));
                SmartPtr<TLinSolv> default_linear_solver = make_sp(new TLinSolv());
                default_linear_solver->set_convergence_check(convCheck);
                default_linear_solver->set_preconditioner(ilu);
                return default_linear_solver;
            }
        }

        template <typename TDomain, typename TAlgebra>
        SmartPtr<ILinearIterator<typename TAlgebra::vector_type>>
        CreatePreconditioner(nlohmann::json &desc, SolverUtil<TDomain, TAlgebra> &solverutil)
        {
            typedef typename TAlgebra::vector_type TVector;
            typedef ILinearIterator<TVector> TPrecond;

            // TODO: Implement preconditioner creation behavior based on 'desc' and 'solverutil'
            // TODO: Check if preset

            nlohmann::json json_default_preconds = json_predefined_defaults::solvers.at("preconditioner");

            SmartPtr<TPrecond> preconditioner;

            SmartPtr<ApproximationSpace<TDomain>> approxSpace = NullSmartPtr();

            if (solverutil.hasComponent("approxSpace"))
            {
                // approxSpace = solverutil.getComponent("approxSpace");
            }

            std::string type = desc["type"];
            if (type == "ilu")
            {
                // create ilu
                typedef ILU<TAlgebra> TILU;
                SmartPtr<TILU> ILU = make_sp(new TILU());

                // configure ilu
                number beta = json_default_preconds["ilu"]["beta"];
                if (desc.contains("beta"))
                {
                    beta = desc["beta"];
                }
                ILU->set_beta(beta);

                number damping = json_default_preconds["ilu"]["damping"];
                if (desc.contains("damping"))
                {
                    damping = desc["damping"];
                }
                ILU->set_damp(damping);

                bool sort = json_default_preconds["ilu"]["sort"];
                if (desc.contains("sort"))
                {
                    sort = desc["sort"];
                }
                ILU->set_sort(sort);

                number sortEps = json_default_preconds["ilu"]["sortEps"];
                if (desc.contains("sortEps"))
                {
                    sortEps = desc["sortEps"];
                }
                ILU->set_sort_eps(sortEps);

                number inversionEps = json_default_preconds["ilu"]["inversionEps"];
                if (desc.contains("inversionEps"))
                {
                    inversionEps = desc["inversionEps"];
                }
                ILU->set_inversion_eps(inversionEps);

                bool consistentInterfaces = json_default_preconds["ilu"]["consistentInterfaces"];
                if (desc.contains("consistentInterfaces"))
                {
                    consistentInterfaces = desc["consistentInterfaces"];
                }
                ILU->enable_consistent_interfaces(consistentInterfaces);

                bool overlap = json_default_preconds["ilu"]["overlap"];
                if (desc.contains("overlap"))
                {
                    overlap = desc["overlap"];
                }
                ILU->enable_overlap(overlap);
                // set ordering (no default value)
                // TODO: ordering = CreateOrdering
                // TODO: precond.set_ordering_algorithm(ordering)

                // Cast ILU to IPreconditioner for return
                preconditioner = ILU.template cast_static<TPrecond>();
            }
            else if (type == "ilut")
            {
                // create ilut
                typedef ILUTPreconditioner<TAlgebra> TILUT;

                number threshold = json_default_preconds["ilut"]["threshold"];
                if (desc["ilut"].contains("threshold"))
                {
                    threshold = desc["ilut"]["threshold"];
                }
                SmartPtr<TILUT> ILUT = make_sp(new TILUT(threshold));

                // TODO: ordering = CreateOrdering
                // TODO: precond.set_ordering_algorithm(ordering)
            }
            else if (type == "jac")
            {
                UG_LOG("creating jacobi\n")
                typedef Jacobi<TAlgebra> TJAC;
                SmartPtr<TJAC> JAC = make_sp(new TJAC());
                UG_LOG("damping default\n")
                number damping = json_default_preconds["jac"]["damping"];
                UG_LOG("damping desc\n")
                if (desc["jac"].contains("damping"))
                {
                    damping = desc["jac"]["damping"];
                }
                UG_LOG("set damping\n")
                JAC->set_damp(damping);
                preconditioner = JAC.template cast_static<TPrecond>();
            }
            else if (type == "gs")
            {
                UG_LOG("creating gauss seidel\n")
                typedef GaussSeidel<TAlgebra> TGS;
                SmartPtr<TGS> GS = make_sp(new TGS());
                UG_LOG("consistentInterfaces default\n")
                bool consistentInterfaces = json_default_preconds["gs"]["consistentInterfaces"];
                UG_LOG("consistentInterfaces desc\n")
                if (desc["gs"].contains("consistentInterfaces"))
                {
                    consistentInterfaces = desc["gs"]["consistentInterfaces"];
                }
                UG_LOG("enable consistentInterfaces\n")
                GS->enable_consistent_interfaces(consistentInterfaces);

                bool overlap = json_default_preconds["gs"]["overlap"];
                if (desc["gs"].contains("overlap"))
                {
                    overlap = desc["gs"]["overlap"];
                }
                UG_LOG("enable overlap\n")
                GS->enable_overlap(overlap);
                preconditioner = GS.template cast_static<TPrecond>();
            }
            else if (type == "sgs")
            {
                UG_LOG("creating symmetric gauss seidel\n")
                typedef SymmetricGaussSeidel<TAlgebra> TSGS;
                SmartPtr<TSGS> SGS = make_sp(new TSGS());
                bool consistentInterfaces = json_default_preconds["sgs"]["consistentInterfaces"];
                if (desc["sgs"].contains("consistentInterfaces"))
                {
                    consistentInterfaces = desc["sgs"]["consistentInterfaces"];
                }
                SGS->enable_consistent_interfaces(consistentInterfaces);

                bool overlap = json_default_preconds["sgs"]["overlap"];
                if (desc["sgs"].contains("overlap"))
                {
                    overlap = desc["sgs"]["overlap"];
                }
                SGS->enable_overlap(overlap);
                preconditioner = SGS.template cast_static<TPrecond>();
            }
            else if (type == "egs")
            {
                UG_LOG("creating element gauss seidel\n")
                typedef ElementGaussSeidel<TDomain, TAlgebra> TEGS;
                SmartPtr<TEGS> EGS = make_sp(new TEGS());
                preconditioner = EGS.template cast_static<TPrecond>();
            }
            else if (type == "cgs")
            { // Component Gauss Seidel in progress
              // UG_LOG("creating component gauss seidel\n")
              // typedef ComponentGaussSeidel<TDomain, TAlgebra> TCGS;
              // SmartPtr<TCGS> CGS = make_sp(new TCGS());
              // preconditioner = CGS.template cast_static<TPrecond>();
            }
            else if (type == "ssc")
            {
                // TODO:Duy create ssc
            }
            else if (type == "gmg")
            {

                UG_LOG("Creating Geometric MultiGrid (GMG)\n");

                typedef AssembledMultiGridCycle<TDomain, TAlgebra> TGMG;
                SmartPtr<ApproximationSpace<TDomain>> approxSpace;
                SmartPtr<TGMG> GMG = make_sp(new TGMG());
                if (!desc["gmg"]["approxSpace"].is_null())
                {
                    SmartPtr<TGMG> GMG = make_sp(new TGMG(approxSpace));
                }
                else
                {
                    // desc["gmg"]["approxSpace"].template cast_dynamic<ApproximationSpace<TDomain>>->print_statistic();
                }
                std::string baseSolverType = json_default_preconds["gmg"]["baseSolver"];
                if (desc["gmg"].contains("baseSolver"))
                {
                    baseSolverType = desc["gmg"]["baseSolver"];
                }
                nlohmann::json baseSolverDesc;
                baseSolverDesc["type"] = baseSolverType;
                SmartPtr<LinearSolver<typename TAlgebra::vector_type>> baseSolver = CreateLinearSolver(baseSolverDesc, solverutil);
                GMG->set_base_solver(baseSolver);

                std::string smootherType = json_default_preconds["gmg"]["smoother"];
                if (desc["gmg"].contains("smoother"))
                {
                    smootherType = desc["gmg"]["smoother"];
                }
                nlohmann::json smootherTypeDesc;
                smootherTypeDesc["type"] = smootherType;
                SmartPtr<ILinearIterator<TVector>> smoother = CreatePreconditioner(smootherTypeDesc, solverutil);
                GMG->set_smoother(smoother);

                number baseLevel = json_default_preconds["gmg"]["baseLevel"];
                if (desc["gmg"].contains("baseLevel"))
                {
                    baseLevel = desc["gmg"]["baseLevel"];
                }
                GMG->set_base_level(baseLevel);

                std::string cycleType = json_default_preconds["gmg"]["cycle"];
                if (desc["gmg"].contains("cycle"))
                {
                    cycleType = desc["gmg"]["cycle"];
                }
                GMG->set_cycle_type(cycleType);
                // Not so sure if it should contain the actual parameter
                number discretization = json_default_preconds["gmg"]["discretization"];
                SmartPtr<ug::IAssemble<TAlgebra>> discretizationDesc;
                if (desc["gmg"].contains("discretization"))
                {
                    number discretization = desc["gmg"]["discretization"];
                    GMG->set_discretization(discretizationDesc);
                }
                else
                {
                }
                number preSmooth = json_default_preconds["gmg"]["preSmooth"];
                if (desc["gmg"].contains("preSmooth"))
                {
                    preSmooth = desc["gmg"]["preSmooth"];
                }
                GMG->set_num_presmooth(preSmooth);

                number postSmooth = json_default_preconds["gmg"]["postSmooth"];
                if (desc["gmg"].contains("postSmooth"))
                {
                    postSmooth = desc["gmg"]["postSmooth"];
                }
                GMG->set_num_postsmooth(postSmooth);

                bool rap = json_default_preconds["gmg"]["rap"];
                if (desc["gmg"].contains("rap"))
                {
                    rap = desc["gmg"]["rap"];
                }
                GMG->set_rap(rap);

                bool rim = json_default_preconds["gmg"]["rim"];
                if (desc["gmg"].contains("rim"))
                {
                    rim = desc["gmg"]["rim"];
                }
                GMG->set_smooth_on_surface_rim(rim);

                bool emulateFullRefined = json_default_preconds["gmg"]["emulateFullRefined"];
                if (desc["gmg"].contains("emulateFullRefined"))
                {
                    emulateFullRefined = desc["gmg"]["emulateFullRefined"];
                }
                GMG->set_emulate_full_refined_grid(emulateFullRefined);

                std::string transfer = json_default_preconds["gmg"]["transfer"];
                if (desc["gmg"].contains("transfer"))
                {
                    std::string transfer = desc["gmg"]["transfer"];
                }
                if (transfer == "std")
                {
                    nlohmann::json transferDesc = desc["gmg"]["transfer"]["std"];
                    nlohmann::json transferDefault = json_default_preconds["transfer"]["std"];
                    typedef StdTransfer<TDomain, TAlgebra> TFD;
                    SmartPtr<TFD> TF = make_sp(new TFD());

                    number restrictionDamp = transferDefault["restrictionDamp"];
                    if (desc["gmg"].contains("restrictionDamp"))
                    {
                        number restrictionDamp = transferDesc["restrictionDamp"];
                    }
                    TF->set_restriction_damping(restrictionDamp);

                    number prolongationDamp = transferDefault["prolongationDamp"];
                    if (desc["gmg"].contains("prolongationDamp"))
                    {
                        number prolongationDamp = transferDesc["prolongationDamp"];
                    }
                    TF->set_prolongation_damping(prolongationDamp);

                    bool optimization = transferDefault["enableP1LagrangeOptimization"];
                    if (desc["gmg"].contains("enableP1LagrangeOptimization"))
                    {
                        bool optimization = transferDesc["enableP1LagrangeOptimization"];
                    }
                    TF->enable_p1_lagrange_optimization(optimization);

                    // TF->set_debug(TF, transferDesc, transferDefault, solverutil);
                    GMG->set_transfer(TF);
                }

                bool debug = json_default_preconds["gmg"]["debug"];
                SmartPtr<ug::GridFunctionDebugWriter<TDomain, TAlgebra>> debugDesc;
                if (desc["gmg"].contains("debug"))
                {
                    bool debug = desc["gmg"]["debug"];
                }
                if (debug == true)
                {
                    GMG->set_debug(debugDesc);
                }
                else
                {
                }

                bool gatheredBaseSolverIfAmbiguous = json_default_preconds["gmg"]["gatheredBaseSolverIfAmbiguous"];
                if (desc["gmg"].contains("gatheredBaseSolverIfAmbiguous"))
                {
                    gatheredBaseSolverIfAmbiguous = desc["gmg"]["gatheredBaseSolverIfAmbiguous"];
                }
                GMG->set_gathered_base_solver_if_ambiguous(gatheredBaseSolverIfAmbiguous);

                std::string mgStats = json_default_preconds["gmg"]["mgStats"];
                if (!desc["gmg"]["mgStats"].is_null())
                {

                    if (mgStats == "standard")
                    {
                        nlohmann::json mgStatsDesc = desc["gmg"]["mgStats"]["standard"];
                        nlohmann::json mgStatsDefault = json_default_preconds["mgStats"]["standard"];
                        typedef MGStats<TDomain, TAlgebra> MGS;
                        SmartPtr<MGS> MGSD = make_sp(new MGS());

                        std::string prefix = mgStatsDefault["filenamePrefix"];
                        if (mgStatsDesc.contains("filenamePrefix"))
                        {
                            std::string prefix = mgStatsDesc["filenamePrefix"];
                        }
                        MGSD->set_filename_prefix(prefix.c_str());

                        bool exitError = mgStatsDefault["exitonError"];
                        if (mgStatsDesc.contains("exitonError"))
                        {
                            bool exitError = mgStatsDesc["exitonError"];
                        }
                        MGSD->set_exit_on_error(exitError);

                        bool errorVec = mgStatsDefault["writeErrVecs"];
                        if (mgStatsDesc.contains("writeErrVecs"))
                        {
                            bool errorVec = mgStatsDefault["writeErrVecs"];
                        }
                        MGSD->set_write_err_vecs(errorVec);

                        bool errorDiff = mgStatsDefault["writeErrDiffs"];
                        if (mgStatsDesc.contains("writeErrDiffs"))
                        {
                            bool errorDiff = mgStatsDesc["writeErrDiffs"];
                        }
                        MGSD->set_write_err_diffs(errorDiff);

                        std::vector<int> activeStage = mgStatsDefault["activeStages"];
                        if (mgStatsDesc.contains("activeStages"))
                        {
                            std::vector<int> activeStage = mgStatsDesc["activeStages"];
                        }
                        MGSD->set_active_stages(activeStage);

                        GMG->set_mg_stats(MGSD);
                    }
                    else
                    {
                    }
                }

                preconditioner = GMG.template cast_static<TPrecond>();
            }
            else if (type == "schur")
            {
                // TODO:Duy create schur
            }

            // return Preconditioner
            return preconditioner;
        }

        template <typename TAlgebra>
        SmartPtr<StandardLineSearch<typename TAlgebra::vector_type>> CreateLineSearch(nlohmann::json &desc)
        {
            // typedef for convenience
            typedef StandardLineSearch<typename TAlgebra::vector_type> line_search_type;
            SmartPtr<line_search_type> ls;

            // load defaults
            nlohmann::json json_default_lineSearch = json_predefined_defaults::solvers["lineSearch"];

            // handle type of line search
            // default
            std::string type = "standard";
            // input type
            if (desc.contains("type"))
            {
                type = desc["type"];
            }

            if (type == "standard")
            {
                // handle parameters of standard line search
                int maxSteps = json_default_lineSearch[type]["maxSteps"];
                if (desc.contains("maxSteps"))
                {
                    maxSteps = desc["maxSteps"];
                }
                number lambdaStart = json_default_lineSearch[type]["lambdaStart"];
                if (desc.contains("lambdaStart"))
                {
                    lambdaStart = desc["lambdaStart"];
                }
                number lambdaReduce = json_default_lineSearch[type]["lambdaReduce"];
                if (desc.contains("lambdaReduce"))
                {
                    lambdaReduce = desc["lambdaReduce"];
                }
                bool acceptBest = json_default_lineSearch[type]["acceptBest"];
                if (desc.contains("acceptBest"))
                {
                    acceptBest = desc["acceptBest"];
                }
                bool checkAll = json_default_lineSearch[type]["checkAll"];
                if (desc.contains("checkAll"))
                {
                    checkAll = desc["checkAll"];
                }
                // create line search with chosen parameters
                ls = make_sp(new line_search_type(maxSteps,
                                                  lambdaStart,
                                                  lambdaReduce,
                                                  acceptBest,
                                                  checkAll));
                // set conditional parameters
                if (desc.contains("verbose"))
                {
                    bool verbose = desc["verbose"];
                    ls->set_verbose(verbose);
                }
                if (desc.contains("suffDesc"))
                {
                    number suffDesc = desc["suffDesc"];
                    ls->set_suff_descent_factor(suffDesc);
                }
                if (desc.contains("maxDefect"))
                {
                    number maxDefect = desc["maxDefect"];
                    ls->set_maximum_defect(maxDefect);
                }
            }
            // force exit if line search is invalid
            CondAbort(ls.invalid(), "Invalid line-search specified: " + type);

            return ls;
        }

        template <typename TDomain, typename TAlgebra>
        void PrepateSolverUtil(nlohmann::json &desc, nlohmann::json &solverutil)
        {

            typedef SolverUtil<TDomain, TAlgebra> TSolverUtil;
            // Create SolverUtil container class
            SmartPtr<TSolverUtil> solv_util = make_sp(new TSolverUtil());

            if (solverutil.contains("ApproxSpace"))
            {
                solv_util->setComponent("ApproxSpace", solverutil["ApproxSpace"]);
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
        template <typename TDomain, typename TAlgebra>
        class SolverUtilFunctionProvider
        {
        public:
            typedef typename TAlgebra::vector_type vector_type;
            typedef typename TAlgebra::matrix_type matrix_type;
            const static int dim = TDomain::dim;

            SolverUtilFunctionProvider() {};

            SmartPtr<ILinearIterator<vector_type>> GetCreatePreconditioner(nlohmann::json &desc, SolverUtil<TDomain, TAlgebra> &solverutil)
            {
                return CreatePreconditioner<TDomain, TAlgebra>(desc, solverutil);
            }

            SmartPtr<StandardLineSearch<vector_type>> GetCreateLineSearch(nlohmann::json &desc)
            {
                return CreateLineSearch<TAlgebra>(desc);
            }
        };
    } // namespace util
} // namespace ug
#endif // UG4_SOLVER_UTIL_H
