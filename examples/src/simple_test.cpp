#include "simple_test.h"

int main()
{
    std::shared_ptr<BasicStates> si = std::make_shared<BasicStates>(2, 1);

    // Declare model variables
    casadi::SX dt = casadi::SX::sym("dt", 1);
    casadi::SX dx = casadi::SX::sym("dx", 2, 1);

    casadi::SX x = casadi::SX::sym("x", 2, 1);
    casadi::SX u = casadi::SX::sym("u", 1);

    // Model equations
    casadi::SX xdot = vertcat((1 - pow(x(1), 2)) * x(0) - x(1) + u, x(0));

    // Objective term
    casadi::SX l = pow(x(0), 2) + pow(x(1), 2) + pow(u, 2);

    casadi::SX x1 = casadi::SX::sym("x1", 2, 1);
    casadi::SX x2 = casadi::SX::sym("x2", 2, 1);

    // Continuous time dynamics (all state variables are Euclidean so we do not need to worry about these manifold operators)
    casadi::Function Fint("Fint", {x, dx, dt}, {dx});
    casadi::Function Fdiff("Fdiff", {x1, x2, dt}, {x2});
    casadi::Function F_state_error("F_state_error", {x1, x2}, {x1 - x2});
    casadi::Function F("F", {x, u}, {xdot});
    casadi::Function L("L", {x, u}, {l});
    casadi::Function Phi("Phi", {x}, {0});

    std::shared_ptr<GeneralProblemData> gp_data = std::make_shared<GeneralProblemData>(Fint, Fdiff, F_state_error, Phi);
    casadi::DM X0 = casadi::DM::zeros(si->nx, 1);
    X0(1, 0) = 1;

    int d = 1;
    int N = 200;
    double T = 10.;
    double h = T / N;

    std::shared_ptr<BasicSequence> seq = std::make_shared<BasicSequence>();
    BasicMode mode;
    seq->addPhase(mode, N, T, F);
    seq->FillPhaseCost(0, L);

    std::shared_ptr<RosenbrockProblemData> problem = std::make_shared<RosenbrockProblemData>(gp_data, seq, si, x, u, dt);
    // std::shared_ptr<ConstraintBuilder<RosenbrockProblemData>> simple_builder =
    //     std::make_shared<SimpleConstraintBuilder<RosenbrockProblemData>>();

    std::vector<std::shared_ptr<ConstraintBuilder<RosenbrockProblemData>>> builders = {};
    std::shared_ptr<DecisionDataBuilder<RosenbrockProblemData>> decision_builder = std::make_shared<SimpleDecisionDataBuilder<RosenbrockProblemData>>();

    // casadi::Dict opts;
    // opts["ipopt.linear_solver"] = "ma97";
    // opts["ipopt.ma97_order"] = "metis";
    // opts["ipopt.fixed_variable_treatment"] = "make_constraint";
    // opts["ipopt.max_iter"] = 250;
    // std::shared_ptr<TrajectoryOpt<RosenbrockProblemData, BasicMode>> traj = std::make_shared<TrajectoryOpt<RosenbrockProblemData, BasicMode>>(problem, seq, builders, decision_builder, false, opts, "ipopt");

    casadi::Dict opts;
    opts["start"] = "hot";
    opts["snopt.Cold Start/Warm Start"] = "Warm";
    opts["snopt.Major iterations limit"] = 250;
    std::shared_ptr<TrajectoryOpt<RosenbrockProblemData, BasicMode>> traj = std::make_shared<TrajectoryOpt<RosenbrockProblemData, BasicMode>>(problem, seq, builders, decision_builder, true, opts, "snopt");

    traj->initFiniteElements(d, X0, casadi::DM::zeros(si->nx, 1));
    casadi::MXVector sol = traj->optimize();

    traj->setInitialGuess(traj->get_w_sol(), traj->get_lam_x_sol(), traj->get_lam_g_sol());
    sol = traj->optimize();

    std::shared_ptr<opt::solution::Solution> solution_interface_ = std::make_shared<galileo::opt::solution::Solution>();
    solution_interface_->UpdateSolution(traj->getSolutionSegments(), traj->get_w_sol(), traj->get_lam_x_sol(), traj->get_lam_g_sol(), traj->get_f_sol());
    solution_interface_->UpdateConstraints(traj->getConstraintDataSegments());

    Eigen::VectorXd query_times = Eigen::VectorXd::LinSpaced(250, 0., seq->getDT());
    Eigen::MatrixXd state_result = Eigen::MatrixXd::Zero(si->nx, query_times.size());
    Eigen::MatrixXd input_result = Eigen::MatrixXd::Zero(si->nu, query_times.size());
    solution_interface_->GetSolution(query_times, state_result, input_result);

    std::vector<std::vector<galileo::opt::constraint_evaluations_t>> constraints = solution_interface_->GetConstraints(query_times, state_result, input_result);

    opt::solution::solution_t solution(query_times, state_result.transpose(), input_result.transpose());

    std::string plot_dir = "../examples/visualization/plots/";
    std::shared_ptr<tools::GNUPlotInterface> plotting_interface = std::make_shared<galileo::tools::GNUPlotInterface>(plot_dir);

    plotting_interface->PlotSolution(solution, {std::make_tuple(0, si->nx)},
                                     {std::make_tuple(0, si->nu)},
                                     {"States"},
                                     {{"x1", "x2"}},
                                     {"Input"},
                                     {{"u"}});

    plotting_interface->PlotConstraints(constraints);

    return 0;
}