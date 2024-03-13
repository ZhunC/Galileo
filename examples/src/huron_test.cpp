#include "huron_test.h"

int main(int argc, char **argv)
{
    double q0[] = {
        0, 0, 1.0627, 0, 0, 0, 1,
        0.0000, 0.0000, -0.3207, 0.7572, -0.4365, 0.0000, 0.0000, 0.0000, -0.3207, 0.7572, -0.4365, 0.0000};

    Eigen::Map<ConfigVector> q0_vec(q0, 19);

    LeggedBody robot(huron_location, num_ees, end_effector_names);
    std::shared_ptr<LeggedRobotStates> si = robot.si;

    casadi::DM X0 = casadi::DM::zeros(si->nx, 1);
    int j = 0;
    for (int i = si->nh + si->ndh; i < si->nh + si->ndh + si->nq; ++i)
    {
        X0(i) = q0_vec[j];
        ++j;
    }

    std::shared_ptr<environment::EnvironmentSurfaces> surfaces = std::make_shared<environment::EnvironmentSurfaces>();
    surfaces->push_back(environment::createInfiniteGround());

    std::vector<int> knot_num = {20, 20, 20, 20, 20};
    std::vector<double> knot_time = {0.25, 0.1667, 0.25, 0.1667, 0.25};
    std::vector<int> mask_vec = {0b11, 0b10, 0b11, 0b01, 0b11}; // static walk

    int num_steps = 1;
    for (int i = 0; i < num_steps; ++i)
    {
        for (size_t j = 0; j < mask_vec.size(); ++j)
        {
            contact::ContactMode mode;
            mode.combination_definition = robot.getContactCombination(mask_vec[j]);
            mode.contact_surfaces.resize(mode.combination_definition.size());
            int k = 0;
            for (auto combo : mode.combination_definition)
            {
                if (combo.second)
                    mode.contact_surfaces[k] = 0;
                else
                    mode.contact_surfaces[k] = environment::NO_SURFACE;
                ++k;
            }
            robot.contact_sequence->addPhase(mode, knot_num[j], knot_time[j]);
        }
    }

    // contact::ContactMode final_mode;
    // final_mode.combination_definition = robot.getContactCombination(0b11);
    // final_mode.contact_surfaces = {0, 0};
    // robot.contact_sequence->addPhase(final_mode, 20, 0.05);

    std::cout << "Filling dynamics" << std::endl;
    robot.fillModeDynamics(true);

    casadi::SX cq0(robot.model.nq);
    pinocchio::casadi::copy(q0_vec, cq0);

    /*ETH Weights*/
    Eigen::VectorXd Q_diag(si->ndx);
    Q_diag << 15., 15., 30., 5., 10., 10.,                          /*Centroidal momentum error weights*/
        0., 0., 0., 0., 0., 0.,                                     /*Rate of Centroidal momentum error weights*/
        500., 500., 500., 200., 200., 200.,                         /*Floating base position and orientation (exponential coordinates) error weights*/
        20., 20., 20., 20., 20., 20., 20., 20., 20., 20., 20., 20., /*Joint position error weights*/
        0., 0., 0., 0., 0., 0.,                                     /*Floating base velocity error weights*/
        0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.;             /*Joint velocity error weights*/
    Eigen::MatrixXd Q_mat = Q_diag.asDiagonal();

    Eigen::VectorXd R_diag(si->nu);
    R_diag << 1e-3, 1e-3, 1e-3, 1e-3, 1e-3, 1e-3,                   /*First contact wrench error weights*/
        1e-3, 1e-3, 1e-3, 1e-3, 1e-3, 1e-3,                         /*Second contact wrench error weights*/
        10., 10., 10., 10., 10., 10., 10., 10., 10., 10., 10., 10.; /*Joint acceleration error weights*/
    Eigen::MatrixXd R_mat = R_diag.asDiagonal();

    casadi::SX Q = casadi::SX::zeros(si->ndx, si->ndx);
    casadi::SX R = casadi::SX::zeros(si->nu, si->nu);

    pinocchio::casadi::copy(Q_mat, Q);
    pinocchio::casadi::copy(R_mat, R);

    casadi::SX target_pos = vertcat(casadi::SXVector{q0[0] + 0., q0[1] + 0., q0[2] + 0.});
    casadi::SX target_rot = casadi::SX::eye(3);

    pinocchio::SE3Tpl<galileo::opt::ADScalar, 0> oMf = robot.cdata.oMf[robot.model.getFrameId("base", pinocchio::BODY)];
    auto rot = oMf.rotation();
    auto pos = oMf.translation();
    casadi::SX crot = casadi::SX::zeros(3, 3);
    casadi::SX cpos = casadi::SX::zeros(3, 1);
    pinocchio::casadi::copy(rot, crot);
    pinocchio::casadi::copy(pos, cpos);

    casadi::SX rot_c = casadi::SX::mtimes(crot.T(), target_rot);

    casadi::SX target_error_casadi = vertcat(casadi::SXVector{cpos - target_pos, casadi::SX::inv_skew(rot_c - rot_c.T()) / 2});

    casadi::SX X_ref = casadi::SX(X0);
    X_ref(casadi::Slice(si->nh + si->ndh, si->nh + si->ndh + 3)) = target_pos;
    casadi::SX X_error = vertcat(casadi::SXVector{robot.cx(casadi::Slice(0, si->nh + si->ndh)) - X_ref(casadi::Slice(0, si->nh + si->ndh)),
                                                  target_error_casadi,
                                                  robot.cx(casadi::Slice(si->nh + si->ndh + si->nqb, si->nx)) - X_ref(casadi::Slice(si->nh + si->ndh + si->nqb, si->nx))});

    casadi::SX U_ref = casadi::SX::zeros(si->nu, 1);
    U_ref(5) = 9.81 * robot.cdata.mass[0] / robot.num_end_effectors_;
    U_ref(11) = 9.81 * robot.cdata.mass[0] / robot.num_end_effectors_;
    casadi::SX u_error = robot.cu - U_ref;

    casadi::Function L("L",
                       {robot.cx, robot.cu},
                       {1. * 0.5 * casadi::SX::dot(X_error, casadi::SX::mtimes(Q, X_error)) +
                        1. * 0.5 * casadi::SX::dot(u_error, casadi::SX::mtimes(R, u_error))});

    casadi::Function Phi("Phi",
                         {robot.cx},
                         {1. * casadi::SX::dot(X_error, casadi::SX::mtimes(Q, X_error))});

    casadi::Dict opts;
    opts["ipopt.linear_solver"] = "ma97";
    opts["ipopt.ma97_order"] = "metis";
    opts["ipopt.fixed_variable_treatment"] = "make_constraint";
    opts["ipopt.max_iter"] = 250;

    std::shared_ptr<GeneralProblemData> gp_data = std::make_shared<GeneralProblemData>(robot.fint, robot.fdif, L, Phi);

    std::shared_ptr<ConstraintBuilder<LeggedRobotProblemData>> friction_cone_constraint_builder =
        std::make_shared<FrictionConeConstraintBuilder<LeggedRobotProblemData>>();

    std::shared_ptr<ConstraintBuilder<LeggedRobotProblemData>> velocity_constraint_builder =
        std::make_shared<VelocityConstraintBuilder<LeggedRobotProblemData>>();

    std::shared_ptr<ConstraintBuilder<LeggedRobotProblemData>> contact_constraint_builder =
        std::make_shared<ContactConstraintBuilder<LeggedRobotProblemData>>();

    std::vector<std::shared_ptr<ConstraintBuilder<LeggedRobotProblemData>>> builders = {velocity_constraint_builder, friction_cone_constraint_builder};
    std::shared_ptr<DecisionDataBuilder<LeggedRobotProblemData>> decision_builder = std::make_shared<LeggedDecisionDataBuilder<LeggedRobotProblemData>>();

    std::shared_ptr<LeggedRobotProblemData> legged_problem_data = std::make_shared<LeggedRobotProblemData>(gp_data, surfaces, robot.contact_sequence, si, std::make_shared<ADModel>(robot.cmodel),
                                                                                                           std::make_shared<ADData>(robot.cdata), robot.getEndEffectors(), robot.cx, robot.cu, robot.cdt, X0);

    TrajectoryOpt<LeggedRobotProblemData, contact::ContactMode> traj(legged_problem_data, robot.contact_sequence, builders, decision_builder, opts);

    traj.initFiniteElements(1, X0);

    casadi::MXVector sol = traj.optimize();

    // std::cout << "Total duration: " << robot.contact_sequence->getDT() << std::endl;
    Eigen::VectorXd new_times = Eigen::VectorXd::LinSpaced(250, 0., robot.contact_sequence->getDT());

    solution_t new_sol = solution_t(new_times);
    traj.getSolution(new_sol);

    Eigen::MatrixXd subMatrix = new_sol.state_result.block(si->nh + si->ndh, 0, si->nq, new_sol.state_result.cols());
    MeshcatInterface meshcat("../examples/visualization/solution_data/");
    meshcat.WriteTimes(new_times, "sol_times.csv");
    meshcat.WriteJointPositions(subMatrix, "sol_states.csv");
    meshcat.WriteMetadata(huron_location, q0_vec, "metadata.csv");

    auto cons = traj.getConstraintViolations(new_sol);

    // Collect the data specific to each end effector
    std::vector<tuple_size_t> wrench_indices;
    std::vector<std::vector<std::string>> wrench_legend_names;
    std::vector<std::string> ee_plot_names;
    for (auto ee : robot.getEndEffectors())
    {
        // auto range = si->frame_id_to_index_range[ee.second->frame_id];
        // wrench_indices.push_back(std::make_tuple(std::get<1>(range) - 1, std::get<1>(range)));
        wrench_indices.push_back(si->frame_id_to_index_range[ee.second->frame_id]);
        if (ee.second->is_6d)
            wrench_legend_names.push_back({"F_{x}", "F_{y}", "F_{z}", "\\tau_{x}", "\\tau_{y}", "\\tau_{z}"});
        else
        {
            wrench_legend_names.push_back({"F_{x}", "F_{y}", "F_{z}"});
            // wrench_legend_names.push_back({"F_{z}"});
        }
        ee_plot_names.push_back("Contact Wrench of " + ee.second->frame_name);
    }

    GNUPlotInterface plotter(new_sol, cons);
    plotter.PlotSolution({std::make_tuple(si->nh + si->ndh, si->nh + si->ndh + 3), std::make_tuple(si->nh + si->ndh + 3, si->nh + si->ndh + si->nqb)},
                         wrench_indices,
                         {"Positions", "Orientations"},
                         {{"x", "y", "z"}, {"qx", "qy", "qz", "qw"}},
                         ee_plot_names,
                         wrench_legend_names);
    plotter.PlotConstraints();

    return 0;
}