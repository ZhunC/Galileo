//For generating go1 dataset
#include "go1_test.h"
#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <sstream>
#include <Eigen/Dense>

std::vector<double> extractVector(const std::string& line) {
    std::vector<double> result;

    // Find the positions of the delimiters
    size_t start = line.find("|(");
    size_t end = line.find(")|vector");

    if (start == std::string::npos || end == std::string::npos || start >= end) {
        std::cerr << "Invalid format" << std::endl;
        return result;
    }

    // Extract the substring containing the vector elements
    std::string vectorStr = line.substr(start + 2, end - start - 2);

    // Use a stringstream to parse the elements
    std::stringstream ss(vectorStr);
    std::string item;

    while (std::getline(ss, item, ',')) {
        result.push_back(std::stod(item));
    }

    return result;
}

int main(int argc, char **argv)
{
    std::vector<std::string> end_effector_names;
    std::vector<int> knot_num;
    std::vector<double> knot_time;
    std::vector<std::vector<galileo::legged::environment::SurfaceID>> contact_surfaces;

    std::vector<double> q0_vec;
    std::vector<double> qf_vec;
    Eigen::VectorXd test_q_diag(24);
    std::ifstream infile("/ros_ws/src/Galileo/resources/go1/Parameters/solver_parameters.txt");
    std::string line;

    if (infile.is_open()) {
        while (std::getline(infile, line)) {
            // Check if the line contains the vector pattern
            if (line.find("cost.Q_diag|") != std::string::npos) {
                std::vector<double> vec_q = extractVector(line);
                test_q_diag = Eigen::Map<Eigen::VectorXd>(vec_q.data(), vec_q.size());
                std::cout << "test_q_diag is " << vec_q << std::endl;

                // Assuming you only need to process the first matching line
                break;
            }
        }
        infile.close();
    } else {
        std::cerr << "Unable to open file" << std::endl;
    }



    galileo::legged::helper::ReadProblemFromParameterFile(problem_parameter_location,
                                                          end_effector_names,
                                                          knot_num,
                                                          knot_time,
                                                          contact_surfaces,
                                                          q0_vec,
                                                          qf_vec);

    galileo::legged::LeggedInterface solver_interface;

    solver_interface.LoadModel(robot_location, end_effector_names);
    solver_interface.LoadParameters(solver_parameter_location);
    solver_interface.SetQDiag(test_q_diag);

    int nx = solver_interface.states()->nx;
    int q_idx = solver_interface.states()->q_index;

    casadi::DM X0;
    std::vector<double> X0_vec = galileo::legged::helper::getXfromq(solver_interface.states()->nx, q_idx, q0_vec);
    galileo::tools::vectorToCasadi(X0_vec, nx, 1, X0);

    casadi::DM Xf;
    std::vector<double> Xf_vec = galileo::legged::helper::getXfromq(solver_interface.states()->nx, q_idx, qf_vec);
    galileo::tools::vectorToCasadi(Xf_vec, nx, 1, Xf);

    solver_interface.addSurface(environment::createInfiniteGround());

    std::vector<uint> mask_vec = galileo::legged::helper::getMaskVectorFromContactSurfaces(contact_surfaces);

    solver_interface.setContactSequence(
        knot_num,
        knot_time,
        mask_vec,
        contact_surfaces);

    solver_interface.Initialize(X0, Xf);
    solver_interface.Update(X0, Xf);

    // Eigen::VectorXd new_times = Eigen::VectorXd::LinSpaced(250, 0., solver_interface.getRobotModel()->contact_sequence->getDT());
    // Eigen::MatrixXd new_states = Eigen::MatrixXd::Zero(solver_interface.states()->nx, new_times.size());
    // Eigen::MatrixXd new_inputs = Eigen::MatrixXd::Zero(solver_interface.states()->nu, new_times.size());
    // solver_interface.GetSolution(new_times, new_states, new_inputs);
    // solver_interface.VisualizeSolutionAndConstraints(new_times, new_states, new_inputs);

    casadi::Function solution = solver_interface.GetTrajectory();
    std::cout << solution << std::endl;

    auto tmp = solution(casadi::DM(0.5));
    std::cout << tmp[0] << std::endl;
    std::cout << tmp[1] << std::endl;
    solution.save("Solution_function.casadi");

    return 0;
}