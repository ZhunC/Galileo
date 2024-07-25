// For generating go1 dataset
#include "go1_test.h"
#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <sstream>
#include <Eigen/Dense>

std::vector<double> extractVector(const std::string &line)
{
    std::vector<double> result;

    // Find the positions of the delimiters
    size_t start = line.find("|(");
    size_t end = line.find(")|vector");

    if (start == std::string::npos || end == std::string::npos || start >= end)
    {
        std::cerr << "Invalid format: vector" << std::endl;
        return result;
    }

    // Extract the substring containing the vector elements
    std::string vectorStr = line.substr(start + 2, end - start - 2);

    // Use a stringstream to parse the elements
    std::stringstream ss(vectorStr);
    std::string item;

    while (std::getline(ss, item, ','))
    {
        result.push_back(std::stod(item));
    }

    return result;
}

double extractDouble(const std::string &line)
{
    // Find the positions of the delimiters
    size_t start = line.find("|");
    size_t end = line.find("|double");

    if (start == std::string::npos || end == std::string::npos || start >= end)
    {
        std::cerr << "Invalid format: double" << std::endl;
        return 0.0; // Return a default value or handle the error appropriately
    }

    // Extract the substring containing the double value
    std::string doubleStr = line.substr(start + 1, end - (start + 1));

    // Convert the extracted string to a double
    double value = std::stod(doubleStr);

    return value;
}

void saveCostParam(const std::string &filename, const galileo::legged::LeggedInterface::CostParameters &CP)
{
    std::ofstream file(filename);
    if (file.is_open())
    {
        file << CP.Q_diag.transpose() << std::endl;
        file << CP.R_diag.transpose() << std::endl;
        file << CP.terminal_weight << std::endl;
        file.close();
    }
    else
    {
        throw std::runtime_error("Unable to open file for writing struct.");
    }
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
    Eigen::VectorXd test_r_diag(24);
    double test_k;

    std::ifstream infile("/ros_ws/src/Galileo/resources/go1/Parameters/solver_parameters.txt");
    std::string line;
    int dataset_size = 3;

    if (infile.is_open())
    {
        while (std::getline(infile, line))
        {
            // Check if the line contains the vector pattern
            if (line.find("cost.Q_diag|") != std::string::npos)
            {
                std::vector<double> vec_q = extractVector(line);
                test_q_diag = Eigen::Map<Eigen::VectorXd>(vec_q.data(), vec_q.size());
            }
            if (line.find("cost.R_diag|") != std::string::npos)
            {
                std::vector<double> vec_r = extractVector(line);
                test_r_diag = Eigen::Map<Eigen::VectorXd>(vec_r.data(), vec_r.size());
            }
            if (line.find("cost.terminal_weight|") != std::string::npos)
            {
                test_k = extractDouble(line);
            }
        }
        infile.close();
    }
    else
    {
        std::cerr << "Unable to open file" << std::endl;
    }

    // Define noise parameters
    double noise_factor = 0.5; // +- 50% noise

    // Random number generator for noise
    std::random_device rd;
    std::mt19937 gen(rd());
    std::normal_distribution<> d(0, 1); // Mean 0, standard deviation 1

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

    for (int i = 0; i < dataset_size; ++i)
    {
        // Vector to hold noisy results
        Eigen::VectorXd noisy_vec_q = test_q_diag;
        Eigen::VectorXd noisy_vec_r = test_r_diag;
        double noisy_k = test_k;

        // Add noise to each entry
        for (int i = 0; i < test_q_diag.size(); ++i)
        {
            double magnitude_q = std::abs(test_q_diag[i]);
            double magnitude_r = std::abs(test_r_diag[i]);
            double noise_q = d(gen) * magnitude_q * noise_factor;
            double noise_r = d(gen) * magnitude_r * noise_factor;
            noisy_vec_q[i] = test_q_diag[i] + noise_q;
            noisy_vec_r[i] = test_q_diag[i] + noise_r;
        }
        double magnitude_k = std::abs(test_k);
        double noise_k = d(gen) * magnitude_k * noise_factor;
        noisy_k = test_k + noise_k;

        // std::cout << "Noisy Vector:\n" << noisy_vec_q.transpose() << std::endl;

        solver_interface.SetQDiag(noisy_vec_q);
        solver_interface.SetRDiag(noisy_vec_r);
        solver_interface.SetK(noisy_k);
        std::cout << "Randomized K, Q, R: " << std::endl;
        std::cout << "K: " << noisy_k << std::endl;
        std::cout << "Q: " << noisy_vec_q.transpose() << std::endl;
        std::cout << "R: " << noisy_vec_r.transpose() << std::endl;

        solver_interface.Update(X0, Xf);

        // Eigen::VectorXd new_times = Eigen::VectorXd::LinSpaced(250, 0., solver_interface.getRobotModel()->contact_sequence->getDT());
        // Eigen::MatrixXd new_states = Eigen::MatrixXd::Zero(solver_interface.states()->nx, new_times.size());
        // Eigen::MatrixXd new_inputs = Eigen::MatrixXd::Zero(solver_interface.states()->nu, new_times.size());
        // solver_interface.GetSolution(new_times, new_states, new_inputs);
        // solver_interface.VisualizeSolutionAndConstraints(new_times, new_states, new_inputs);
        galileo::legged::LeggedInterface::CostParameters costParam = solver_interface.getCostParameters();
        casadi::Function solution = solver_interface.GetTrajectoryFunction();

        std::string base_save_path = "../dataset/Solution_function_";
        std::string base_struct_path = "../dataset/CostParam_";
        std::string solution_save_path = base_save_path + std::to_string(i) + ".casadi";
        std::string costparam_save_path = base_struct_path + std::to_string(i) + ".txt";

        solution.save(solution_save_path);
        saveCostParam(costparam_save_path, costParam);
        std::cout << "Results saved in: " << costparam_save_path << std::endl;
        std::cout << solution << std::endl;

        casadi::DMVector tmp = solution(casadi::DM(0.5));
        std::cout << "States: " << tmp[0] << std::endl;
        std::cout << "Inputs: " << tmp[1] << std::endl;
    }

    return 0;
}