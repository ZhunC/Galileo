// For generating go1 dataset

#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <sstream>
#include <Eigen/Dense>
#include <chrono>
#include "go1_cross_guess.h"



double readDoubleFromFile(const std::string& filename) {
    double read_value = 0.0;
    std::ifstream file(filename);
    std::string line;

    if (file.is_open()) {
        while (std::getline(file, line)) {
            std::size_t pos = line.find("Execution time:");
            if (pos != std::string::npos) {
                std::string numberStr = line.substr(pos + 15); 
                std::stringstream ss(numberStr);
                ss >> read_value;
                break;
            }
        }
        file.close();
    } else {
        std::cerr << "Unable to open file " << filename << std::endl;
    }

    return read_value;
}

void readProblemFromFile(const std::string& filename, Problem& problem) {
    std::ifstream file(filename);
    std::string line;

    if (file.is_open()) {
        // Read the first line into Q
        if (std::getline(file, line)) {
            std::stringstream ss(line);
            std::vector<double> Q;
            double value;
            while (ss >> value) {
                Q.push_back(value);
            }
            problem.Q = Eigen::Map<Eigen::VectorXd>(Q.data(), Q.size());
        }

        // Read the second line into R
        if (std::getline(file, line)) {
            std::stringstream ss(line);
            std::vector<double> R;
            double value;
            while (ss >> value) {
                R.push_back(value);
            }
            problem.R = Eigen::Map<Eigen::VectorXd>(R.data(), R.size());
        }

        // Read the third line into K
        if (std::getline(file, line)) {
            std::stringstream ss(line);
            ss >> problem.K;
        }

        file.close();
    } else {
        std::cerr << "Unable to open file " << filename << std::endl;
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

    int dataset_size = 5;
    int numProblems = dataset_size;
    int numSolutions = dataset_size;


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

    std::vector<Problem> problems(numProblems);

    for (int i = 0; i < numProblems; ++i) {
        std::string paramFilename = data_location + "CostParam_" + std::to_string(i) + ".txt";
        std::string solutionFilename = data_location + "Solution_function_" + std::to_string(i) + ".casadi";
        std::string timeFilename = data_location + "recorded_time_" + std::to_string(i) + ".txt";
        
        problems[i].ID = i; // Assign the ID
        readProblemFromFile(paramFilename, problems[i]);
        problems[i].solution = casadi::Function::load(solutionFilename);
        problems[i].solveTime = readDoubleFromFile(timeFilename);
    }
        


    for (int problemIdx = 0; problemIdx < numProblems; ++problemIdx) {
        for (int solutionIdx = 0; solutionIdx < numSolutions; ++solutionIdx) {
            
            
            solver_interface.SetQDiag(problems[problemIdx].Q);
            solver_interface.SetRDiag(problems[problemIdx].R);
            solver_interface.SetK(problems[problemIdx].K);
            solver_interface.SetInitialGuess(problems[problemIdx].solution);


            auto start = std::chrono::system_clock::now();
            solver_interface.Update(X0, Xf);
            auto end = std::chrono::system_clock::now();
            std::chrono::duration<double, std::milli> elapsed = end - start;


            // Print the timing information (for demonstration purposes)
            std::cout << "Problem " << problemIdx << " with solution " << solutionIdx
                    << " took " << elapsed.count() << " miliseconds, whereas the original problem took "
                    <<  problems[problemIdx].solveTime << " miliseconds." << std::endl;
        }
    }

    return 0;
}