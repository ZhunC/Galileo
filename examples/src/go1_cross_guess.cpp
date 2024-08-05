// For generating go1 dataset

#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <sstream>
#include <Eigen/Dense>
#include <chrono>
#include "go1_cross_guess.h"
#include <numeric>
#include <casadi/casadi.hpp>
#include <typeinfo>



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

void appendToCSV(const std::string& filename, const DataRow& row) {
    std::ofstream file(filename, std::ios::app);

    // Check if the file is open
    if (!file.is_open()) {
        std::cerr << "Could not open the file!" << std::endl;
        return;
    }

    // Write the data
    file << row.problemID << ','
         << row.solutionID << ','
         << row.originalSolvetime << ','
         << row.currentSolvetime << ','
         << row.similarityScore << '\n';

    // Close the file
    file.close();
}

int main(int argc, char **argv)
{
    std::vector<std::string> end_effector_names;
    std::vector<int> knot_num;
    std::vector<double> knot_time;
    std::vector<std::vector<galileo::legged::environment::SurfaceID>> contact_surfaces;

    std::vector<double> q0_vec;
    std::vector<double> qf_vec;

    int dataset_size = 100;
    int numProblems = dataset_size;
    int numSolutions = dataset_size;
    // int numSolutions = 5;
    DataRow output; 

    



    galileo::legged::helper::ReadProblemFromParameterFile(problem_parameter_location,
                                                          end_effector_names,
                                                          knot_num,
                                                          knot_time,
                                                          contact_surfaces,
                                                          q0_vec,
                                                          qf_vec);

    

    galileo::legged::LeggedInterface solver_interface;


    double sum = std::accumulate(knot_time.begin(), knot_time.end(), 0.0);
    int interpolation_num = 20;
    double step = sum / (interpolation_num - 1);
    double total_percent_error = 0.0;
    //std::cout << "step is " << step << std::endl;
    //std::cout << "knot time is " << knot_time << std::endl;


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
        // std::cout << solutionFilename << std::endl;
        
        problems[i].ID = i; // Assign the ID
        readProblemFromFile(paramFilename, problems[i]);
        problems[i].solution = casadi::Function::load(solutionFilename);
        problems[i].solveTime = readDoubleFromFile(timeFilename);
    }
        
    // std::cout << problems[3].solution(casadi::DM(0.4)) << std::endl;
    // std::cout << problems[4].solution(casadi::DM(0.4)) << std::endl;
    // std::cout << solutionFilename << std::endl;

    // int testidx = 0;
    // std::cout << problems[testidx].solution(casadi::DM(0.5)) << std::endl;

    for (int problemIdx = 0; problemIdx < numProblems; ++problemIdx) {
        
        solver_interface.SetQDiag(problems[problemIdx].Q);
        solver_interface.SetRDiag(problems[problemIdx].R);
        solver_interface.SetK(problems[problemIdx].K);


        for (int solutionIdx = 0; solutionIdx < numSolutions; ++solutionIdx) {

            double interpolation_point;
            casadi::DMVector prev_sol;
            casadi::DMVector curr_sol;
            casadi::DM diff;
            casadi::DM abs_prev_sol;

            // solutionIdx = testidx;
            // std::cout << "check 1" << std::endl;

            solver_interface.SetInitialGuess(problems[solutionIdx].solution);

            // std::cout << "check 2" << std::endl;


            auto start = std::chrono::system_clock::now();
            solver_interface.Update(X0, Xf);
            std::cout << "check 3" << std::endl;
            auto end = std::chrono::system_clock::now();
            std::chrono::duration<double, std::milli> elapsed = end - start;
            casadi::Function curr_solution_func = solver_interface.GetTrajectoryFunction();


            // Print the timing information (for demonstration purposes)
            std::cout << "Problem " << problemIdx << " with solution " << solutionIdx
                    << " took " << elapsed.count() << " miliseconds, whereas the original problem took "
                    <<  problems[problemIdx].solveTime << " miliseconds." << std::endl;



            for (int k = 0; k < interpolation_num; ++k) {
                interpolation_point = k * step + 1e-5;
                prev_sol = problems[problemIdx].solution(casadi::DM(interpolation_point));
                std::cout << interpolation_point << std::endl;
                // std::cout << k << std::endl;
                curr_sol = curr_solution_func(casadi::DM(interpolation_point));
                diff = casadi::DM::abs(prev_sol[0] - curr_sol[0]);
                //std::cout << prev_sol << std::endl;
                std::cout << curr_sol << std::endl;

                abs_prev_sol = casadi::DM::abs(prev_sol[0]);
                for (int l = 0; l < diff.size1(); ++l) {
                    total_percent_error += diff(l, 0).scalar() / abs_prev_sol(l, 0).scalar();
                    //std::cout << total_percent_error << std::endl;
                }
            }
            double mean_per_diff = total_percent_error / (interpolation_num * diff.size1()); // this tells you how different the solutions are
            
            output.problemID = problemIdx;
            output.solutionID = solutionIdx;
            output.originalSolvetime = problems[problemIdx].solveTime;
            output.currentSolvetime =  elapsed.count();
            output.similarityScore = mean_per_diff; 

            appendToCSV(csv_location, output);


            std::cout << "The two solutions have a mean percentage difference of: " << mean_per_diff << std::endl;



            // std::cout << "The solution is of size: " << diff.size1() << ", " << diff.size2() << std::endl;
            
            /*
            casadi::DM dm_interpolation_points(interpolation_points);
            
            casadi::Function prev_sol_interpolate = problems[problemIdx].solution.map(interpolation_num);

            casadi::DM prev_solution_values = prev_sol_interpolate(dm_interpolation_points.T());
            std::cout << "here" << std::endl;
            casadi::DM abs_prev_sol = casadi::DM::abs(prev_solution_values);
            casadi::DM curr_solution_values = solution(dm_interpolation_points);
            casadi::DM abs_diff = casadi::DM::abs(prev_solution_values - curr_solution_values);

            double total_percent_error;
            for (int i = 0; i < abs_diff.size1(); ++i) {
                for (int j = 0; j < abs_diff.size2(); ++j) {
                    total_percent_error += abs_diff(i, j).scalar() / abs_prev_sol(i, j).scalar();
                }
            }
            double mean_per_diff = total_percent_error / interpolation_num; // this tells you how different the solutions are
            

            std::cout << "The two solutions have a mean absolute difference of:" << mean_per_diff << std::endl;

            */
        }
    }

    return 0;
}