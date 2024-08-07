#include "go1_cross_guess.h"
#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <sstream>

std::vector<double> readDoublesFromFile(const std::string& filePath) {
    std::vector<double> numbers;
    std::ifstream file(filePath);

    if (!file.is_open()) {
        std::cerr << "Error opening file: " << filePath << std::endl;
        return numbers;
    }

    std::string line;
    if (std::getline(file, line)) {
        std::stringstream ss(line);
        std::string number;
        while (std::getline(ss, number, ',')) {
            try {
                numbers.push_back(std::stod(number));
            } catch (const std::invalid_argument& e) {
                std::cerr << "Invalid number: " << number << std::endl;
            }
        }
    }

    file.close();
    return numbers;
}

double readTimeFromFile(const std::string &filename)
{
    double read_value = 0.0;
    std::ifstream file(filename);
    std::string line;

    if (file.is_open())
    {
        while (std::getline(file, line))
        {
            std::size_t pos = line.find("Execution time:");
            if (pos != std::string::npos)
            {
                std::string numberStr = line.substr(pos + 15);
                std::stringstream ss(numberStr);
                ss >> read_value;
                break;
            }
        }
        file.close();
    }
    else
    {
        std::cerr << "Unable to open file " << filename << std::endl;
    }

    return read_value;
}

double readCostFromFile(const std::string &filename)
{
    double read_value = 0.0;
    std::ifstream file(filename);
    std::string line;

    if (file.is_open())
    {
        while (std::getline(file, line))
        {
            std::size_t pos = line.find("Total cost is:");
            if (pos != std::string::npos)
            {
                std::string numberStr = line.substr(pos + 14);
                std::stringstream ss(numberStr);
                ss >> read_value;
                break;
            }
        }
        file.close();
    }
    else
    {
        std::cerr << "Unable to open file " << filename << std::endl;
    }

    return read_value;
}

double readFirstNumberFromFile(const std::string &filename) {
    double number = 0.0;
    std::ifstream file(filename);
    std::string line;

    if (file.is_open()) {
        while (std::getline(file, line)) {
            std::stringstream ss(line);
            std::string word;
            while (ss >> word) {
                // Check if the word contains a number
                std::size_t pos = 0;
                while (pos < word.size() && !std::isdigit(word[pos]) && word[pos] != '-' && word[pos] != '.') {
                    ++pos;
                }
                if (pos < word.size()) {
                    // Try to convert the remaining part of the word to a number
                    try {
                        number = std::stod(word.substr(pos));
                        return number; // Return the first number found
                    } catch (const std::invalid_argument &) {
                        // Continue if conversion fails
                    }
                }
            }
        }
        file.close();
    } else {
        std::cerr << "Unable to open file " << filename << std::endl;
    }

    return number;
}

void readProblemFromFile(const std::string& filename, Problem &problem) {
    std::ifstream file(filename);
    if (!file.is_open()) {
        std::cerr << "Unable to open file" << std::endl;
        return;
    }

    std::string line;
    std::vector<double> Q_values, R_values;

    // Read first line
    if (std::getline(file, line)) {
        std::stringstream ss(line);
        std::string number;
        while (std::getline(ss, number, ',')) {
            Q_values.push_back(std::stod(number));
        }
    }

    // Read second line
    if (std::getline(file, line)) {
        std::stringstream ss(line);
        std::string number;
        while (std::getline(ss, number, ',')) {
            R_values.push_back(std::stod(number));
        }
    }

    // Read third line
    if (std::getline(file, line)) {
        problem.K = std::stod(line);
    }

    // Convert vectors to Eigen::VectorXd
    problem.Q = Eigen::VectorXd::Map(Q_values.data(), Q_values.size());
    problem.R = Eigen::VectorXd::Map(R_values.data(), R_values.size());

    file.close();
}

void appendToCSV(const std::string &filename, const DataRow &row)
{
    std::ofstream file(filename, std::ios::app);

    // Check if the file is open
    if (!file.is_open())
    {
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

    int dataset_size = 1;
    int numProblems = dataset_size;
    int numSolutions = dataset_size;
    DataRow output;

    double sum = std::accumulate(knot_time.begin(), knot_time.end(), 0.0);
    int interpolation_num = 20;
    double step = sum / (interpolation_num - 1);

    std::vector<Problem> problems(numProblems);

    for (int i = 0; i < numProblems; ++i)
    {
        std::string paramFilename = data_location + "CostParam_" + std::to_string(i) + ".txt";
        std::string solutionFilename = data_location + "Solution_function_" + std::to_string(i) + ".casadi";
        std::string timeFilename = data_location + "recorded_time_" + std::to_string(i) + ".txt";
        std::string costFilename = data_location + "recorded_cost_" + std::to_string(i) + ".txt";

        std::string wFilename = data_location + "primal_w_solution_" + std::to_string(i) + ".txt";
        std::string lam_xFilename = data_location + "dual_lam_x_solution_" + std::to_string(i) + ".txt";
        std::string lam_gFilename = data_location + "dual_lam_g_solution_" + std::to_string(i) + ".txt";

        problems[i].ID = i; // Assign the ID
        readProblemFromFile(paramFilename, problems[i]);
        problems[i].solution = casadi::Function::load(solutionFilename);
        problems[i].solveTime = readFirstNumberFromFile(timeFilename);
        problems[i].cost = readFirstNumberFromFile(costFilename);

        std::vector<double> w_vec = readDoublesFromFile(wFilename);
        std::vector<double> lam_x_vec = readDoublesFromFile(lam_xFilename);
        std::vector<double> lam_g_vec = readDoublesFromFile(lam_gFilename);

        tools::vectorToCasadi(w_vec, w_vec.size(), 1, problems[i].w);
        tools::vectorToCasadi(lam_x_vec, lam_x_vec.size(), 1, problems[i].lam_x);
        tools::vectorToCasadi(lam_g_vec, lam_g_vec.size(), 1, problems[i].lam_g);
    }

    for (int problemIdx = 0; problemIdx < numProblems; ++problemIdx)
    {
        // std::cout << "Q: " << problems[problemIdx].Q << std::endl;
        // std::cout << "R: " << problems[problemIdx].R << std::endl;
        // std::cout << "K: " << problems[problemIdx].K << std::endl;
        solver_interface.SetQDiag(problems[problemIdx].Q);
        solver_interface.SetRDiag(problems[problemIdx].R);
        solver_interface.SetK(problems[problemIdx].K);

        for (int solutionIdx = 0; solutionIdx < numSolutions; ++solutionIdx)
        {

            double interpolation_point;
            casadi::DMVector prev_sol;
            casadi::DMVector curr_sol;
            casadi::DM diff;
            casadi::DM abs_prev_sol;

            // solver_interface.SetInitialGuess(problems[solutionIdx].solution);
            solver_interface.SetInitialGuess(problems[problemIdx].w, problems[problemIdx].lam_x, problems[problemIdx].lam_g);

            auto start = std::chrono::system_clock::now();
            solver_interface.Update(X0, Xf);
            auto end = std::chrono::system_clock::now();
            std::chrono::duration<double, std::milli> elapsed = end - start;
            casadi::Function curr_solution_func = solver_interface.GetTrajectoryFunction();
            double curr_cost = solver_interface.GetFSol().get_elements()[0];

            // Print the timing information (for demonstration purposes)
            std::cout << "Problem " << problemIdx << " with solution " << solutionIdx
                      << " took " << elapsed.count() << " miliseconds, whereas the original problem took "
                      << problems[problemIdx].solveTime << " miliseconds." << std::endl;

            double total_percent_error = 0.0;

            for (int k = 0; k < interpolation_num; ++k)
            {
                interpolation_point = k * step + 1e-5;
                prev_sol = problems[problemIdx].solution(casadi::DM(interpolation_point));
                curr_sol = curr_solution_func(casadi::DM(interpolation_point));
                diff = casadi::DM::abs(prev_sol[0] - curr_sol[0]);

                abs_prev_sol = casadi::DM::abs(prev_sol[0]);
                for (int l = 0; l < diff.size1(); ++l)
                {
                    if (abs_prev_sol(l, 0).scalar() == 0.0)
                    {
                        continue;
                    }
                    total_percent_error += diff(l, 0).scalar() / abs_prev_sol(l, 0).scalar();
                }
            }
            double mean_per_diff = total_percent_error / (interpolation_num * diff.size1()); // this tells you how different the solutions are
            if (abs(total_percent_error) <= 1e-3)
            {
                mean_per_diff = 0.0;
            }

            output.problemID = problemIdx;
            output.solutionID = solutionIdx;
            output.originalSolvetime = problems[problemIdx].solveTime;
            output.currentSolvetime = elapsed.count();
            output.similarityScore = mean_per_diff;

            appendToCSV(csv_location, output);

            std::cout << "The two solutions have a mean percentage difference of: " << mean_per_diff << std::endl;
            std::cout << "Original solution had cost of " << problems[problemIdx].cost << " and the new solution has cost of " << curr_cost << std::endl;

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