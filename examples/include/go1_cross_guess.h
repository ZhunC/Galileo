// For generating go1 dataset

#include "galileo/legged-model/LeggedBody.h"
#include "galileo/legged-model/LeggedRobotProblemData.h"
#include "galileo/legged-model/LeggedInterface.h"
#include "galileo/opt/TrajectoryOpt.h"
#include "galileo/tools/MeshcatInterface.h"
#include "galileo/tools/GNUPlotInterface.h"
#include "galileo/math/Quat2Euler.h"

#include "galileo/legged-model/LeggedModelHelpers.h"

using namespace galileo;
using namespace legged;
using namespace opt;
using namespace constraints;
using namespace tools;
using namespace math;

struct Problem {
    int ID;
    Eigen::VectorXd Q;
    Eigen::VectorXd R;
    double K;
    double solveTime;
    double cost;
    casadi::Function solution;

    casadi::DM w;
    casadi::DM lam_x;
    casadi::DM lam_g;
};

struct DataRow {
    int problemID;
    int solutionID;
    double originalSolvetime;
    double currentSolvetime;
    double originalCost;
    double currentCost;
    double similarityScore;
};

const std::string robot_location = "../resources/go1/urdf/go1.urdf";
// const std::string solver_parameter_location = "../resources/go1/Parameters/solver_parameters.txt";
const std::string solver_parameter_location = "../resources/go1/Parameters/solver_parameters_snopt.txt";
const std::string problem_parameter_location = "../resources/go1/Parameters/problem_parameters.txt";
const std::string data_location = "../dataset/";
const std::string csv_location = "../dataset/results.csv";