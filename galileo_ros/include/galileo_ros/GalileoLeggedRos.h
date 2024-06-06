#pragma once

#include <galileo/legged-model/LeggedInterface.h>
#include <galileo/legged-model/LeggedModelHelpers.h>
#include <galileo/math/OrientationConversions.h>

#include <ros/ros.h>
#include <ros/package.h>
#include <std_msgs/String.h>
#include <std_msgs/Float64.h>
#include <std_msgs/Float64MultiArray.h>
#include <std_srvs/Trigger.h>

#include "galileo_ros/RobotModel.h"
#include "galileo_ros/ParameterFileLocation.h"
#include "galileo_ros/ContactSequence.h"
#include "galileo_ros/ContactMode.h"
#include "galileo_ros/ContactPhase.h"
#include "galileo_ros/EnvironmentSurface.h"
#include "galileo_ros/GalileoCommand.h"
#include "galileo_ros/SolutionRequest.h"
#include "galileo_ros/WriteSolPlotCons.h"
#include "galileo_ros/InitState.h"

namespace galileo
{
    namespace legged
    {
        class GalileoLeggedRos : public LeggedInterface
        {
        public:
            GalileoLeggedRos(std::shared_ptr<ros::NodeHandle> nh, std::string solver_id, std::string sol_data_dir = "", std::string plot_dir = "") : LeggedInterface(sol_data_dir, plot_dir), nh_(nh), solver_id_(solver_id)
            {
                InitPublishers();
                InitSubscribers();
                InitServices();
            }

            ~GalileoLeggedRos() {}

        private:
            // Initialize ROS

            // Initialize ROS publishers
            void InitPublishers();
            // Initialize ROS subscribers
            void InitSubscribers();
            // Initialize ROS services
            void InitServices();

            // Callbacks
            // Callback for model location subscriber, loads the model and sets the end effectors
            void ModelLocationCallback(const galileo_ros::RobotModel::ConstPtr &msg);
            // Callback for parameter location subscriber, loads the parameters
            void ParameterLocationCallback(const galileo_ros::ParameterFileLocation::ConstPtr &msg);
            // Callback for contact sequence subscriber, sets the contact sequence
            void ContactSequenceCallback(const galileo_ros::ContactSequence::ConstPtr &msg);
            // Callback for environment surface subscriber, adds an environment surface to the library of surfaces
            void SurfaceCallback(const galileo_ros::EnvironmentSurface::ConstPtr &msg);

            // Callback for initialization service, checks if the solver can be initialized
            bool InitStateServiceCallback(galileo_ros::InitState::Request &req, galileo_ros::InitState::Response &res);

            void GeneralCommandCallback(const galileo_ros::GalileoCommand::ConstPtr &msg);

            // Callback for initialization command subscriber, initializes the solver
            void InitializationCallback(const galileo_ros::GalileoCommand::ConstPtr &msg);
            // Callback for update command subscriber, updates the solver
            void UpdateCallback(const galileo_ros::GalileoCommand::ConstPtr &msg);

            bool GetSolutionCallback(galileo_ros::SolutionRequest::Request &req, galileo_ros::SolutionRequest::Response &res);

            bool WriteSolPlotConsCallback(galileo_ros::WriteSolPlotCons::Request &req, galileo_ros::WriteSolPlotCons::Response &res);

            std::shared_ptr<ros::NodeHandle> nh_;
            std::string solver_id_;
            math::OrientationDefinition internal_orientation_representation_ = math::OrientationDefinition::Quaternion;
            math::OrientationDefinition published_orientation_representation_ = math::OrientationDefinition::Quaternion;

            ros::Subscriber model_location_subscriber_;     // Subscriber for model file location and end effector names
            ros::Subscriber parameter_location_subscriber_; // Subscriber for parameter value file location
            ros::Subscriber contact_sequence_subscriber_;   // Subscriber for setting the "contact sequence"
            ros::Subscriber surface_subscriber_;            // Subscriber for adding an environment surface. Unimplemented

            ros::Subscriber command_subscriber_; // Subscriber for galileo commands comprising of an initial and final state, as well as an "update" function

            ros::ServiceServer init_state_service_;          // Service for checking if the solver can be initialized
            ros::ServiceServer get_solution_service_;        // Service for getting the solution at a set of query times
            ros::ServiceServer write_sol_plot_cons_service_; // Service for writing the solution and plotting the constraints
        };

    }
}
