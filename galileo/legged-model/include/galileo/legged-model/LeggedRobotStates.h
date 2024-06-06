#pragma once

#include "galileo/opt/States.h"
#include "galileo/math/OrientationDefinition.h"
#include "galileo/legged-model/EndEffector.h"
#include <pinocchio/autodiff/casadi.hpp>
#include "pinocchio/multibody/model.hpp"

namespace galileo
{
    namespace legged
    {
        /**
         * @brief Simple slicer class for getting state variables.
         *
         */
        class LeggedRobotStates : public opt::States
        {
        public:
            /**
             * @brief Construct a new Legged Robot State object.
             *
             * @param nq_ The number of position variables.
             * @param nv_ The number of velocity variables.
             * @param ees The end effectors.
             * @param orientation_def The internal representation of the orientation.
             */
            LeggedRobotStates(int nq_, int nv_, legged::contact::RobotEndEffectors ees, math::OrientationDefinition orientation_def);

            /**
             * @brief Get momenta: nh x 1.
             *
             * @tparam Sym The type of the input
             * @param cx The input
             * @return const Sym The momenta
             */
            template <class Sym>
            const Sym get_ch(const Sym &cx);

            /**
             * @brief Get momenta delta: nh x 1.
             *
             * @tparam Sym The type of the input
             * @param cdx The input
             * @return const Sym The momenta delta
             */
            template <class Sym>
            const Sym get_ch_d(const Sym &cdx);

            /**
             * @brief Get q: nq x 1.
             *
             * @tparam Sym The type of the input
             * @param cx The input
             * @return const Sym The position
             */
            template <class Sym>
            const Sym get_q(const Sym &cx);

            /**
             * @brief Get q delta: nv x 1.
             *
             * @tparam Sym The type of the input
             * @param cdx The input
             * @return const Sym The position delta
             */
            template <class Sym>
            const Sym get_q_d(const Sym &cdx);

            /**
             * @brief Get qb: nqb x 1.
             *
             * @tparam Sym The type of the input
             * @param cx The input
             * @return const Sym The base coordinates
             */
            template <class Sym>
            const Sym get_qb(const Sym &cx);

            /**
             * @brief Get qbp: nqbp x 1.
             *
             * @tparam Sym The type of the input
             * @param cx The input
             * @return const Sym The base position
             */
            template <class Sym>
            const Sym get_qbp(const Sym &cx);

            /**
             * @brief Get qbo: nqbo x 1.
             *
             * @tparam Sym The type of the input
             * @param cx The input
             * @return const Sym The base orientation
             */
            template <class Sym>
            const Sym get_qbo(const Sym &cx);

            /**
             * @brief Get qj: (nq - nqb) x 1.
             *
             * @tparam Sym The type of the input
             * @param cx The input
             * @return const Sym The joint position
             */
            template <class Sym>
            const Sym get_qj(const Sym &cx);

            /**
             * @brief Get f: 3 x 1.
             *
             * @tparam Sym The type of the input
             * @param u The input
             * @param ee_id The end effector id
             * @return const Sym The force
             */
            template <class Sym>
            const Sym get_f(const Sym &u, pinocchio::FrameIndex ee_id);

            /**
             * @brief Get tau: 3 x 1.
             *
             * @tparam Sym The type of the input
             * @param u The input
             * @param ee_id The end effector id
             * @return const Sym The torque
             */
            template <class Sym>
            const Sym get_tau(const Sym &u, pinocchio::FrameIndex ee_id);

            /**
             * @brief Get the contact wrench: 6 x 1.
             *
             * @tparam Sym The type of the input
             * @param u The input
             * @param ee_id The end effector id
             * @return const Sym The contact wrench
             */
            template <class Sym>
            const Sym get_wrench(const Sym &u, pinocchio::FrameIndex ee_id);

            /**
             * @brief Get all wrenches: nF x 1.
             *
             * @tparam Sym The type of the input
             * @param u The input
             * @return const Sym The wrenches
             */
            template <class Sym>
            const Sym get_all_wrenches(const Sym &u);

            /**
             * @brief Get vju: nvju x 1.
             *
             * @tparam Sym The type of the input
             * @param u The input
             * @return const Sym The joint velocity input
             */
            template <class Sym>
            const Sym get_vju(const Sym &u);

            /**
             * @brief Get the general forces: 3 x 1.
             *
             * @tparam Sym The type of the input
             * @param u_general The input
             * @return const Sym The forces
             */
            template <class Sym>
            const Sym get_general_forces(const Sym &u_general);

            /**
             * @brief Get the general torques: 3 x 1.
             *
             * @tparam Sym The type of the input
             * @param u_general The input
             * @return const Sym The torques
             */
            template <class Sym>
            const Sym get_general_torques(const Sym &u_general);

            /**
             * @brief Get the general joint velocities: nvju x 1.
             *
             * @tparam Sym The type of the input
             * @param u_general The input
             * @return const Sym The joint velocities
             */
            template <class Sym>
            const Sym get_general_joint_velocities(const Sym &u_general);

            /**
             * @brief Map between frame id and its index in the input vector.
             * 
             */
            std::map<pinocchio::FrameIndex, std::tuple<int, int>> frame_id_to_index_range;

            /**
             * @brief The internal representation of the orientation.
             * 
             */
            math::OrientationDefinition orientation_definition;

            /*-----------------------------
            Sizes of certain state elements
            -----------------------------*/

            /**
             * @brief Momenta space dimension.
             *
             */
            static const int nh = 6;

            /**
             * @brief Momenta time derivative offset.
             *
             */
            static const int ndh = 6;

            /**
             * @brief Number of coordinates for the base (number of position variables + orientation variables).
             *
             */
            int nqb = 0;

            /**
             * @brief Number of position variables for the base.
             *
             */
            static const int nqbp = 3;

            /**
             * @brief Number of orientation variables for the base.
             *
             */
            int nqbo = 0;

            /**
             * @brief Number of velocity coordinates for the base.
             *
             */
            static const int nvb = 6;

            /**
             * @brief Number of position variables.
             *
             */
            int nq = 0;

            /**
             * @brief Number of velocity variables.
             *
             */
            int nv = 0;

            /**
             * @brief Number of joint velocity inputs.
             *
             */
            int nvju = 0;

            /**
             * @brief Number of wrench variables
             *
             */
            int nF = 0;

            /**
             * @brief Summed wrenches followed by the joint velocities.
             *
             */
            int nu_general = 0;

            /*----------------------------------------
            Starting indices of certain state elements
            ----------------------------------------*/

            /**
             * @brief Starting index of the momenta.
             * 
             */
            int h_index = 0;

            /**
             * @brief Starting index of the position variables.
             * 
             */
            int q_index = 0;

            /**
             * @brief Starting index of the variables describing the orientation.
             * 
             */
            int qbo_index = 0;

            /**
             * @brief Starting index of the joint position variables.
             * 
             */
            int qj_index = 0;

            /**
             * @brief Starting index of the generalized force variables.
             * 
             */
            int general_force_index = 0;

            /**
             * @brief Starting index of the generalized torque variables.
             * 
             */
            int general_torque_index = 0;

            /**
             * @brief Starting index of the generalized joint velocity inputs.
             * 
             */
            int general_vju_index = 0;
        };

        typedef double Scalar;
        typedef casadi::SX ADScalar;

        typedef pinocchio::ModelTpl<Scalar> Model;
        typedef Model::Data Data;

        typedef pinocchio::ModelTpl<ADScalar> ADModel;
        typedef ADModel::Data ADData;

        typedef Model::ConfigVectorType ConfigVector;
        typedef Model::TangentVectorType TangentVector;

        typedef ADModel::ConfigVectorType ConfigVectorAD;
        typedef ADModel::TangentVectorType TangentVectorAD;
    }
}