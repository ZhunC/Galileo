#include "galileo/legged-model/LeggedRobotStates.h"

namespace galileo
{
    namespace legged
    {
        legged::LeggedRobotStates::LeggedRobotStates(int nq_, int nv_, legged::contact::RobotEndEffectors ees)
        {
            this->nq = nq_;
            this->nv = nv_;
            this->nx = this->nh + this->ndh + this->nq + this->nv;
            this->ndx = this->nh + this->ndh + 2 * this->nv;
            this->nvju = this->nv - this->nvb;
            this->nu = this->nvju;
            for (auto ee : ees)
            {
                if (ee.second->is_6d)
                {
                    this->frame_id_to_index_range[ee.second->frame_id] = std::make_tuple(this->nF, this->nF + 6);
                    this->nF += 6;
                }
                else
                {
                    this->frame_id_to_index_range[ee.second->frame_id] = std::make_tuple(this->nF, this->nF + 3);
                    this->nF += 3;
                }
            }
            this->nu += this->nF;
            this->nu_general = 6 + this->nvju;

            this->h_index = 0;
            this->dh_index = this->h_index + this->nh;
            this->q_index = this->dh_index + this->ndh;
            this->qj_index = this->q_index + this->nqb;
            this->v_index = this->q_index + this->nq;
            this->vj_index = this->v_index + this->nvb;
            
            this->general_force_index = 0;
            this->general_torque_index = this->general_force_index + 3;
            this->general_vju_index = this->general_torque_index + 3;
        }

        template <typename Sym>
        const Sym LeggedRobotStates::get_ch(const Sym &cx)
        {
            casadi_assert(cx.size1() == this->nx, "Invalid input size");
            return cx(casadi::Slice(0, this->nh));
        }

        template const casadi::DM LeggedRobotStates::get_ch<casadi::DM>(const casadi::DM &cx);
        template const casadi::SX LeggedRobotStates::get_ch<casadi::SX>(const casadi::SX &cx);
        template const casadi::MX LeggedRobotStates::get_ch<casadi::MX>(const casadi::MX &cx);

        template <typename Sym>
        const Sym LeggedRobotStates::get_ch_d(const Sym &cdx)
        {
            casadi_assert(cdx.size1() == this->ndx, "Invalid input size");
            return cdx(casadi::Slice(0, this->nh));
        }

        template const casadi::DM LeggedRobotStates::get_ch_d<casadi::DM>(const casadi::DM &cdx);
        template const casadi::SX LeggedRobotStates::get_ch_d<casadi::SX>(const casadi::SX &cdx);
        template const casadi::MX LeggedRobotStates::get_ch_d<casadi::MX>(const casadi::MX &cdx);

        template <typename Sym>
        const Sym LeggedRobotStates::get_cdh(const Sym &cx)
        {
            casadi_assert(cx.size1() == this->nx, "Invalid input size");
            return cx(casadi::Slice(this->nh, this->nh + this->ndh));
        }

        template const casadi::DM LeggedRobotStates::get_cdh<casadi::DM>(const casadi::DM &cx);
        template const casadi::SX LeggedRobotStates::get_cdh<casadi::SX>(const casadi::SX &cx);
        template const casadi::MX LeggedRobotStates::get_cdh<casadi::MX>(const casadi::MX &cx);

        template <typename Sym>
        const Sym LeggedRobotStates::get_cdh_d(const Sym &cdx)
        {
            casadi_assert(cdx.size1() == this->ndx, "Invalid input size");
            return cdx(casadi::Slice(this->nh, this->nh + this->ndh));
        }

        template const casadi::DM LeggedRobotStates::get_cdh_d<casadi::DM>(const casadi::DM &cdx);
        template const casadi::SX LeggedRobotStates::get_cdh_d<casadi::SX>(const casadi::SX &cdx);
        template const casadi::MX LeggedRobotStates::get_cdh_d<casadi::MX>(const casadi::MX &cdx);

        template <typename Sym>
        const Sym LeggedRobotStates::get_q(const Sym &cx)
        {
            casadi_assert(cx.size1() == this->nx, "Invalid input size");
            return cx(casadi::Slice(this->nh + this->ndh, this->nh + this->ndh + this->nq));
        }

        template const casadi::DM LeggedRobotStates::get_q<casadi::DM>(const casadi::DM &cx);
        template const casadi::SX LeggedRobotStates::get_q<casadi::SX>(const casadi::SX &cx);
        template const casadi::MX LeggedRobotStates::get_q<casadi::MX>(const casadi::MX &cx);

        template <typename Sym>
        const Sym LeggedRobotStates::get_q_d(const Sym &cdx)
        {
            casadi_assert(cdx.size1() == this->ndx, "Invalid input size");
            return cdx(casadi::Slice(this->nh + this->ndh, this->nh + this->ndh + this->nv));
        }

        template const casadi::DM LeggedRobotStates::get_q_d<casadi::DM>(const casadi::DM &cdx);
        template const casadi::SX LeggedRobotStates::get_q_d<casadi::SX>(const casadi::SX &cdx);
        template const casadi::MX LeggedRobotStates::get_q_d<casadi::MX>(const casadi::MX &cdx);

        template <typename Sym>
        const Sym LeggedRobotStates::get_qj(const Sym &cx)
        {
            casadi_assert(cx.size1() == this->nx, "Invalid input size");
            return cx(casadi::Slice(this->nh + this->ndh + this->nqb, this->nh + this->ndh + this->nq));
        }

        template const casadi::DM LeggedRobotStates::get_qj<casadi::DM>(const casadi::DM &cx);
        template const casadi::SX LeggedRobotStates::get_qj<casadi::SX>(const casadi::SX &cx);
        template const casadi::MX LeggedRobotStates::get_qj<casadi::MX>(const casadi::MX &cx);

        template <typename Sym>
        const Sym LeggedRobotStates::get_v(const Sym &cx)
        {
            casadi_assert(cx.size1() == this->nx, "Invalid input size");
            return cx(casadi::Slice(this->nh + this->ndh + this->nq, this->nx));
        }

        template const casadi::DM LeggedRobotStates::get_v<casadi::DM>(const casadi::DM &cx);
        template const casadi::SX LeggedRobotStates::get_v<casadi::SX>(const casadi::SX &cx);
        template const casadi::MX LeggedRobotStates::get_v<casadi::MX>(const casadi::MX &cx);

        template <typename Sym>
        const Sym LeggedRobotStates::get_v_d(const Sym &cdx)
        {
            casadi_assert(cdx.size1() == this->ndx, "Invalid input size");
            return cdx(casadi::Slice(this->nh + this->ndh + this->nv, this->ndx));
        }

        template const casadi::DM LeggedRobotStates::get_v_d<casadi::DM>(const casadi::DM &cdx);
        template const casadi::SX LeggedRobotStates::get_v_d<casadi::SX>(const casadi::SX &cdx);
        template const casadi::MX LeggedRobotStates::get_v_d<casadi::MX>(const casadi::MX &cdx);

        template <typename Sym>
        const Sym LeggedRobotStates::get_vj(const Sym &cx)
        {
            casadi_assert(cx.size1() == this->nx, "Invalid input size");
            return cx(casadi::Slice(this->nh + this->ndh + this->nq + this->nvb, this->nx));
        }

        template const casadi::DM LeggedRobotStates::get_vj<casadi::DM>(const casadi::DM &cx);
        template const casadi::SX LeggedRobotStates::get_vj<casadi::SX>(const casadi::SX &cx);
        template const casadi::MX LeggedRobotStates::get_vj<casadi::MX>(const casadi::MX &cx);

        template <typename Sym>
        const Sym LeggedRobotStates::get_f(const Sym &u, pinocchio::FrameIndex ee_id)
        {
            casadi_assert(u.size1() == this->nu, "Invalid input size");
            Sym tmp = u(casadi::Slice(std::get<0>(this->frame_id_to_index_range[ee_id]), std::get<1>(this->frame_id_to_index_range[ee_id])));
            return tmp(casadi::Slice(0, 3));
        }

        template const casadi::DM LeggedRobotStates::get_f<casadi::DM>(const casadi::DM &u, pinocchio::FrameIndex ee_id);
        template const casadi::SX LeggedRobotStates::get_f<casadi::SX>(const casadi::SX &u, pinocchio::FrameIndex ee_id);
        template const casadi::MX LeggedRobotStates::get_f<casadi::MX>(const casadi::MX &u, pinocchio::FrameIndex ee_id);

        template <typename Sym>
        const Sym LeggedRobotStates::get_tau(const Sym &u, pinocchio::FrameIndex ee_id)
        {
            casadi_assert(u.size1() == this->nu, "Invalid input size");
            Sym tmp = u(casadi::Slice(std::get<0>(this->frame_id_to_index_range[ee_id]), std::get<1>(this->frame_id_to_index_range[ee_id])));
            if (tmp.size1() < 6)
            {
                std::cout << "tau does not exist for ee " << ee_id << "\n";
                return tmp;
            }
            else
                return tmp(casadi::Slice(3, 6));
        }

        template const casadi::DM LeggedRobotStates::get_tau<casadi::DM>(const casadi::DM &u, pinocchio::FrameIndex ee_id);
        template const casadi::SX LeggedRobotStates::get_tau<casadi::SX>(const casadi::SX &u, pinocchio::FrameIndex ee_id);
        template const casadi::MX LeggedRobotStates::get_tau<casadi::MX>(const casadi::MX &u, pinocchio::FrameIndex ee_id);

        template <typename Sym>
        const Sym LeggedRobotStates::get_wrench(const Sym &u, pinocchio::FrameIndex ee_id)
        {
            casadi_assert(u.size1() == this->nu, "Invalid input size");
            return u(casadi::Slice(std::get<0>(this->frame_id_to_index_range[ee_id]), std::get<1>(this->frame_id_to_index_range[ee_id])));
        }

        template const casadi::DM LeggedRobotStates::get_wrench<casadi::DM>(const casadi::DM &u, pinocchio::FrameIndex ee_id);
        template const casadi::SX LeggedRobotStates::get_wrench<casadi::SX>(const casadi::SX &u, pinocchio::FrameIndex ee_id);
        template const casadi::MX LeggedRobotStates::get_wrench<casadi::MX>(const casadi::MX &u, pinocchio::FrameIndex ee_id);

        template <typename Sym>
        const Sym LeggedRobotStates::get_all_wrenches(const Sym &u)
        {
            casadi_assert(u.size1() == this->nu, "Invalid input size");
            return u(casadi::Slice(0, this->nvju));
        }

        template const casadi::DM LeggedRobotStates::get_all_wrenches<casadi::DM>(const casadi::DM &u);
        template const casadi::SX LeggedRobotStates::get_all_wrenches<casadi::SX>(const casadi::SX &u);
        template const casadi::MX LeggedRobotStates::get_all_wrenches<casadi::MX>(const casadi::MX &u);

        template <typename Sym>
        const Sym LeggedRobotStates::get_vju(const Sym &u)
        {
            casadi_assert(u.size1() == this->nu, "Invalid input size");
            return u(casadi::Slice(this->nF, this->nu));
        }

        template const casadi::DM LeggedRobotStates::get_vju<casadi::DM>(const casadi::DM &u);
        template const casadi::SX LeggedRobotStates::get_vju<casadi::SX>(const casadi::SX &u);
        template const casadi::MX LeggedRobotStates::get_vju<casadi::MX>(const casadi::MX &u);

        template <typename Sym>
        const Sym LeggedRobotStates::get_general_forces(const Sym &u_general)
        {
            casadi_assert(u_general.size1() == this->nu_general, "Invalid input size");
            return u_general(casadi::Slice(0, 3));
        }

        template const casadi::DM LeggedRobotStates::get_general_forces<casadi::DM>(const casadi::DM &u_general);
        template const casadi::SX LeggedRobotStates::get_general_forces<casadi::SX>(const casadi::SX &u_general);
        template const casadi::MX LeggedRobotStates::get_general_forces<casadi::MX>(const casadi::MX &u_general);

        template <typename Sym>
        const Sym LeggedRobotStates::get_general_torques(const Sym &u_general)
        {
            casadi_assert(u_general.size1() == this->nu_general, "Invalid input size");
            return u_general(casadi::Slice(3, 6));
        }

        template const casadi::DM LeggedRobotStates::get_general_torques<casadi::DM>(const casadi::DM &u_general);
        template const casadi::SX LeggedRobotStates::get_general_torques<casadi::SX>(const casadi::SX &u_general);
        template const casadi::MX LeggedRobotStates::get_general_torques<casadi::MX>(const casadi::MX &u_general);

        template <typename Sym>
        const Sym LeggedRobotStates::get_general_joint_velocities(const Sym &u_general)
        {
            casadi_assert(u_general.size1() == this->nu_general, "Invalid input size");
            return u_general(casadi::Slice(6, this->nu_general));
        }

        template const casadi::DM LeggedRobotStates::get_general_joint_velocities<casadi::DM>(const casadi::DM &u_general);
        template const casadi::SX LeggedRobotStates::get_general_joint_velocities<casadi::SX>(const casadi::SX &u_general);
        template const casadi::MX LeggedRobotStates::get_general_joint_velocities<casadi::MX>(const casadi::MX &u_general);
    }
}