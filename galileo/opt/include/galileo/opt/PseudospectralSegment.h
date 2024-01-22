#pragma once

#include "galileo/opt/Segment.h"

using namespace casadi;

namespace galileo
{
    namespace opt
    {

        /**
         * @brief Helper class for storing polynomial information.
         *
         */
        class LagrangePolynomial
        {
        public:
            /**
             * @brief Construct a new Lagrange Polynomial object.
             *
             */
            LagrangePolynomial(){};

            /**
             * @brief Construct a new Lagrange Polynomial object. Compute and store the coefficients for a given degree and collocation scheme.
             *
             * @param d_ Degree of the polynomial
             * @param scheme Collocation scheme: "radau" or "legendre"
             */
            LagrangePolynomial(int d_, const std::string &scheme = "radau");

            /**
             * @brief Perform symbolic Lagrange Interpolation, which, given a time from the Lagrange time scale, interpolates terms to find the value at time t.
             *
             * @param t Time to interpolate at
             * @param terms Terms at knot points to use for interpolation
             * @return const Scalar Resultant expression for the symbolic interpolated value
             */
            template <class Scalar>
            Scalar lagrange_interpolation(double t, const std::vector<Scalar> terms) const;

            /**
             * @brief Degree of the polynomial.
             *
             */
            int d;

            /**
             * @brief The roots of the polynomial.
             *
             */
            std::vector<double> tau_root;

            /**
             * @brief Quadrature coefficients.
             *
             */
            std::vector<double> B;

            /**
             * @brief Collocation coeffficients.
             *
             */
            std::vector<std::vector<double>> C;

            /**
             * @brief Continuity coefficients.
             *
             */
            std::vector<double> D;
        };

        /**
         * @brief PseudospectalSegment class.
         *
         */
        template <class ProblemData>
        class PseudospectralSegment : public Segment<ProblemData>
        {
        public:
            /**
             * @brief Construct a new Pseudospectral Segment object.
             *
             * @param d Polynomial degree
             * @param knot_num_ Number of knots in the segment
             * @param h_ Period of each knot segment
             * @param global_times_ Vector of glbal times
             * @param st_m_ Pointer to the state indices helper
             * @param Fint_ Integrator function
             * @param Fdif_ Difference function
             */
            PseudospectralSegment(int d, int knot_num_, double h_, std::shared_ptr<casadi::DM> global_times_, std::shared_ptr<States> st_m_, Function &Fint_, Function &Fdif_, DM x0_global_, MX x0_local_, Function &F, Function &L, std::vector<std::shared_ptr<ConstraintData>> G, std::shared_ptr<DecisionData> Wx, std::shared_ptr<DecisionData> Wu, MX &J0, MXVector &w, MXVector &g);

            /**
             * @brief Initialize the relevant expressions.
             *
             * @param d Polynomial degree
             */
            void initialize_expression_variables(int d);

            /**
             * @brief Initialize the vector of local times which constraints are evaluated at.
             *
             */
            void initialize_local_time_vector();

            /**
             * @brief Initialize the vector of times which coincide to the decision variables U occur at.
             *
             */
            void initialize_u_time_vector();

            /**
             * @brief Create all the knot segments.
             *
             * @param x0_global Global starting state to integrate from (used for initial guess)
             * @param x0_local Local starting state to integrate from
             *
             */
            void initialize_knot_segments(DM x0_global, MX x0_local);

            /**
             * @brief Build the function graph.
             *
             * @param F Function for the system dynamics
             * @param L Integrated cost
             * @param G Vector of constraint data
             * @param Wx Decision bound and initial guess data for the state
             * @param Wu Decision bound and initial guess data for the input
             */
            void initialize_expression_graph(Function &F, Function &L, std::vector<std::shared_ptr<ConstraintData>> G, std::shared_ptr<DecisionData> Wx, std::shared_ptr<DecisionData> Wu);

            /**
             * @brief Evaluate the expressions with the actual decision variables.
             *
             * @param J0 Accumulated cost so far
             * @param w Decision variable vector to fill
             * @param g Constraint vector to fill
             */
            void evaluate_expression_graph(MX &J0, MXVector &w, MXVector &g);

            /**
             * @brief Extract the solution from the decision variable vector.
             *
             * @param w Decision variable vector
             * @return MX Solution values
             */
            MXVector extract_solution(MX &w) const;

            /**
             * @brief Get the initial state.
             *
             * @return MX The initial state
             */
            MX get_initial_state() const;

            /**
             * @brief Get the initial state deviant.
             *
             * @return MX The initial state deviant
             */
            MX get_initial_state_deviant() const;

            /**
             * @brief Get the final state deviant.
             *
             * @return MX The final state deviant
             */
            MX get_final_state_deviant() const;

            /**
             * @brief Get the actual final state.
             *
             * @return MX The final state.
             */
            MX get_final_state() const;

            /**
             * @brief Get the global times vector.
             *
             * @return std::shared_ptr<casadi::DM> The global times vector
             */
            std::shared_ptr<casadi::DM> get_global_times() const;

            /**
             * @brief Fills the lower bounds on decision variable (lbw) and upper bounds on decision variable (ubw) vectors with values.
             *
             * This function takes in two vectors, lbw and ubw, and fills them with values.
             * The filled values represent the constraint lower and upper bounds on decision variables.
             *
             * @param lbw The vector to be filled with lower bound values on decision variables.
             * @param ubw The vector to be filled with upper bound values on decision variables.
             */
            void fill_lbw_ubw(std::vector<double> &lbw, std::vector<double> &ubw);

            /**
             * @brief Fills the lower bounds on general constraints (lbg) and upper bounds on general constraints (ubg) vectors with values.
             *
             * This function takes in two vectors, lbg and ubg, and fills them with values.
             * The filled values represent the general constraint lower and upper bounds.
             *
             * @param lbg The vector to be filled with general lower bound values.
             * @param ubg The vector to be filled with general upper bound values.
             */
            void fill_lbg_ubg(std::vector<double> &lbg, std::vector<double> &ubg);

            /**
             * @brief Fills the initial guess vector (w0) with values.
             *
             * This function takes in a vector, w0, and fills it with values.
             * The filled values represent the initial guess for the decision variables.
             *
             * @param w0 The vector to be filled with initial guess values.
             */
            void fill_w0(std::vector<double> &w0) const;

            /**
             * @brief Returns the starting and ending index in w.
             *
             * @return tuple_size_t The range of indices
             */
            tuple_size_t get_range_idx_decision_variables() const;

            /**
             * @brief Returns the starting and ending index in g (call after evaluate_expression_graph!).
             *
             * @return tuple_size_t The range of indices
             */
            tuple_size_t get_range_idx_constraint_expressions() const;

            /**
             * @brief Returns the starting and ending index in lbg/ubg.
             *
             * @return tuple_size_t The range of indices
             */
            tuple_size_t get_range_idx_constraint_bounds() const;

            /**
             * @brief Returns the starting and ending index in lbw/ubw.
             *
             * @return tuple_size_t The range of indices
             */
            tuple_size_t get_range_idx_decision_bounds() const;

        private:
            /**
             * @brief Helper function to process a vector of type MX.
             *
             * This function takes a vector of type MX and performs some processing on it.
             * It creates a temporary vector by copying the input vector and removing the last element.
             * The modified vector is then returned.
             *
             * @param vec The input vector of type MX.
             * @return The processed vector of type MX.
             */
            MX processVector(MXVector &vec) const;

            /**
             * @brief Process the offset vector by removing the first element and concatenating the remaining elements horizontally.
             *
             * @param vec The input offset vector.
             * @return The processed offset vector.
             */
            MX processOffsetVector(MXVector &vec) const;

            /**
             * @brief Collocation state decision variables.
             *
             */
            MXVector dXc_var_vec;

            /**
             * @brief Collocation input decision variables.
             *
             */
            MXVector U_var_vec;

            /**
             * @brief Collocation input decision expressions at the state collocation points.
             *
             * Decision variables of control and state are potentially approximated by different degree polynomials.
             *
             */
            SXVector U_at_c_vec;

            /**
             * @brief Collocation state decision expressions at the collocation points.
             *
             */
            SXVector x_at_c_vec;

            /**
             * @brief Knot point deviants state decision variables.
             *
             */
            MXVector dX0_var_vec;

            /**
             * @brief Knot point state expressions (integral functions of the deviants).
             *
             */
            MXVector X0_var_vec;

            /**
             * @brief Function map for extracting the solution from the ocp solution vector.
             *
             */
            Function sol_map_func;

            /**
             * @brief Solution function.
             *
             */
            Function get_sol_func;

            /**
             * @brief Implicit discrete-time function map. This function map returns the vector of collocation equations
                necessary to match the derivative defect between the approximated dynamics and actual system
                dynamics.
             *
             */
            Function collocation_constraint_map;

            /**
             * @brief Implicit discrete-time function map. The map which matches the approximated final state expression with the initial
                state of the next segment.
             *
             */
            Function xf_constraint_map;

            /**
             * @brief Implicit discrete-time function map. The accumulated cost across all the knot segments found using quadrature rules.
             *
             */
            Function q_cost_fold;

            /**
             * @brief User defined constraints, which are functions with certain bounds associated with them.
             *
             */
            std::vector<Function> general_constraint_maps;

            /**
             * @brief Lower bounds associated with the general constraint maps.
             *
             */
            DM general_lbg;

            /**
             * @brief Upper bounds associated with the general constraint maps.
             *
             */
            DM general_ubg;

            /**
             * @brief Lower bounds associated with the decision variable constraint maps.
             *
             */
            DM general_lbw;

            /**
             * @brief Upper bounds associated with the decision variable constraint maps.
             *
             */
            DM general_ubw;

            /**
             * @brief Initial guess for associated with this segment
             */
            DM w0;

            /**
             * @brief Integrator function.
             *
             */
            Function Fint;

            /**
             * @brief Difference function.
             *
             */
            Function Fdif;

            /**
             * @brief Collocation states used to build the expression graphs.
             *
             */
            SXVector dXc;

            /**
             * @brief Collocation inputs used to build the expression graphs.
             *
             */
            SXVector Uc;

            /**
             * @brief Knot states deviants used to build the expression graphs.
             *
             */
            SX dX0;

            /**
             * @brief Knot states used to build the expression graphs.
             *
             */
            SX X0;

            /**
             * @brief Accumulator expression used to build the expression graphs.
             *
             */
            SX Lc;

            /**
             * @brief Input polynomial. Helper object to store polynomial information for the input.
             *
             */
            LagrangePolynomial U_poly;

            /**
             * @brief State polynomial. Helper object to store polynomial information for the state.
             *
             */
            LagrangePolynomial dX_poly;

            /**
             * @brief Helper for indexing the state variables.
             *
             */
            std::shared_ptr<States> st_m;

            /**
             * @brief Number of knot segments.
             *
             */
            int knot_num;

            /**
             * @brief Vector of global times including knot points.
             *
             */
            std::shared_ptr<DM> global_times;

            /**
             * @brief Vector of local times including knot points. Note that this coincides with the times for the decision variables of x.
             *
             */
            DM local_times;

            /**
             * @brief Vector of unique times for this segment, including terminal but not including initial.
             * Used for bound constraint evaluation, so that we don't overconstrain variables which occur at the same time
             * e.g, x0 and xf of adjacent segments
             *
             * In the frame of the global times (NOT the local times)
             *
             */
            DM collocation_times;

            /**
             * @brief Vector of times for the knot points for this segment.
             *
             * In the frame of the global times (NOT the local times)
             */
            DM knot_times;

            /**
             * @brief Vector of times for the decision variables of u.
             */
            DM u_times;
        };

        LagrangePolynomial::LagrangePolynomial(int d_, const std::string &scheme)
        {
            this->d = d_;
            /*Choose collocation points*/
            this->tau_root = collocation_points(this->d, scheme);
            this->tau_root.insert(this->tau_root.begin(), 0);

            /*Coefficients of the quadrature function*/
            this->B.resize(this->d + 1);

            /*Coefficients of the collocation equation*/
            this->C.resize(this->d + 1);
            for (int j = 0; j < this->d + 1; ++j)
                this->C[j].resize(this->d + 1);

            /*Coefficients of the continuity equation*/
            this->D.resize(this->d + 1);

            /*For all collocation points*/
            for (int j = 0; j < this->d + 1; ++j)
            {
                /*Construct Lagrange polynomials to get the polynomial basis at the collocation point*/
                Polynomial p = 1;
                for (int r = 0; r < this->d + 1; ++r)
                {
                    if (r != j)
                    {
                        p *= Polynomial(-this->tau_root[r], 1) / (this->tau_root[j] - this->tau_root[r]);
                    }
                }
                /*Evaluate the polynomial at the final time to get the coefficients of the continuity equation*/
                this->D[j] = p(1.0);

                /*Evaluate the time derivative of the polynomial at all collocation points to get the coefficients of the continuity equation*/
                Polynomial dp = p.derivative();
                for (int r = 0; r < this->d + 1; ++r)
                {
                    this->C[j][r] = dp(this->tau_root[r]);
                }
                Polynomial pint = p.anti_derivative();
                this->B[j] = pint(1.0);
            }
        }

        template <class Scalar>
        Scalar LagrangePolynomial::lagrange_interpolation(double t, const std::vector<Scalar> terms) const
        {
            assert((t >= 0.0) && (t <= 1.0) && "t must be in the range [0,1]");

            SX result = 0;
            for (int j = 0; j < this->d; ++j)
            {
                SX term = terms[j];
                for (int r = 0; r < this->d + 1; ++r)
                {
                    if (r != j)
                    {
                        term *= (t - this->tau_root[r]) / (this->tau_root[j] - this->tau_root[r]);
                    }
                }
                result += term;
            }
            return result;
        }

        template <class ProblemData>
        PseudospectralSegment<ProblemData>::PseudospectralSegment(int d, int knot_num_, double h_, std::shared_ptr<casadi::DM> global_times_, std::shared_ptr<States> st_m_, Function &Fint_, Function &Fdif_, DM x0_global_, MX x0_local_, Function &F, Function &L, std::vector<std::shared_ptr<ConstraintData>> G, std::shared_ptr<DecisionData> Wx, std::shared_ptr<DecisionData> Wu, MX &J0, MXVector &w, MXVector &g)
        {
            assert(d > 0 && d < 10 && "d must be greater than 0 and less than 10");
            assert(h_ > 0 && "h must be a positive duration");
            assert(Fint_.n_in() == 3 && "Fint must have 3 inputs");
            assert(Fint_.n_out() == 1 && "Fint must have 1 output");
            assert(Fdif_.n_in() == 3 && "Fdif must have 3 inputs");
            assert(Fdif_.n_out() == 1 && "Fdif must have 1 output");

            Fint_.assert_size_in(0, st_m_->nx, 1);
            Fint_.assert_size_in(1, st_m_->ndx, 1);
            Fint_.assert_size_in(2, 1, 1);
            Fint_.assert_size_out(0, st_m_->nx, 1);

            Fdif_.assert_size_in(0, st_m_->nx, 1);
            Fdif_.assert_size_in(1, st_m_->nx, 1);
            Fdif_.assert_size_in(2, 1, 1);
            Fdif_.assert_size_out(0, st_m_->ndx, 1);

            this->knot_num = knot_num_;
            this->h = h_;
            this->global_times = global_times_;
            this->st_m = st_m_;
            this->Fint = Fint_;
            this->Fdif = Fdif_;
            this->T = (this->knot_num) * this->h;

            this->initialize_expression_variables(d);

            this->initialize_local_time_vector();
            this->initialize_u_time_vector();
            this->initialize_knot_segments(x0_global_, x0_local_);
            this->initialize_expression_graph(F, L, G, Wx, Wu);
            this->evaluate_expression_graph(J0, w, g);
        }

        template <class ProblemData>
        void PseudospectralSegment<ProblemData>::initialize_expression_variables(int d)
        {
            this->dXc.clear();
            this->Uc.clear();

            this->dX_poly = LagrangePolynomial(d);
            this->U_poly = LagrangePolynomial(d - 1);

            for (int j = 0; j < this->dX_poly.d; ++j)
            {
                this->dXc.push_back(SX::sym("dXc_" + std::to_string(j), this->st_m->ndx, 1));
                if (j < this->U_poly.d)
                {
                    this->Uc.push_back(SX::sym("Uc_" + std::to_string(j), this->st_m->nu, 1));
                }
            }
            this->dX0 = SX::sym("dX0", this->st_m->ndx, 1);
            this->X0 = SX::sym("X0", this->st_m->nx, 1);
            this->Lc = SX::sym("Lc", 1, 1);
        }

        template <class ProblemData>
        void PseudospectralSegment<ProblemData>::initialize_local_time_vector()
        {
            this->local_times = casadi::DM::zeros(this->knot_num * (this->dX_poly.d + 1), 1);
            this->collocation_times = casadi::DM::zeros(this->knot_num * (this->dX_poly.d), 1);
            this->knot_times = casadi::DM::zeros(this->knot_num + 1, 1);
            int i = 0;
            int unique_i = 0;
            int j = 0;
            double kh = 0.0;
            for (int k = 0; k < this->knot_num; ++k)
            {
                kh = k * this->h;
                this->local_times(i, 0) = kh;
                this->knot_times(k, 0) = kh;
                ++i;
                for (j = 0; j < this->dX_poly.d; ++j)
                {
                    this->local_times(i, 0) = kh + this->dX_poly.tau_root[j + 1] * this->h;
                    if (i > 0 && unique_i < this->collocation_times.size1())
                    {
                        this->collocation_times(unique_i, 0) = this->local_times(i, 0);
                        ++unique_i;
                    }
                    ++i;
                }
            }
            this->knot_times(this->knot_times.size1() - 2, 0) = this->T - this->h;
            this->knot_times(this->knot_times.size1() - 1, 0) = this->T;

            double start_time = 0.0;
            if (global_times != nullptr)
            {
                start_time = (*this->global_times)(global_times->size1() - 1, 0).get_elements()[0];
                auto tmp = this->local_times + start_time;
                this->global_times = std::make_shared<DM>(vertcat(*this->global_times, tmp));
            }
            else
            {
                this->global_times = std::make_shared<DM>(this->local_times);
            }
            this->collocation_times += start_time;
            this->knot_times += start_time;
        }

        template <class ProblemData>
        void PseudospectralSegment<ProblemData>::initialize_u_time_vector()
        {
            this->u_times = casadi::DM::zeros(this->knot_num * (this->U_poly.d), 1);
            int j = 0;
            double kh = 0.0;
            this->u_times(0, 0) = kh;
            int i = 1;
            for (int k = 0; k < this->knot_num; ++k)
            {
                kh = k * this->h;
                for (j = 0; j < this->U_poly.d; ++j)
                {
                    if (i < this->u_times.size1())
                        this->u_times(i, 0) = kh + this->U_poly.tau_root[j + 1] * this->h;
                    else
                        break;
                    ++i;
                }
            }
        }

        template <class ProblemData>
        void PseudospectralSegment<ProblemData>::initialize_knot_segments(DM x0_global_, MX x0_local_)
        {
            this->x0_global = x0_global_;
            this->x0_local = x0_local_;
            assert(this->x0_local.size1() == this->st_m->nx && this->x0_local.size2() == 1 && "x0 must be a column vector of size nx");

            this->dXc_var_vec.clear();
            this->U_var_vec.clear();
            this->dX0_var_vec.clear();
            this->X0_var_vec.clear();
            this->U_at_c_vec.clear();
            this->x_at_c_vec.clear();
            for (int k = 0; k < this->knot_num; ++k)
            {
                this->dXc_var_vec.push_back(MX::sym("dXc_" + std::to_string(k), this->st_m->ndx * this->dX_poly.d, 1));
                this->U_var_vec.push_back(MX::sym("U_" + std::to_string(k), this->st_m->nu * this->U_poly.d, 1));
            }

            for (int k = 0; k < this->knot_num + 1; ++k)
            {
                this->dX0_var_vec.push_back(MX::sym("dX0_" + std::to_string(k), this->st_m->ndx, 1));
                this->X0_var_vec.push_back(this->Fint(MXVector{this->x0_local, this->dX0_var_vec[k], 1.0}).at(0));
            }
        }

        template <class ProblemData>
        void PseudospectralSegment<ProblemData>::initialize_expression_graph(Function &F, Function &L, std::vector<std::shared_ptr<ConstraintData>> G, std::shared_ptr<DecisionData> Wx, std::shared_ptr<DecisionData> Wu)
        {
            assert(F.n_in() == 2 && "F must have 2 inputs");
            assert(F.n_out() == 1 && "F must have 1 output");

            assert(L.n_in() == 2 && "L must have 2 inputs");
            assert(L.n_out() == 1 && "L must have 1 output");

            F.assert_size_in(0, this->st_m->nx, 1);
            F.assert_size_in(1, this->st_m->nu, 1);
            F.assert_size_out(0, this->st_m->ndx, 1);

            L.assert_size_in(0, this->st_m->nx, 1);
            L.assert_size_in(1, this->st_m->nu, 1);
            L.assert_size_out(0, 1, 1);

            /*Collocation equations*/
            SXVector eq;
            /*State at the end of the collocation interval*/
            SX dXf = this->dX_poly.D[0] * this->dX0;
            /*Cost at the end of the collocation interval*/
            SX Qf = 0;
            /*Actual state at collocation points*/
            SXVector x_at_c;
            /*U interpolated at the dx polynomial collocation points*/
            SXVector u_at_c;
            SXVector tmp_x;
            SXVector tmp_dx;
            tmp_x.push_back(this->X0);
            tmp_dx.push_back(this->dX0);

            for (int j = 0; j < this->dX_poly.d; ++j)
            {
                double dt_j = (this->dX_poly.tau_root[j + 1] - this->dX_poly.tau_root[j]) * this->h;
                /*Expression for the state derivative at the collocation point*/
                SX dxp = this->dX_poly.C[0][j + 1] * this->dX0;
                for (int r = 0; r < this->dX_poly.d; ++r)
                {
                    dxp += this->dX_poly.C[r + 1][j + 1] * this->dXc[r];
                }
                /*dXc must exist in a Euclidean space, but we need x_c in order to evaluate the objective. Fint can simply return dXc[j] if the states are already Euclidean*/
                SX x_c = this->Fint(SXVector{this->X0, this->dXc[j], dt_j}).at(0);
                SX u_c = this->U_poly.lagrange_interpolation(this->dX_poly.tau_root[j], this->Uc);
                // SX u_c = this->Uc[0];

                x_at_c.push_back(x_c);
                u_at_c.push_back(u_c);
                tmp_x.push_back(x_c);
                tmp_dx.push_back(this->dXc[j]);

                /*Append collocation equations*/
                eq.push_back(this->h * F(SXVector{x_c, u_c}).at(0) - dxp);

                /*Add cost contribution*/
                SXVector L_out = L(SXVector{x_c, u_c});
                /*This is fine as long as the cost is not related to the Lie Group elements. See the state integrator and dX for clarity*/
                Qf += this->dX_poly.B[j + 1] * L_out.at(0) * this->h;
                // Qf += this->U_poly.B(j + 1) * L_out.at(1) * this->h;

                dXf += this->dX_poly.D[j + 1] * this->dXc[j];
            }

            casadi::Dict opts;
            // opts["jit"] = true;
            // opts["jit_options.flags"] = "-O3";
            // opts["jit_options.compiler"] = "gcc";
            // opts["compiler"] = "shell";

            auto collocation_constraint = Function("feq",
                                                   SXVector{this->X0, vertcat(this->dXc), this->dX0, vertcat(this->Uc)},
                                                   SXVector{vertcat(eq)}, opts);
            // collocation_constraint_original.generate("feq");
            // int flag1 = system("gcc -fPIC -shared -O3 feq.c -o feq.so");
            // casadi_assert(flag1==0, "Compilation failed");
            // auto collocation_constraint = external("feq");

            // auto collocation_constraint_adj1 = collocation_constraint_original.reverse(1);
            // collocation_constraint_adj1.generate("adj1_feq");
            // int flag2 = system("gcc -fPIC -shared -O3 adj1_feq.c -o adj1_feq.so");
            // casadi_assert(flag2==0, "Compilation failed");
            // auto collocation_constraint_adj = external("adj1_feq");

            auto xf_constraint = Function("fxf",
                                          SXVector{this->X0, vertcat(this->dXc), this->dX0, vertcat(this->Uc)},
                                          SXVector{dXf}, opts);
            // xf_constraint.generate("xf_constraint");
            // int flag2 = system("gcc -fPIC -shared -O3 xf_constraint.c -o xf_constraint.so");
            // casadi_assert(flag2==0, "Compilation failed");
            // xf_constraint = external("xf_constraint");

            auto q_cost = Function("fxq", SXVector{this->Lc, this->X0, vertcat(this->dXc), this->dX0, vertcat(this->Uc)},
                                   SXVector{this->Lc + Qf}, opts);
            // q_cost.generate("q_cost");
            // int flag3 = system("gcc -fPIC -shared -O3 q_cost.c -o q_cost.so");
            // casadi_assert(flag3==0, "Compilation failed");
            // q_cost = external("q_cost");

            /*Implicit discrete-time equations*/
            this->collocation_constraint_map = collocation_constraint.map(this->knot_num, "openmp");
            /*When you evaluate this map, subtract by the knot points list offset by 1 to be correct*/
            this->xf_constraint_map = xf_constraint.map(this->knot_num, "openmp");
            this->q_cost_fold = q_cost.fold(this->knot_num);

            this->sol_map_func = Function("sol_map",
                                          SXVector{this->X0, vertcat(this->dXc), this->dX0, vertcat(this->Uc)},
                                          SXVector{horzcat(tmp_x), horzcat(u_at_c)})
                                     .map(this->knot_num, "serial");

            casadi_int N = this->collocation_constraint_map.size1_out(0) * this->collocation_constraint_map.size2_out(0) +
                           this->xf_constraint_map.size1_out(0) * this->xf_constraint_map.size2_out(0);
            auto tmp = N;

            std::vector<tuple_size_t> ranges_G;

            /*Map the constraint to each collocation point, and then map the mapped constraint to each knot segment*/
            for (std::size_t i = 0; i < G.size(); ++i)
            {
                auto g_data = G[i];

                assert(g_data->G.n_in() == 2 && "G must have 2 inputs");
                g_data->G.assert_size_in(0, this->st_m->nx, 1);
                g_data->G.assert_size_in(1, this->st_m->nu, 1);
                /*TODO: Add assertions to check the bounds functions here!!!*/

                if (g_data->global)
                {
                    auto tmap = Function("fg",
                                         SXVector{this->X0, vertcat(this->dXc), this->dX0, vertcat(this->Uc)},
                                         SXVector{vertcat(g_data->G.map(this->dX_poly.d, "serial")((SXVector{horzcat(x_at_c), horzcat(u_at_c)})))})
                                    .map(this->knot_num, "serial");
                    this->general_constraint_maps.push_back(tmap);
                    ranges_G.push_back(tuple_size_t(N, N + tmap.size1_out(0) * tmap.size2_out(0)));
                    N += tmap.size1_out(0) * tmap.size2_out(0);
                }
                else
                {
                    for (int k = 0; k < g_data->apply_at.rows(); ++k)
                    {
                    }
                }
            }

            this->general_lbg.resize(N, 1);
            this->general_ubg.resize(N, 1);
            this->general_lbg(Slice(0, tmp)) = DM::zeros(tmp, 1);
            this->general_ubg(Slice(0, tmp)) = DM::zeros(tmp, 1);

            for (casadi_int i = 0; i < G.size(); ++i)
            {
                auto g_data = G[i];
                if (g_data->global)
                {
                    this->general_lbg(Slice(std::get<0>(ranges_G[i])), std::get<1>(ranges_G[i])) =
                        vertcat(g_data->lower_bound.map(this->knot_num * (this->dX_poly.d), "serial")(this->collocation_times));
                    this->general_ubg(Slice(std::get<0>(ranges_G[i])), std::get<1>(ranges_G[i])) =
                        vertcat(g_data->upper_bound.map(this->knot_num * (this->dX_poly.d), "serial")(this->collocation_times));
                }
                else
                {
                }
            }

            auto Ndxknot = this->st_m->ndx * (this->knot_num + 1);
            auto Ndx = this->st_m->ndx * (this->dX_poly.d + 1) * this->knot_num + this->st_m->ndx;
            auto Ndxcol = Ndx - Ndxknot;

            auto Nu = this->st_m->nu * this->U_poly.d * this->knot_num;
            this->w0 = DM::zeros(Ndx + Nu, 1);
            this->general_lbw = -DM::inf(Ndx + Nu, 1);
            this->general_ubw = DM::inf(Ndx + Nu, 1);

            /*Transform initial guess for x to an initial guess for dx, using f_dif, the inverse of f_int*/

            MX xkg_sym = casadi::MX::sym("xkg", this->st_m->nx, 1);
            MX xckg_sym = casadi::MX::sym("Xckg", this->st_m->nx * this->dX_poly.d, 1);

            if (!Wx->initial_guess.is_null())
            {
                auto xg = vertcat(Wx->initial_guess.map(this->knot_num + 1, "serial")(this->knot_times));
                Function dxg_func = Function("xg_fun", MXVector{xkg_sym}, MXVector{this->Fdif(MXVector{this->x0_global, xkg_sym, 1.0}).at(0)})
                                        .map(this->knot_num + 1, "serial");
                this->w0(Slice(0, Ndxknot)) = reshape(dxg_func(DMVector{xg}).at(0), Ndxknot, 1);
                /*The transformation of xc to dxc is a slightly less trivial. While x_k = fint(x0_init, dx_k), for xc_k, we have xc_k = fint(x_k, dxc_k) which is equivalent to xc_k = fint(fint(x0_init, dx_k), dxc_k).
                Thus, dxc_k = fdif(fint(x0_init, dx_k), xc_k)). This could be done with maps like above, but it is not necessary.*/
                auto xc_g = vertcat(Wx->initial_guess.map((this->dX_poly.d) * this->knot_num, "serial")(this->collocation_times));
                for (casadi_int i = 0; i < this->knot_num; ++i)
                {
                    auto xk = xg(Slice(i * this->st_m->nx, (i + 1) * this->st_m->nx));
                    auto xck = xc_g(Slice(i * this->st_m->nx * this->dX_poly.d, (i + 1) * this->st_m->nx * this->dX_poly.d));
                    for (casadi_int j = 0; j < this->dX_poly.d; ++j)
                    {
                        this->w0(Slice(Ndxknot + i * this->st_m->ndx * this->dX_poly.d + j * this->st_m->ndx, Ndxknot + i * this->st_m->ndx * this->dX_poly.d + (j + 1) * this->st_m->ndx)) =
                            reshape(this->Fdif(DMVector{xk, xck(Slice(j * this->st_m->nx, (j + 1) * this->st_m->nx)), this->h}).at(0), this->st_m->ndx, 1);
                    }
                }
            }

            if (!Wx->lower_bound.is_null() && !Wx->upper_bound.is_null())
            {
                this->general_lbw(Slice(0, Ndxknot)) = reshape(vertcat(Wx->lower_bound.map(this->knot_num + 1, "serial")(this->knot_times)), Ndxknot, 1);
                this->general_ubw(Slice(0, Ndxknot)) = reshape(vertcat(Wx->upper_bound.map(this->knot_num + 1, "serial")(this->knot_times)), Ndxknot, 1);
                this->general_lbw(Slice(Ndxknot, Ndx)) = reshape(vertcat(Wx->lower_bound.map((this->dX_poly.d) * this->knot_num, "serial")(this->collocation_times)), Ndxcol, 1);
                this->general_ubw(Slice(Ndxknot, Ndx)) = reshape(vertcat(Wx->upper_bound.map((this->dX_poly.d) * this->knot_num, "serial")(this->collocation_times)), Ndxcol, 1);
            }

            if (!Wu->initial_guess.is_null())
            {
                this->w0(Slice(Ndx, Ndx + Nu)) = reshape(vertcat(Wu->initial_guess.map(this->U_poly.d * this->knot_num, "serial")(this->u_times)), Nu, 1);
            }

            if (!Wu->lower_bound.is_null() && !Wu->upper_bound.is_null())
            {
                this->general_lbw(Slice(Ndx, Ndx + Nu)) = reshape(vertcat(Wu->lower_bound.map(this->U_poly.d * this->knot_num, "serial")(this->u_times)), Nu, 1);
                this->general_ubw(Slice(Ndx, Ndx + Nu)) = reshape(vertcat(Wu->upper_bound.map(this->U_poly.d * this->knot_num, "serial")(this->u_times)), Nu, 1);
            }
        }

        template <class ProblemData>
        MX PseudospectralSegment<ProblemData>::processVector(MXVector &vec) const
        {
            MXVector temp = vec;
            temp.pop_back();
            return horzcat(temp);
        }

        template <class ProblemData>
        MX PseudospectralSegment<ProblemData>::processOffsetVector(MXVector &vec) const
        {
            MXVector temp = vec;
            temp.erase(temp.begin());
            return horzcat(temp);
        }

        template <class ProblemData>
        void PseudospectralSegment<ProblemData>::evaluate_expression_graph(MX &J0, MXVector &w, MXVector &g)
        {
            assert(J0.size1() == 1 && J0.size2() == 1 && "J0 must be a scalar");

            MXVector result;

            MX xs = this->processVector(this->X0_var_vec);
            MX dxs = this->processVector(this->dX0_var_vec);
            MX dxcs = horzcat(this->dXc_var_vec);
            MX us = horzcat(this->U_var_vec);
            MX xs_offset = this->processOffsetVector(this->X0_var_vec);
            MX dxs_offset = this->processOffsetVector(this->dX0_var_vec);

            MXVector solmap_restult = this->sol_map_func(MXVector{xs, dxcs, dxs, us});
            MX all_xs = solmap_restult.at(0);
            MX all_us = solmap_restult.at(1);

            /*This section cannot get much faster, it is bounded by the time to evaluate the constraint*/
            MX col_con_mat = this->collocation_constraint_map(MXVector{xs, dxcs, dxs, us}).at(0);
            MX xf_con_mat = this->xf_constraint_map(MXVector{xs, dxcs, dxs, us}).at(0);
            dxs_offset = reshape(dxs_offset, dxs_offset.size1() * dxs_offset.size2(), 1);

            result.push_back(reshape(col_con_mat, col_con_mat.size1() * col_con_mat.size2(), 1));
            result.push_back(reshape(xf_con_mat, xf_con_mat.size1() * xf_con_mat.size2(), 1) -
                             dxs_offset);

            for (std::size_t i = 0; i < this->general_constraint_maps.size(); ++i)
            {
                MX g_con_mat = this->general_constraint_maps[i](MXVector{xs, dxcs, dxs, us}).at(0);
                result.push_back(reshape(g_con_mat, g_con_mat.size1() * g_con_mat.size2(), 1));
            }

            MX cost = this->q_cost_fold(MXVector{J0, xs, dxcs, dxs, us}).at(0);
            J0 = cost;
            /*where g of this segment starts*/
            std::size_t g_size = g.size();
            /*Use std::move to avoid copying the vectors. Reserve space for g in advance outside of PseudospectralSegment.*/
            g.insert(g.end(), std::make_move_iterator(result.begin()), std::make_move_iterator(result.end()));
            this->g_range = tuple_size_t(g_size, g.size());

            /*where w of this segment starts*/
            std::size_t w_size = w.size();
            /*Use std::move to avoid copying the vectors. Reserve space for w in advance outside of PseudospectralSegment.*/
            w.insert(w.end(), std::make_move_iterator(this->dX0_var_vec.begin()), std::make_move_iterator(this->dX0_var_vec.end()));
            w.insert(w.end(), std::make_move_iterator(this->dXc_var_vec.begin()), std::make_move_iterator(this->dXc_var_vec.end()));
            w.insert(w.end(), std::make_move_iterator(this->U_var_vec.begin()), std::make_move_iterator(this->U_var_vec.end()));

            this->w_range = tuple_size_t(w_size, std::accumulate(w.begin() + w_size, w.end(), 0, [](int sum, const MX &item)
                                                                 { return sum + item.size1() * item.size2(); }));
            this->get_sol_func = Function("func",
                                          MXVector({vertcat(MXVector(w.begin() + w_size, w.begin() + w.size()))}),
                                          MXVector({all_xs, all_us}));
        }

        template <class ProblemData>
        MXVector PseudospectralSegment<ProblemData>::extract_solution(MX &w) const
        {
            return this->get_sol_func(MXVector{w(Slice(casadi_int(std::get<0>(this->w_range)), casadi_int(std::get<1>(this->w_range))))});
        }

        template <class ProblemData>
        MX PseudospectralSegment<ProblemData>::get_initial_state_deviant() const
        {
            return this->dX0_var_vec.front();
        }

        template <class ProblemData>
        MX PseudospectralSegment<ProblemData>::get_initial_state() const
        {
            return this->X0_var_vec.front();
        }

        template <class ProblemData>
        MX PseudospectralSegment<ProblemData>::get_final_state_deviant() const
        {
            return this->dX0_var_vec.back();
        }

        template <class ProblemData>
        MX PseudospectralSegment<ProblemData>::get_final_state() const
        {
            return this->X0_var_vec.back();
        }

        template <class ProblemData>
        std::shared_ptr<casadi::DM> PseudospectralSegment<ProblemData>::get_global_times() const
        {
            return this->global_times;
        }

        template <class ProblemData>
        void PseudospectralSegment<ProblemData>::fill_lbw_ubw(std::vector<double> &lbw, std::vector<double> &ubw)
        {
            /*where lb/ub of this segment starts*/
            auto bw_size = lbw.size();
            std::vector<double> element_access1 = this->general_lbw.get_elements();
            std::vector<double> element_access2 = this->general_ubw.get_elements();

            lbw.insert(lbw.end(), element_access1.begin(), element_access1.end());
            ubw.insert(ubw.end(), element_access2.begin(), element_access2.end());
            this->lbw_ubw_range = tuple_size_t(bw_size, lbw.size());
        }

        template <class ProblemData>
        void PseudospectralSegment<ProblemData>::fill_lbg_ubg(std::vector<double> &lbg, std::vector<double> &ubg)
        {
            /*where lb/ub of this segment starts*/
            auto bg_size = lbg.size();
            std::vector<double> element_access1 = this->general_lbg.get_elements();
            std::vector<double> element_access2 = this->general_ubg.get_elements();

            lbg.insert(lbg.end(), element_access1.begin(), element_access1.end());
            ubg.insert(ubg.end(), element_access2.begin(), element_access2.end());

            this->lbg_ubg_range = tuple_size_t(bg_size, lbg.size());
        }

        template <class ProblemData>
        void PseudospectralSegment<ProblemData>::fill_w0(std::vector<double> &all_w0) const
        {
            std::vector<double> element_access1 = this->w0.get_elements();
            all_w0.insert(all_w0.end(), element_access1.begin(), element_access1.end());
        }

        template <class ProblemData>
        tuple_size_t PseudospectralSegment<ProblemData>::get_range_idx_decision_variables() const
        {
            return this->w_range;
        }

        template <class ProblemData>
        tuple_size_t PseudospectralSegment<ProblemData>::get_range_idx_constraint_expressions() const
        {
            return this->g_range;
        }

        template <class ProblemData>
        tuple_size_t PseudospectralSegment<ProblemData>::get_range_idx_constraint_bounds() const
        {
            return this->lbg_ubg_range;
        }

        template <class ProblemData>
        tuple_size_t PseudospectralSegment<ProblemData>::get_range_idx_decision_bounds() const
        {
            return this->lbw_ubw_range;
        }

    }
}