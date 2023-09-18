// Generated by rstantools.  Do not edit by hand.

#ifndef MODELS_HPP
#define MODELS_HPP
#define STAN__SERVICES__COMMAND_HPP
#include <rstan/rstaninc.hpp>
// Code generated by Stan version 2.21.0
#include <stan/model/model_header.hpp>
namespace model_multinom_namespace {
using std::istream;
using std::string;
using std::stringstream;
using std::vector;
using stan::io::dump;
using stan::math::lgamma;
using stan::model::prob_grad;
using namespace stan::math;
static int current_statement_begin__;
stan::io::program_reader prog_reader__() {
    stan::io::program_reader reader;
    reader.add_event(0, 0, "start", "model_multinom");
    reader.add_event(54, 52, "end", "model_multinom");
    return reader;
}
#include <stan_meta_header.hpp>
class model_multinom
  : public stan::model::model_base_crtp<model_multinom> {
private:
        int n_y;
        int N;
        int p;
        int n_grp;
        matrix_d X;
        matrix_d Y;
        std::vector<int> grp;
        int has_int;
        double prior_sigma;
        double prior_beta;
        double prior_int;
        vector_d ones_y;
        vector_d int_col;
public:
    model_multinom(stan::io::var_context& context__,
        std::ostream* pstream__ = 0)
        : model_base_crtp(0) {
        ctor_body(context__, 0, pstream__);
    }
    model_multinom(stan::io::var_context& context__,
        unsigned int random_seed__,
        std::ostream* pstream__ = 0)
        : model_base_crtp(0) {
        ctor_body(context__, random_seed__, pstream__);
    }
    void ctor_body(stan::io::var_context& context__,
                   unsigned int random_seed__,
                   std::ostream* pstream__) {
        typedef double local_scalar_t__;
        boost::ecuyer1988 base_rng__ =
          stan::services::util::create_rng(random_seed__, 0);
        (void) base_rng__;  // suppress unused var warning
        current_statement_begin__ = -1;
        static const char* function__ = "model_multinom_namespace::model_multinom";
        (void) function__;  // dummy to suppress unused var warning
        size_t pos__;
        (void) pos__;  // dummy to suppress unused var warning
        std::vector<int> vals_i__;
        std::vector<double> vals_r__;
        local_scalar_t__ DUMMY_VAR__(std::numeric_limits<double>::quiet_NaN());
        (void) DUMMY_VAR__;  // suppress unused var warning
        try {
            // initialize data block variables from context__
            current_statement_begin__ = 2;
            context__.validate_dims("data initialization", "n_y", "int", context__.to_vec());
            n_y = int(0);
            vals_i__ = context__.vals_i("n_y");
            pos__ = 0;
            n_y = vals_i__[pos__++];
            check_greater_or_equal(function__, "n_y", n_y, 0);
            current_statement_begin__ = 3;
            context__.validate_dims("data initialization", "N", "int", context__.to_vec());
            N = int(0);
            vals_i__ = context__.vals_i("N");
            pos__ = 0;
            N = vals_i__[pos__++];
            check_greater_or_equal(function__, "N", N, 0);
            current_statement_begin__ = 4;
            context__.validate_dims("data initialization", "p", "int", context__.to_vec());
            p = int(0);
            vals_i__ = context__.vals_i("p");
            pos__ = 0;
            p = vals_i__[pos__++];
            check_greater_or_equal(function__, "p", p, 0);
            current_statement_begin__ = 5;
            context__.validate_dims("data initialization", "n_grp", "int", context__.to_vec());
            n_grp = int(0);
            vals_i__ = context__.vals_i("n_grp");
            pos__ = 0;
            n_grp = vals_i__[pos__++];
            check_greater_or_equal(function__, "n_grp", n_grp, 0);
            current_statement_begin__ = 7;
            validate_non_negative_index("X", "N", N);
            validate_non_negative_index("X", "p", p);
            context__.validate_dims("data initialization", "X", "matrix_d", context__.to_vec(N,p));
            X = Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic>(N, p);
            vals_r__ = context__.vals_r("X");
            pos__ = 0;
            size_t X_j_2_max__ = p;
            size_t X_j_1_max__ = N;
            for (size_t j_2__ = 0; j_2__ < X_j_2_max__; ++j_2__) {
                for (size_t j_1__ = 0; j_1__ < X_j_1_max__; ++j_1__) {
                    X(j_1__, j_2__) = vals_r__[pos__++];
                }
            }
            current_statement_begin__ = 8;
            validate_non_negative_index("Y", "N", N);
            validate_non_negative_index("Y", "n_y", n_y);
            context__.validate_dims("data initialization", "Y", "matrix_d", context__.to_vec(N,n_y));
            Y = Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic>(N, n_y);
            vals_r__ = context__.vals_r("Y");
            pos__ = 0;
            size_t Y_j_2_max__ = n_y;
            size_t Y_j_1_max__ = N;
            for (size_t j_2__ = 0; j_2__ < Y_j_2_max__; ++j_2__) {
                for (size_t j_1__ = 0; j_1__ < Y_j_1_max__; ++j_1__) {
                    Y(j_1__, j_2__) = vals_r__[pos__++];
                }
            }
            current_statement_begin__ = 9;
            validate_non_negative_index("grp", "N", N);
            context__.validate_dims("data initialization", "grp", "int", context__.to_vec(N));
            grp = std::vector<int>(N, int(0));
            vals_i__ = context__.vals_i("grp");
            pos__ = 0;
            size_t grp_k_0_max__ = N;
            for (size_t k_0__ = 0; k_0__ < grp_k_0_max__; ++k_0__) {
                grp[k_0__] = vals_i__[pos__++];
            }
            size_t grp_i_0_max__ = N;
            for (size_t i_0__ = 0; i_0__ < grp_i_0_max__; ++i_0__) {
                check_greater_or_equal(function__, "grp[i_0__]", grp[i_0__], 1);
                check_less_or_equal(function__, "grp[i_0__]", grp[i_0__], n_grp);
            }
            current_statement_begin__ = 11;
            context__.validate_dims("data initialization", "has_int", "int", context__.to_vec());
            has_int = int(0);
            vals_i__ = context__.vals_i("has_int");
            pos__ = 0;
            has_int = vals_i__[pos__++];
            check_greater_or_equal(function__, "has_int", has_int, 0);
            check_less_or_equal(function__, "has_int", has_int, 1);
            current_statement_begin__ = 12;
            context__.validate_dims("data initialization", "prior_sigma", "double", context__.to_vec());
            prior_sigma = double(0);
            vals_r__ = context__.vals_r("prior_sigma");
            pos__ = 0;
            prior_sigma = vals_r__[pos__++];
            check_greater_or_equal(function__, "prior_sigma", prior_sigma, 0);
            current_statement_begin__ = 13;
            context__.validate_dims("data initialization", "prior_beta", "double", context__.to_vec());
            prior_beta = double(0);
            vals_r__ = context__.vals_r("prior_beta");
            pos__ = 0;
            prior_beta = vals_r__[pos__++];
            check_greater_or_equal(function__, "prior_beta", prior_beta, 0);
            current_statement_begin__ = 14;
            context__.validate_dims("data initialization", "prior_int", "double", context__.to_vec());
            prior_int = double(0);
            vals_r__ = context__.vals_r("prior_int");
            pos__ = 0;
            prior_int = vals_r__[pos__++];
            check_greater_or_equal(function__, "prior_int", prior_int, 0);
            // initialize transformed data variables
            current_statement_begin__ = 18;
            validate_non_negative_index("ones_y", "n_y", n_y);
            ones_y = Eigen::Matrix<double, Eigen::Dynamic, 1>(n_y);
            stan::math::fill(ones_y, DUMMY_VAR__);
            stan::math::assign(ones_y,rep_vector(1, n_y));
            current_statement_begin__ = 19;
            validate_non_negative_index("int_col", "N", N);
            int_col = Eigen::Matrix<double, Eigen::Dynamic, 1>(N);
            stan::math::fill(int_col, DUMMY_VAR__);
            stan::math::assign(int_col,rep_vector(has_int, N));
            // execute transformed data statements
            // validate transformed data
            // validate, set parameter ranges
            num_params_r__ = 0U;
            param_ranges_i__.clear();
            current_statement_begin__ = 23;
            validate_non_negative_index("intercept", "n_y", n_y);
            num_params_r__ += n_y;
            current_statement_begin__ = 24;
            validate_non_negative_index("beta", "p", p);
            validate_non_negative_index("beta", "n_y", n_y);
            num_params_r__ += (p * n_y);
            current_statement_begin__ = 25;
            validate_non_negative_index("u", "n_grp", n_grp);
            validate_non_negative_index("u", "n_y", n_y);
            num_params_r__ += (n_grp * n_y);
            current_statement_begin__ = 27;
            validate_non_negative_index("sigma_grp", "n_y", n_y);
            num_params_r__ += n_y;
            current_statement_begin__ = 28;
            validate_non_negative_index("L", "n_y", n_y);
            validate_non_negative_index("L", "n_y", n_y);
            num_params_r__ += ((n_y * (n_y - 1)) / 2);
        } catch (const std::exception& e) {
            stan::lang::rethrow_located(e, current_statement_begin__, prog_reader__());
            // Next line prevents compiler griping about no return
            throw std::runtime_error("*** IF YOU SEE THIS, PLEASE REPORT A BUG ***");
        }
    }
    ~model_multinom() { }
    void transform_inits(const stan::io::var_context& context__,
                         std::vector<int>& params_i__,
                         std::vector<double>& params_r__,
                         std::ostream* pstream__) const {
        typedef double local_scalar_t__;
        stan::io::writer<double> writer__(params_r__, params_i__);
        size_t pos__;
        (void) pos__; // dummy call to supress warning
        std::vector<double> vals_r__;
        std::vector<int> vals_i__;
        current_statement_begin__ = 23;
        if (!(context__.contains_r("intercept")))
            stan::lang::rethrow_located(std::runtime_error(std::string("Variable intercept missing")), current_statement_begin__, prog_reader__());
        vals_r__ = context__.vals_r("intercept");
        pos__ = 0U;
        validate_non_negative_index("intercept", "n_y", n_y);
        context__.validate_dims("parameter initialization", "intercept", "row_vector_d", context__.to_vec(n_y));
        Eigen::Matrix<double, 1, Eigen::Dynamic> intercept(n_y);
        size_t intercept_j_1_max__ = n_y;
        for (size_t j_1__ = 0; j_1__ < intercept_j_1_max__; ++j_1__) {
            intercept(j_1__) = vals_r__[pos__++];
        }
        try {
            writer__.row_vector_unconstrain(intercept);
        } catch (const std::exception& e) {
            stan::lang::rethrow_located(std::runtime_error(std::string("Error transforming variable intercept: ") + e.what()), current_statement_begin__, prog_reader__());
        }
        current_statement_begin__ = 24;
        if (!(context__.contains_r("beta")))
            stan::lang::rethrow_located(std::runtime_error(std::string("Variable beta missing")), current_statement_begin__, prog_reader__());
        vals_r__ = context__.vals_r("beta");
        pos__ = 0U;
        validate_non_negative_index("beta", "p", p);
        validate_non_negative_index("beta", "n_y", n_y);
        context__.validate_dims("parameter initialization", "beta", "matrix_d", context__.to_vec(p,n_y));
        Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> beta(p, n_y);
        size_t beta_j_2_max__ = n_y;
        size_t beta_j_1_max__ = p;
        for (size_t j_2__ = 0; j_2__ < beta_j_2_max__; ++j_2__) {
            for (size_t j_1__ = 0; j_1__ < beta_j_1_max__; ++j_1__) {
                beta(j_1__, j_2__) = vals_r__[pos__++];
            }
        }
        try {
            writer__.matrix_unconstrain(beta);
        } catch (const std::exception& e) {
            stan::lang::rethrow_located(std::runtime_error(std::string("Error transforming variable beta: ") + e.what()), current_statement_begin__, prog_reader__());
        }
        current_statement_begin__ = 25;
        if (!(context__.contains_r("u")))
            stan::lang::rethrow_located(std::runtime_error(std::string("Variable u missing")), current_statement_begin__, prog_reader__());
        vals_r__ = context__.vals_r("u");
        pos__ = 0U;
        validate_non_negative_index("u", "n_grp", n_grp);
        validate_non_negative_index("u", "n_y", n_y);
        context__.validate_dims("parameter initialization", "u", "matrix_d", context__.to_vec(n_grp,n_y));
        Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> u(n_grp, n_y);
        size_t u_j_2_max__ = n_y;
        size_t u_j_1_max__ = n_grp;
        for (size_t j_2__ = 0; j_2__ < u_j_2_max__; ++j_2__) {
            for (size_t j_1__ = 0; j_1__ < u_j_1_max__; ++j_1__) {
                u(j_1__, j_2__) = vals_r__[pos__++];
            }
        }
        try {
            writer__.matrix_unconstrain(u);
        } catch (const std::exception& e) {
            stan::lang::rethrow_located(std::runtime_error(std::string("Error transforming variable u: ") + e.what()), current_statement_begin__, prog_reader__());
        }
        current_statement_begin__ = 27;
        if (!(context__.contains_r("sigma_grp")))
            stan::lang::rethrow_located(std::runtime_error(std::string("Variable sigma_grp missing")), current_statement_begin__, prog_reader__());
        vals_r__ = context__.vals_r("sigma_grp");
        pos__ = 0U;
        validate_non_negative_index("sigma_grp", "n_y", n_y);
        context__.validate_dims("parameter initialization", "sigma_grp", "vector_d", context__.to_vec(n_y));
        Eigen::Matrix<double, Eigen::Dynamic, 1> sigma_grp(n_y);
        size_t sigma_grp_j_1_max__ = n_y;
        for (size_t j_1__ = 0; j_1__ < sigma_grp_j_1_max__; ++j_1__) {
            sigma_grp(j_1__) = vals_r__[pos__++];
        }
        try {
            writer__.vector_lb_unconstrain(0, sigma_grp);
        } catch (const std::exception& e) {
            stan::lang::rethrow_located(std::runtime_error(std::string("Error transforming variable sigma_grp: ") + e.what()), current_statement_begin__, prog_reader__());
        }
        current_statement_begin__ = 28;
        if (!(context__.contains_r("L")))
            stan::lang::rethrow_located(std::runtime_error(std::string("Variable L missing")), current_statement_begin__, prog_reader__());
        vals_r__ = context__.vals_r("L");
        pos__ = 0U;
        validate_non_negative_index("L", "n_y", n_y);
        validate_non_negative_index("L", "n_y", n_y);
        context__.validate_dims("parameter initialization", "L", "matrix_d", context__.to_vec(n_y,n_y));
        Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> L(n_y, n_y);
        size_t L_j_2_max__ = n_y;
        size_t L_j_1_max__ = n_y;
        for (size_t j_2__ = 0; j_2__ < L_j_2_max__; ++j_2__) {
            for (size_t j_1__ = 0; j_1__ < L_j_1_max__; ++j_1__) {
                L(j_1__, j_2__) = vals_r__[pos__++];
            }
        }
        try {
            writer__.cholesky_factor_corr_unconstrain(L);
        } catch (const std::exception& e) {
            stan::lang::rethrow_located(std::runtime_error(std::string("Error transforming variable L: ") + e.what()), current_statement_begin__, prog_reader__());
        }
        params_r__ = writer__.data_r();
        params_i__ = writer__.data_i();
    }
    void transform_inits(const stan::io::var_context& context,
                         Eigen::Matrix<double, Eigen::Dynamic, 1>& params_r,
                         std::ostream* pstream__) const {
      std::vector<double> params_r_vec;
      std::vector<int> params_i_vec;
      transform_inits(context, params_i_vec, params_r_vec, pstream__);
      params_r.resize(params_r_vec.size());
      for (int i = 0; i < params_r.size(); ++i)
        params_r(i) = params_r_vec[i];
    }
    template <bool propto__, bool jacobian__, typename T__>
    T__ log_prob(std::vector<T__>& params_r__,
                 std::vector<int>& params_i__,
                 std::ostream* pstream__ = 0) const {
        typedef T__ local_scalar_t__;
        local_scalar_t__ DUMMY_VAR__(std::numeric_limits<double>::quiet_NaN());
        (void) DUMMY_VAR__;  // dummy to suppress unused var warning
        T__ lp__(0.0);
        stan::math::accumulator<T__> lp_accum__;
        try {
            stan::io::reader<local_scalar_t__> in__(params_r__, params_i__);
            // model parameters
            current_statement_begin__ = 23;
            Eigen::Matrix<local_scalar_t__, 1, Eigen::Dynamic> intercept;
            (void) intercept;  // dummy to suppress unused var warning
            if (jacobian__)
                intercept = in__.row_vector_constrain(n_y, lp__);
            else
                intercept = in__.row_vector_constrain(n_y);
            current_statement_begin__ = 24;
            Eigen::Matrix<local_scalar_t__, Eigen::Dynamic, Eigen::Dynamic> beta;
            (void) beta;  // dummy to suppress unused var warning
            if (jacobian__)
                beta = in__.matrix_constrain(p, n_y, lp__);
            else
                beta = in__.matrix_constrain(p, n_y);
            current_statement_begin__ = 25;
            Eigen::Matrix<local_scalar_t__, Eigen::Dynamic, Eigen::Dynamic> u;
            (void) u;  // dummy to suppress unused var warning
            if (jacobian__)
                u = in__.matrix_constrain(n_grp, n_y, lp__);
            else
                u = in__.matrix_constrain(n_grp, n_y);
            current_statement_begin__ = 27;
            Eigen::Matrix<local_scalar_t__, Eigen::Dynamic, 1> sigma_grp;
            (void) sigma_grp;  // dummy to suppress unused var warning
            if (jacobian__)
                sigma_grp = in__.vector_lb_constrain(0, n_y, lp__);
            else
                sigma_grp = in__.vector_lb_constrain(0, n_y);
            current_statement_begin__ = 28;
            Eigen::Matrix<local_scalar_t__, Eigen::Dynamic, Eigen::Dynamic> L;
            (void) L;  // dummy to suppress unused var warning
            if (jacobian__)
                L = in__.cholesky_factor_corr_constrain(n_y, lp__);
            else
                L = in__.cholesky_factor_corr_constrain(n_y);
            // transformed parameters
            current_statement_begin__ = 32;
            validate_non_negative_index("lsft", "N", N);
            validate_non_negative_index("lsft", "n_y", n_y);
            Eigen::Matrix<local_scalar_t__, Eigen::Dynamic, Eigen::Dynamic> lsft(N, n_y);
            stan::math::initialize(lsft, DUMMY_VAR__);
            stan::math::fill(lsft, DUMMY_VAR__);
            // transformed parameters block statements
            {
            current_statement_begin__ = 34;
            validate_non_negative_index("linpred", "N", N);
            validate_non_negative_index("linpred", "n_y", n_y);
            Eigen::Matrix<local_scalar_t__, Eigen::Dynamic, Eigen::Dynamic> linpred(N, n_y);
            stan::math::initialize(linpred, DUMMY_VAR__);
            stan::math::fill(linpred, DUMMY_VAR__);
            current_statement_begin__ = 35;
            validate_non_negative_index("Sigma", "n_y", n_y);
            validate_non_negative_index("Sigma", "n_y", n_y);
            Eigen::Matrix<local_scalar_t__, Eigen::Dynamic, Eigen::Dynamic> Sigma(n_y, n_y);
            stan::math::initialize(Sigma, DUMMY_VAR__);
            stan::math::fill(Sigma, DUMMY_VAR__);
            stan::math::assign(Sigma,diag_pre_multiply(sigma_grp, L));
            current_statement_begin__ = 36;
            stan::math::assign(linpred, add(add(multiply(int_col, intercept), multiply(X, beta)), transpose(multiply(Sigma, transpose(stan::model::rvalue(u, stan::model::cons_list(stan::model::index_multi(grp), stan::model::nil_index_list()), "u"))))));
            current_statement_begin__ = 39;
            stan::math::assign(lsft, subtract(linpred, rep_matrix(stan::math::log(multiply(stan::math::exp(linpred), ones_y)), n_y)));
            }
            // validate transformed parameters
            const char* function__ = "validate transformed params";
            (void) function__;  // dummy to suppress unused var warning
            current_statement_begin__ = 32;
            size_t lsft_j_1_max__ = N;
            size_t lsft_j_2_max__ = n_y;
            for (size_t j_1__ = 0; j_1__ < lsft_j_1_max__; ++j_1__) {
                for (size_t j_2__ = 0; j_2__ < lsft_j_2_max__; ++j_2__) {
                    if (stan::math::is_uninitialized(lsft(j_1__, j_2__))) {
                        std::stringstream msg__;
                        msg__ << "Undefined transformed parameter: lsft" << "(" << j_1__ << ", " << j_2__ << ")";
                        stan::lang::rethrow_located(std::runtime_error(std::string("Error initializing variable lsft: ") + msg__.str()), current_statement_begin__, prog_reader__());
                    }
                }
            }
            // model body
            current_statement_begin__ = 44;
            lp_accum__.add(sum(elt_multiply(Y, lsft)));
            current_statement_begin__ = 46;
            lp_accum__.add(normal_log<propto__>(intercept, 0, prior_int));
            current_statement_begin__ = 47;
            lp_accum__.add(normal_log<propto__>(to_vector(beta), 0, prior_beta));
            current_statement_begin__ = 48;
            lp_accum__.add(std_normal_log<propto__>(to_vector(u)));
            current_statement_begin__ = 50;
            lp_accum__.add(gamma_log<propto__>(sigma_grp, 2.0, (2.0 / prior_sigma)));
            current_statement_begin__ = 51;
            lp_accum__.add(lkj_corr_cholesky_log<propto__>(L, 2.0));
        } catch (const std::exception& e) {
            stan::lang::rethrow_located(e, current_statement_begin__, prog_reader__());
            // Next line prevents compiler griping about no return
            throw std::runtime_error("*** IF YOU SEE THIS, PLEASE REPORT A BUG ***");
        }
        lp_accum__.add(lp__);
        return lp_accum__.sum();
    } // log_prob()
    template <bool propto, bool jacobian, typename T_>
    T_ log_prob(Eigen::Matrix<T_,Eigen::Dynamic,1>& params_r,
               std::ostream* pstream = 0) const {
      std::vector<T_> vec_params_r;
      vec_params_r.reserve(params_r.size());
      for (int i = 0; i < params_r.size(); ++i)
        vec_params_r.push_back(params_r(i));
      std::vector<int> vec_params_i;
      return log_prob<propto,jacobian,T_>(vec_params_r, vec_params_i, pstream);
    }
    void get_param_names(std::vector<std::string>& names__) const {
        names__.resize(0);
        names__.push_back("intercept");
        names__.push_back("beta");
        names__.push_back("u");
        names__.push_back("sigma_grp");
        names__.push_back("L");
        names__.push_back("lsft");
    }
    void get_dims(std::vector<std::vector<size_t> >& dimss__) const {
        dimss__.resize(0);
        std::vector<size_t> dims__;
        dims__.resize(0);
        dims__.push_back(n_y);
        dimss__.push_back(dims__);
        dims__.resize(0);
        dims__.push_back(p);
        dims__.push_back(n_y);
        dimss__.push_back(dims__);
        dims__.resize(0);
        dims__.push_back(n_grp);
        dims__.push_back(n_y);
        dimss__.push_back(dims__);
        dims__.resize(0);
        dims__.push_back(n_y);
        dimss__.push_back(dims__);
        dims__.resize(0);
        dims__.push_back(n_y);
        dims__.push_back(n_y);
        dimss__.push_back(dims__);
        dims__.resize(0);
        dims__.push_back(N);
        dims__.push_back(n_y);
        dimss__.push_back(dims__);
    }
    template <typename RNG>
    void write_array(RNG& base_rng__,
                     std::vector<double>& params_r__,
                     std::vector<int>& params_i__,
                     std::vector<double>& vars__,
                     bool include_tparams__ = true,
                     bool include_gqs__ = true,
                     std::ostream* pstream__ = 0) const {
        typedef double local_scalar_t__;
        vars__.resize(0);
        stan::io::reader<local_scalar_t__> in__(params_r__, params_i__);
        static const char* function__ = "model_multinom_namespace::write_array";
        (void) function__;  // dummy to suppress unused var warning
        // read-transform, write parameters
        Eigen::Matrix<double, 1, Eigen::Dynamic> intercept = in__.row_vector_constrain(n_y);
        size_t intercept_j_1_max__ = n_y;
        for (size_t j_1__ = 0; j_1__ < intercept_j_1_max__; ++j_1__) {
            vars__.push_back(intercept(j_1__));
        }
        Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> beta = in__.matrix_constrain(p, n_y);
        size_t beta_j_2_max__ = n_y;
        size_t beta_j_1_max__ = p;
        for (size_t j_2__ = 0; j_2__ < beta_j_2_max__; ++j_2__) {
            for (size_t j_1__ = 0; j_1__ < beta_j_1_max__; ++j_1__) {
                vars__.push_back(beta(j_1__, j_2__));
            }
        }
        Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> u = in__.matrix_constrain(n_grp, n_y);
        size_t u_j_2_max__ = n_y;
        size_t u_j_1_max__ = n_grp;
        for (size_t j_2__ = 0; j_2__ < u_j_2_max__; ++j_2__) {
            for (size_t j_1__ = 0; j_1__ < u_j_1_max__; ++j_1__) {
                vars__.push_back(u(j_1__, j_2__));
            }
        }
        Eigen::Matrix<double, Eigen::Dynamic, 1> sigma_grp = in__.vector_lb_constrain(0, n_y);
        size_t sigma_grp_j_1_max__ = n_y;
        for (size_t j_1__ = 0; j_1__ < sigma_grp_j_1_max__; ++j_1__) {
            vars__.push_back(sigma_grp(j_1__));
        }
        Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> L = in__.cholesky_factor_corr_constrain(n_y);
        size_t L_j_2_max__ = n_y;
        size_t L_j_1_max__ = n_y;
        for (size_t j_2__ = 0; j_2__ < L_j_2_max__; ++j_2__) {
            for (size_t j_1__ = 0; j_1__ < L_j_1_max__; ++j_1__) {
                vars__.push_back(L(j_1__, j_2__));
            }
        }
        double lp__ = 0.0;
        (void) lp__;  // dummy to suppress unused var warning
        stan::math::accumulator<double> lp_accum__;
        local_scalar_t__ DUMMY_VAR__(std::numeric_limits<double>::quiet_NaN());
        (void) DUMMY_VAR__;  // suppress unused var warning
        if (!include_tparams__ && !include_gqs__) return;
        try {
            // declare and define transformed parameters
            current_statement_begin__ = 32;
            validate_non_negative_index("lsft", "N", N);
            validate_non_negative_index("lsft", "n_y", n_y);
            Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> lsft(N, n_y);
            stan::math::initialize(lsft, DUMMY_VAR__);
            stan::math::fill(lsft, DUMMY_VAR__);
            // do transformed parameters statements
            {
            current_statement_begin__ = 34;
            validate_non_negative_index("linpred", "N", N);
            validate_non_negative_index("linpred", "n_y", n_y);
            Eigen::Matrix<local_scalar_t__, Eigen::Dynamic, Eigen::Dynamic> linpred(N, n_y);
            stan::math::initialize(linpred, DUMMY_VAR__);
            stan::math::fill(linpred, DUMMY_VAR__);
            current_statement_begin__ = 35;
            validate_non_negative_index("Sigma", "n_y", n_y);
            validate_non_negative_index("Sigma", "n_y", n_y);
            Eigen::Matrix<local_scalar_t__, Eigen::Dynamic, Eigen::Dynamic> Sigma(n_y, n_y);
            stan::math::initialize(Sigma, DUMMY_VAR__);
            stan::math::fill(Sigma, DUMMY_VAR__);
            stan::math::assign(Sigma,diag_pre_multiply(sigma_grp, L));
            current_statement_begin__ = 36;
            stan::math::assign(linpred, add(add(multiply(int_col, intercept), multiply(X, beta)), transpose(multiply(Sigma, transpose(stan::model::rvalue(u, stan::model::cons_list(stan::model::index_multi(grp), stan::model::nil_index_list()), "u"))))));
            current_statement_begin__ = 39;
            stan::math::assign(lsft, subtract(linpred, rep_matrix(stan::math::log(multiply(stan::math::exp(linpred), ones_y)), n_y)));
            }
            if (!include_gqs__ && !include_tparams__) return;
            // validate transformed parameters
            const char* function__ = "validate transformed params";
            (void) function__;  // dummy to suppress unused var warning
            // write transformed parameters
            if (include_tparams__) {
                size_t lsft_j_2_max__ = n_y;
                size_t lsft_j_1_max__ = N;
                for (size_t j_2__ = 0; j_2__ < lsft_j_2_max__; ++j_2__) {
                    for (size_t j_1__ = 0; j_1__ < lsft_j_1_max__; ++j_1__) {
                        vars__.push_back(lsft(j_1__, j_2__));
                    }
                }
            }
            if (!include_gqs__) return;
        } catch (const std::exception& e) {
            stan::lang::rethrow_located(e, current_statement_begin__, prog_reader__());
            // Next line prevents compiler griping about no return
            throw std::runtime_error("*** IF YOU SEE THIS, PLEASE REPORT A BUG ***");
        }
    }
    template <typename RNG>
    void write_array(RNG& base_rng,
                     Eigen::Matrix<double,Eigen::Dynamic,1>& params_r,
                     Eigen::Matrix<double,Eigen::Dynamic,1>& vars,
                     bool include_tparams = true,
                     bool include_gqs = true,
                     std::ostream* pstream = 0) const {
      std::vector<double> params_r_vec(params_r.size());
      for (int i = 0; i < params_r.size(); ++i)
        params_r_vec[i] = params_r(i);
      std::vector<double> vars_vec;
      std::vector<int> params_i_vec;
      write_array(base_rng, params_r_vec, params_i_vec, vars_vec, include_tparams, include_gqs, pstream);
      vars.resize(vars_vec.size());
      for (int i = 0; i < vars.size(); ++i)
        vars(i) = vars_vec[i];
    }
    std::string model_name() const {
        return "model_multinom";
    }
    void constrained_param_names(std::vector<std::string>& param_names__,
                                 bool include_tparams__ = true,
                                 bool include_gqs__ = true) const {
        std::stringstream param_name_stream__;
        size_t intercept_j_1_max__ = n_y;
        for (size_t j_1__ = 0; j_1__ < intercept_j_1_max__; ++j_1__) {
            param_name_stream__.str(std::string());
            param_name_stream__ << "intercept" << '.' << j_1__ + 1;
            param_names__.push_back(param_name_stream__.str());
        }
        size_t beta_j_2_max__ = n_y;
        size_t beta_j_1_max__ = p;
        for (size_t j_2__ = 0; j_2__ < beta_j_2_max__; ++j_2__) {
            for (size_t j_1__ = 0; j_1__ < beta_j_1_max__; ++j_1__) {
                param_name_stream__.str(std::string());
                param_name_stream__ << "beta" << '.' << j_1__ + 1 << '.' << j_2__ + 1;
                param_names__.push_back(param_name_stream__.str());
            }
        }
        size_t u_j_2_max__ = n_y;
        size_t u_j_1_max__ = n_grp;
        for (size_t j_2__ = 0; j_2__ < u_j_2_max__; ++j_2__) {
            for (size_t j_1__ = 0; j_1__ < u_j_1_max__; ++j_1__) {
                param_name_stream__.str(std::string());
                param_name_stream__ << "u" << '.' << j_1__ + 1 << '.' << j_2__ + 1;
                param_names__.push_back(param_name_stream__.str());
            }
        }
        size_t sigma_grp_j_1_max__ = n_y;
        for (size_t j_1__ = 0; j_1__ < sigma_grp_j_1_max__; ++j_1__) {
            param_name_stream__.str(std::string());
            param_name_stream__ << "sigma_grp" << '.' << j_1__ + 1;
            param_names__.push_back(param_name_stream__.str());
        }
        size_t L_j_2_max__ = n_y;
        size_t L_j_1_max__ = n_y;
        for (size_t j_2__ = 0; j_2__ < L_j_2_max__; ++j_2__) {
            for (size_t j_1__ = 0; j_1__ < L_j_1_max__; ++j_1__) {
                param_name_stream__.str(std::string());
                param_name_stream__ << "L" << '.' << j_1__ + 1 << '.' << j_2__ + 1;
                param_names__.push_back(param_name_stream__.str());
            }
        }
        if (!include_gqs__ && !include_tparams__) return;
        if (include_tparams__) {
            size_t lsft_j_2_max__ = n_y;
            size_t lsft_j_1_max__ = N;
            for (size_t j_2__ = 0; j_2__ < lsft_j_2_max__; ++j_2__) {
                for (size_t j_1__ = 0; j_1__ < lsft_j_1_max__; ++j_1__) {
                    param_name_stream__.str(std::string());
                    param_name_stream__ << "lsft" << '.' << j_1__ + 1 << '.' << j_2__ + 1;
                    param_names__.push_back(param_name_stream__.str());
                }
            }
        }
        if (!include_gqs__) return;
    }
    void unconstrained_param_names(std::vector<std::string>& param_names__,
                                   bool include_tparams__ = true,
                                   bool include_gqs__ = true) const {
        std::stringstream param_name_stream__;
        size_t intercept_j_1_max__ = n_y;
        for (size_t j_1__ = 0; j_1__ < intercept_j_1_max__; ++j_1__) {
            param_name_stream__.str(std::string());
            param_name_stream__ << "intercept" << '.' << j_1__ + 1;
            param_names__.push_back(param_name_stream__.str());
        }
        size_t beta_j_2_max__ = n_y;
        size_t beta_j_1_max__ = p;
        for (size_t j_2__ = 0; j_2__ < beta_j_2_max__; ++j_2__) {
            for (size_t j_1__ = 0; j_1__ < beta_j_1_max__; ++j_1__) {
                param_name_stream__.str(std::string());
                param_name_stream__ << "beta" << '.' << j_1__ + 1 << '.' << j_2__ + 1;
                param_names__.push_back(param_name_stream__.str());
            }
        }
        size_t u_j_2_max__ = n_y;
        size_t u_j_1_max__ = n_grp;
        for (size_t j_2__ = 0; j_2__ < u_j_2_max__; ++j_2__) {
            for (size_t j_1__ = 0; j_1__ < u_j_1_max__; ++j_1__) {
                param_name_stream__.str(std::string());
                param_name_stream__ << "u" << '.' << j_1__ + 1 << '.' << j_2__ + 1;
                param_names__.push_back(param_name_stream__.str());
            }
        }
        size_t sigma_grp_j_1_max__ = n_y;
        for (size_t j_1__ = 0; j_1__ < sigma_grp_j_1_max__; ++j_1__) {
            param_name_stream__.str(std::string());
            param_name_stream__ << "sigma_grp" << '.' << j_1__ + 1;
            param_names__.push_back(param_name_stream__.str());
        }
        size_t L_j_1_max__ = ((n_y * (n_y - 1)) / 2);
        for (size_t j_1__ = 0; j_1__ < L_j_1_max__; ++j_1__) {
            param_name_stream__.str(std::string());
            param_name_stream__ << "L" << '.' << j_1__ + 1;
            param_names__.push_back(param_name_stream__.str());
        }
        if (!include_gqs__ && !include_tparams__) return;
        if (include_tparams__) {
            size_t lsft_j_2_max__ = n_y;
            size_t lsft_j_1_max__ = N;
            for (size_t j_2__ = 0; j_2__ < lsft_j_2_max__; ++j_2__) {
                for (size_t j_1__ = 0; j_1__ < lsft_j_1_max__; ++j_1__) {
                    param_name_stream__.str(std::string());
                    param_name_stream__ << "lsft" << '.' << j_1__ + 1 << '.' << j_2__ + 1;
                    param_names__.push_back(param_name_stream__.str());
                }
            }
        }
        if (!include_gqs__) return;
    }
}; // model
}  // namespace
typedef model_multinom_namespace::model_multinom stan_model;
#ifndef USING_R
stan::model::model_base& new_model(
        stan::io::var_context& data_context,
        unsigned int seed,
        std::ostream* msg_stream) {
  stan_model* m = new stan_model(data_context, seed, msg_stream);
  return *m;
}
#endif
#endif
