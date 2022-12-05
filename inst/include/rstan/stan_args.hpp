/*
 * ORIGINALLY FROM rstan PACKAGE, `inst/include/rstan` DIRECTORY
 * STRIPPED DOWN TO CODE NEEDED FOR OPTIMIZATION ONLY
 * (c) STAN DEVELOPMENT TEAM 2022
 */

#ifndef RSTAN__STAN_ARGS_HPP
#define RSTAN__STAN_ARGS_HPP


#include <Rcpp.h>
// #include <R.h>
// #include <Rinternals.h>

#include <algorithm>
#include <rstan/io/r_ostream.hpp>
#include <stan/version.hpp>
#include <boost/lexical_cast.hpp>

namespace rstan {

  namespace {
    /*
     * Get an element of Rcpp::List by name. If not found, set it to a
     * default value.
     * @param lst The list to look for elements
     * @param n The name of an element of interest
     * @param t Where to save the element
     * @param v0 The default value if not found in the list
     */
    template <class T>
    bool get_rlist_element(const Rcpp::List& lst, const char* n, T& t, const T& v0) {
      bool b = lst.containsElementNamed(n);
      if (b)  t = Rcpp::as<T>(const_cast<Rcpp::List&>(lst)[n]);
      else  t = T(v0);
      return b;
    }

    template <class T>
    bool get_rlist_element(const Rcpp::List& lst, const char* n, T& t) {
      bool b = lst.containsElementNamed(n);
      if (b) t = Rcpp::as<T>(const_cast<Rcpp::List&>(lst)[n]);
      return b;
    }

    template <>
    bool get_rlist_element(const Rcpp::List& lst, const char* n, SEXP& t) {
      bool b = lst.containsElementNamed(n);
      if (b) t = const_cast<Rcpp::List&>(lst)[n];
      return b;
    }

    inline unsigned int sexp2seed(SEXP seed) {
      if (TYPEOF(seed) == STRSXP)
        return boost::lexical_cast<unsigned int>(Rcpp::as<std::string>(seed));
      return Rcpp::as<unsigned int>(seed);
    }

    void write_comment(std::ostream& o) {
      o << "#" << std::endl;
    }

    template <typename M>
    void write_comment(std::ostream& o, const M& msg) {
      o << "# " << msg << std::endl;
    }

    template <typename K, typename V>
    void write_comment_property(std::ostream& o, const K& key, const V& val) {
      o << "# " << key << "=" << val << std::endl;
    }

    /**
     * Find the index of an element in a vector.
     * @param v the vector in which an element are searched.
     * @param e the element that we are looking for.
     * @return If e is in v, return the index (0 to size - 1);
     *  otherwise, return the size.
     */

    template <class T, class T2>
    size_t find_index(const std::vector<T>& v, const T2& e) {
      return std::distance(v.begin(), std::find(v.begin(), v.end(), T(e)));
    }
  }

  enum optim_algo_t { Newton = 1, BFGS = 3, LBFGS = 4};
  enum stan_args_method_t { OPTIM = 2 };

  /**
   *
   */
  class stan_args {
  private:
    unsigned int random_seed;
    unsigned int chain_id;
    std::string init;
    SEXP init_list;
    double init_radius;
    // FIXME(syclik): remove `enable_random_init`
    bool enable_random_init; // enable randomly partially specifying inits
    std::string sample_file; // the file for outputting the samples
    bool append_samples;
    bool sample_file_flag; // true: write out to a file; false, do not
    stan_args_method_t method;
    std::string diagnostic_file;
    bool diagnostic_file_flag;
    union {
      struct {
        int iter; // default to 2000
        int refresh; // default to 100
        optim_algo_t algorithm; // Newton, (L)BFGS
        bool save_iterations; // default to false
        double init_alpha; // default to 0.001, for (L)BFGS
        double tol_obj; // default to 1e-12, for (L)BFGS
        double tol_grad; // default to 1e-8, for (L)BFGS
        double tol_param; // default to 1e-8, for (L)BFGS
        double tol_rel_obj; // default to 1e4, for (L)BFGS
        double tol_rel_grad; // default to 1e7, for (L)BFGS
        int history_size; // default to 5, for LBFGS only
      } optim;
    } ctrl;

  private:
    void validate_args() {
      if (init_radius < 0) {
        std::stringstream msg;
        msg << "Invalid value for parameter init_r (found "
            << init_radius << "; require >= 0).";
        throw std::invalid_argument(msg.str());
      }
      switch (method) {
        case OPTIM:
          if (ctrl.optim.init_alpha < 0) {
            std::stringstream msg;
            msg << "Invalid adaptation parameter (found init_alpha="
                << ctrl.optim.init_alpha << "; require init_alpha > 0).";
            throw std::invalid_argument(msg.str());
          }
          break;
      }
    }

  public:
    stan_args(const Rcpp::List& in) : init_list(R_NilValue) {

      std::string t_str;
      SEXP t_sexp;
      bool b;
      get_rlist_element(in, "chain_id", chain_id, static_cast<unsigned int>(1));
      get_rlist_element(in, "append_samples", append_samples, false);
      b = get_rlist_element(in, "method", t_str);
      if (!b) method = OPTIM;
      else {
        if ("optim" == t_str)  method = OPTIM;
        else method = OPTIM;
      }

      sample_file_flag = get_rlist_element(in, "sample_file", sample_file);
      diagnostic_file_flag = get_rlist_element(in, "diagnostic_file", diagnostic_file);
      b = get_rlist_element(in, "seed", t_sexp);
      if (b) random_seed = sexp2seed(t_sexp);
      else random_seed = std::time(0);

      int calculated_thin;
      get_rlist_element(in, "control", t_sexp, R_NilValue);
      Rcpp::List ctrl_lst(t_sexp);

      switch (method) {
        case OPTIM:
          get_rlist_element(in, "iter", ctrl.optim.iter, 2000);
          if (get_rlist_element(in, "algorithm", t_str)) {
            if ("BFGS" == t_str)  ctrl.optim.algorithm = BFGS;
            else if ("Newton" == t_str)  ctrl.optim.algorithm = Newton;
            else if ("LBFGS" == t_str)  ctrl.optim.algorithm = LBFGS;
            else {
              std::stringstream msg;
              msg << "Invalid value for parameter algorithm (found "
                  << t_str << "; require (L)BFGS or Newton).";
              throw std::invalid_argument(msg.str());
            }
          } else {
            ctrl.optim.algorithm = LBFGS;
          }

          if (!get_rlist_element(in, "refresh", ctrl.optim.refresh)) {
            ctrl.optim.refresh = ctrl.optim.iter / 100;
            if (ctrl.optim.refresh < 1) ctrl.optim.refresh = 1;
          }

          get_rlist_element(in, "init_alpha", ctrl.optim.init_alpha, 0.001);
          get_rlist_element(in, "tol_obj", ctrl.optim.tol_obj, 1e-12);
          get_rlist_element(in, "tol_grad", ctrl.optim.tol_grad, 1e-8);
          get_rlist_element(in, "tol_param", ctrl.optim.tol_param, 1e-8);
          get_rlist_element(in, "tol_rel_obj", ctrl.optim.tol_rel_obj, 1e4);
          get_rlist_element(in, "tol_rel_grad", ctrl.optim.tol_rel_grad, 1e7);
          get_rlist_element(in, "save_iterations", ctrl.optim.save_iterations, true);
          get_rlist_element(in, "history_size", ctrl.optim.history_size, static_cast<int>(5));
          break;
      }

      if (get_rlist_element(in, "init", t_sexp)) {
        switch (TYPEOF(t_sexp)) {
          case STRSXP: init = Rcpp::as<std::string>(t_sexp); break;
          case VECSXP: init = "user"; init_list = t_sexp; break;
          default: init = "random";
        }
      } else {
        init = "random";
      }
      get_rlist_element(in, "init_r", init_radius, 2.0);
      if (0 >= init_radius)  init = "0";
      if (init == "0") init_radius = 0;
      get_rlist_element(in, "enable_random_init", enable_random_init, true);
      validate_args();
    }

    /**
     * return all the arguments used as an R list
     * @return An R list containing all the arguments for a chain.
     */
    SEXP stan_args_to_rlist() const {
      std::map<std::string, SEXP> args;
      std::map<std::string, SEXP> ctrl_args;
      std::stringstream ss;
      ss << random_seed;
      args["random_seed"] = Rcpp::wrap(ss.str());
      args["chain_id"] = Rcpp::wrap(chain_id);
      args["init"] = Rcpp::wrap(init);
      args["init_list"] = init_list;
      args["init_radius"] = Rcpp::wrap(init_radius);
      args["enable_random_init"] = Rcpp::wrap(enable_random_init);
      args["append_samples"] = Rcpp::wrap(append_samples);
      if (sample_file_flag)
        args["sample_file"] = Rcpp::wrap(sample_file);
      if (diagnostic_file_flag)
        args["diagnostic_file_flag"] = Rcpp::wrap(diagnostic_file);

      std::string sampler_t;
      switch (method) {
        case OPTIM:
          args["method"] = Rcpp::wrap("optim");
          args["iter"] = Rcpp::wrap(ctrl.optim.iter);
          args["refresh"] = Rcpp::wrap(ctrl.optim.refresh);
          args["save_iterations"] = Rcpp::wrap(ctrl.optim.save_iterations);
          switch (ctrl.optim.algorithm) {
            case Newton: args["algorithm"] = Rcpp::wrap("Newton"); break;
            case LBFGS: args["algorithm"] = Rcpp::wrap("LBFGS");
                        args["init_alpha"] = Rcpp::wrap(ctrl.optim.init_alpha);
                        args["tol_param"] = Rcpp::wrap(ctrl.optim.tol_param);
                        args["tol_obj"] = Rcpp::wrap(ctrl.optim.tol_obj);
                        args["tol_grad"] = Rcpp::wrap(ctrl.optim.tol_grad);
                        args["tol_rel_obj"] = Rcpp::wrap(ctrl.optim.tol_rel_obj);
                        args["tol_rel_grad"] = Rcpp::wrap(ctrl.optim.tol_rel_grad);
                        args["history_size"] = Rcpp::wrap(ctrl.optim.history_size);
                        break;
            case BFGS: args["algorithm"] = Rcpp::wrap("BFGS");
                       args["init_alpha"] = Rcpp::wrap(ctrl.optim.init_alpha);
                       args["tol_param"] = Rcpp::wrap(ctrl.optim.tol_param);
                       args["tol_obj"] = Rcpp::wrap(ctrl.optim.tol_obj);
                       args["tol_grad"] = Rcpp::wrap(ctrl.optim.tol_grad);
                       args["tol_rel_obj"] = Rcpp::wrap(ctrl.optim.tol_rel_obj);
                       args["tol_rel_grad"] = Rcpp::wrap(ctrl.optim.tol_rel_grad);
                       break;
          }
          break;
      }
      return Rcpp::wrap(args);
    }

    inline const std::string& get_sample_file() const {
      return sample_file;
    }
    inline bool get_sample_file_flag() const {
      return sample_file_flag;
    }
    inline bool get_diagnostic_file_flag() const {
      return diagnostic_file_flag;
    }
    inline const std::string& get_diagnostic_file() const {
      return diagnostic_file;
    }

    void set_random_seed(unsigned int seed) {
      random_seed = seed;
    }

    inline unsigned int get_random_seed() const {
      return random_seed;
    }

    inline bool get_append_samples() const {
      return append_samples;
    }
    inline stan_args_method_t get_method() const {
      return method;
    }
    inline int get_refresh() const {
      switch (method) {
        case OPTIM: return ctrl.optim.refresh;
      }
      return 0;
    }
    inline int get_iter() const {
      switch (method) {
        case OPTIM: return ctrl.optim.iter;
      }
      return 0;
    }

    inline optim_algo_t get_ctrl_optim_algorithm() const {
      return ctrl.optim.algorithm;
    }
    inline int get_ctrl_optim_refresh() const {
      return ctrl.optim.refresh;
    }
    inline bool get_ctrl_optim_save_iterations() const {
      return ctrl.optim.save_iterations;
    }
    inline double get_ctrl_optim_init_alpha() const {
      return ctrl.optim.init_alpha;
    }
    inline double get_ctrl_optim_tol_obj() const {
      return ctrl.optim.tol_obj;
    }
    inline double get_ctrl_optim_tol_grad() const {
      return ctrl.optim.tol_grad;
    }
    inline double get_ctrl_optim_tol_param() const {
      return ctrl.optim.tol_param;
    }
    inline double get_ctrl_optim_tol_rel_obj() const {
      return ctrl.optim.tol_rel_obj;
    }
    inline double get_ctrl_optim_tol_rel_grad() const {
      return ctrl.optim.tol_rel_grad;
    }
    inline int get_ctrl_optim_history_size() const {
      return ctrl.optim.history_size;
    }
    inline unsigned int get_chain_id() const {
      return chain_id;
    }
    inline double get_init_radius() const {
      return init_radius;
    }
    inline bool get_enable_random_init() const {
      return enable_random_init;
    }
    const std::string& get_init() const {
      return init;
    }
    SEXP get_init_list() const {
      return init_list;
    }

    void write_args_as_comment(std::ostream& ostream) const {
      write_comment_property(ostream,"init",init);
      write_comment_property(ostream,"enable_random_init",enable_random_init);
      write_comment_property(ostream,"seed",random_seed);
      write_comment_property(ostream,"chain_id",chain_id);
      write_comment_property(ostream,"iter",get_iter());
      switch (method) {
        case OPTIM:
          write_comment_property(ostream,"refresh",ctrl.optim.refresh);
          write_comment_property(ostream,"save_iterations",ctrl.optim.save_iterations);
          switch (ctrl.optim.algorithm) {
            case Newton: write_comment_property(ostream,"algorithm", "Newton"); break;
            case BFGS: write_comment_property(ostream,"algorithm", "BFGS");
                       write_comment_property(ostream,"init_alpha", ctrl.optim.init_alpha);
                       write_comment_property(ostream,"tol_obj", ctrl.optim.tol_obj);
                       write_comment_property(ostream,"tol_grad", ctrl.optim.tol_grad);
                       write_comment_property(ostream,"tol_param", ctrl.optim.tol_param);
                       write_comment_property(ostream,"tol_rel_obj", ctrl.optim.tol_rel_obj);
                       write_comment_property(ostream,"tol_rel_grad", ctrl.optim.tol_rel_grad);
                       break;
            case LBFGS: write_comment_property(ostream,"algorithm", "LBFGS");
                       write_comment_property(ostream,"init_alpha", ctrl.optim.init_alpha);
                       write_comment_property(ostream,"tol_obj", ctrl.optim.tol_obj);
                       write_comment_property(ostream,"tol_grad", ctrl.optim.tol_grad);
                       write_comment_property(ostream,"tol_param", ctrl.optim.tol_param);
                       write_comment_property(ostream,"tol_rel_obj", ctrl.optim.tol_rel_obj);
                       write_comment_property(ostream,"tol_rel_grad", ctrl.optim.tol_rel_grad);
                       write_comment_property(ostream,"history_size", ctrl.optim.history_size);
                       break;
          }
      }
      if (sample_file_flag)
        write_comment_property(ostream,"sample_file",sample_file);
      if (diagnostic_file_flag)
        write_comment_property(ostream,"diagnostic_file",diagnostic_file);
      write_comment_property(ostream,"append_samples",append_samples);
      write_comment(ostream);
    }
  };
}

#endif

