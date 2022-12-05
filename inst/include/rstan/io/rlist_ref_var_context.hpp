#ifndef RSTAN__IO__RLIST_REF_VAR_CONTEXT_HPP
#define RSTAN__IO__RLIST_REF_VAR_CONTEXT_HPP

#include <cstddef>
#include <stdexcept>
#include <map>
#include <vector>
#include <utility>
#include <string>
#include <sstream>
#include <iostream>
#include <locale>

#include <boost/throw_exception.hpp>
#include <stan/io/var_context.hpp>
#include <Rcpp.h>

namespace rstan {

  namespace io {

    namespace {
      /**
      size_t product(std::vector<size_t> dims) {
        size_t y = 1U;
        for (size_t i = 0; i < dims.size(); ++i)
          y *= dims[i];
        return y;
      }
      */

      template <class T1, class T2>
      void  T1v_to_T2v(const std::vector<T1>& v,
                       std::vector<T2>& v2) {
        v2.resize(0);
        for (typename std::vector<T1>::const_iterator it = v.begin();
             it != v.end();
             ++it) {
          v2.push_back(static_cast<T2>(*it));
        }
      }
    }

    /**
     * Represents named arrays with dimensions.
     *
     * A rlist_ref_var_context implements var_context from
     * a R list (Rcpp::List) --- named variables (scalar,
     * vector, array) with dimensions.  The values for
     * an array are typed to double or int.  However,
     * it is R's job to pass data with correct types
     * though R uses double as the default atomic type.
     * Instead of copying the data as in
     * <code>rlist_var_context</code>, data are
     * obtained by references.
     *
     * <p>The dimensions and values of variables
     * may be accessed by name.
     *
     * <p> Difference from <code>rlist_var_context</code>:
     * here, when the object is constructed, we only
     * keep in the information of dimensions (the same as
     * in rlist_var_context) but the data are not copied. The
     * data are accessed directly from the R list.
     * T
     */
    class rlist_ref_var_context : public stan::io::var_context {
    private:
      Rcpp::List rlist_;
      std::map<std::string, std::vector<size_t> > vars_r_dim_;
      std::map<std::string, std::vector<size_t> > vars_i_dim_;
      std::vector<double> const empty_vec_r_;
      std::vector<int> const empty_vec_i_;
      std::vector<size_t> const empty_vec_ui_;
      /**
       * Return <code>true</code> if this rlist_ref_var_context contains the
       * specified variable name is defined in the real values. This method
       * returns <code>false</code> if the values are all integers.
       *
       * @param name Variable name to test.
       * @return <code>true</code> if the variable exists in the
       * real values of the rlist_ref_var_context.
       */
      bool contains_r_only(const std::string& name) const {
        return vars_r_dim_.count(name) > 0;
      }
    public:

      /**
       * Construct a rlist_ref__var_context object from the specified R list.
       *
       * @param in Input of R list (represented by SEXP)
       * from which to read.
       */
      rlist_ref_var_context(SEXP in): rlist_(in) {

        if (0 == rlist_.size()) return;
        std::vector<std::string> varnames
          = Rcpp::as<std::vector<std::string> >(rlist_.names());
        for (int i = 0; i < rlist_.size(); i++) {
          SEXP ee = rlist_[i];
          SEXP dim = Rf_getAttrib(ee, R_DimSymbol);
          R_len_t eelen = Rf_length(ee);

          // Note that in R, the default is real, which causes problems as we
          // are thinking they are integers, but Rf_isInteger would return
          // FALSE. One solution is use L suffix. Anyway, the R code should
          // change the data to integers when it is convertible.

          typedef std::map<std::string, std::vector<size_t> >::value_type
            psd_v_t; // a Pair of String and Dimensions

          if (Rf_isInteger(ee)) {
            if (Rf_length(dim) > 0) {
              std::vector<size_t> d;
              T1v_to_T2v(Rcpp::as<std::vector<unsigned int> >(dim), d);
              vars_i_dim_.insert(psd_v_t(varnames[i], d));
            } else {
              if (1 == eelen) {  // scalar
                vars_i_dim_.insert(psd_v_t(varnames[i], empty_vec_ui_));
              } else { // vector
                vars_i_dim_.insert(psd_v_t(varnames[i],
                                           std::vector<size_t>(1, eelen)));
              }
            }
          } else if(Rf_isNumeric(ee)) {
            if (Rf_length(dim) > 0) {
              std::vector<size_t> d;
              T1v_to_T2v(Rcpp::as<std::vector<unsigned int> >(dim), d);
              vars_r_dim_.insert(psd_v_t(varnames[i], d));
            } else {
              if (1 == eelen) {  // scalar
                vars_r_dim_.insert(psd_v_t(varnames[i], empty_vec_ui_));
              } else {
                vars_r_dim_.insert(psd_v_t(varnames[i],
                                           std::vector<size_t>(1, eelen)));
              }
            }
          } else {
            continue; // ignore non-numeric data
            // or should we report error?
          }
        }
      }

      /**
       * Return <code>true</code> if this rlist_ref_var_context contains the
       * specified variable name is defined. This method returns
       * <code>true</code> even if the values are all integers.
       *
       * @param name Variable name to test.
       * @return <code>true</code> if the variable exists.
       */
      bool contains_r(const std::string& name) const {
        return contains_r_only(name) || contains_i(name);
      }

      /**
       * Return <code>true</code> if this rlist_ref_var_context contains an integer
       * valued array with the specified name.
       *
       * @param name Variable name to test.
       * @return <code>true</code> if the variable name has an integer
       * array value.
       */
      bool contains_i(const std::string& name) const {
        return vars_i_dim_.count(name) > 0;
      }

      /**
       * Return the double values for the variable with the specified
       * name or null.
       *
       * @param name Name of variable.
       * @return Values of variable.
       */
      std::vector<double> vals_r(const std::string& name) const {
        if (contains_r(name)) {
          SEXP ee = rlist_[name];
          return Rcpp::as<std::vector<double> >(ee);
        }
        return empty_vec_r_;
      }

      /**
       * Return the dimensions for the double variable with the specified
       * name.
       *
       * @param name Name of variable.
       * @return Dimensions of variable.
       */
      std::vector<size_t> dims_r(const std::string& name) const {
        if (contains_r_only(name)) {
          return vars_r_dim_.find(name)->second;
        } else if (contains_i(name)) {
          return vars_i_dim_.find(name)->second;
        }
        return empty_vec_ui_;
      }

      /**
       * Return the integer values for the variable with the specified
       * name.
       *
       * @param name Name of variable.
       * @return Values.
       */
      std::vector<int> vals_i(const std::string& name) const {
        if (contains_i(name)) {
          SEXP ee = rlist_[name];
          return Rcpp::as<std::vector<int> >(ee);
        }
        return empty_vec_i_;
      }

      /**
       * Return the dimensions for the integer variable with the specified
       * name.
       *
       * @param name Name of variable.
       * @return Dimensions of variable.
       */
      std::vector<size_t> dims_i(const std::string& name) const {
        if (contains_i(name)) {
          return vars_i_dim_.find(name)->second;
        }
        return empty_vec_ui_;
      }

      /**
       * Return a list of the names of the floating point variables in
       * the rlist_ref_var_context.
       *
       * @param names Vector to store the list of names in.
       */
      virtual void names_r(std::vector<std::string>& names) const {
        names.resize(0);
        for (std::map<std::string, std::vector<size_t> >
                 ::const_iterator it = vars_r_dim_.begin();
             it != vars_r_dim_.end(); ++it)
          names.push_back((*it).first);
      }

      /**
       * Return a list of the names of the integer variables in
       * the rlist_ref_var_context.
       *
       * @param names Vector to store the list of names in.
       */
      virtual void names_i(std::vector<std::string>& names) const {
        names.resize(0);
        for (std::map<std::string, std::vector<size_t> >
                 ::const_iterator it = vars_i_dim_.begin();
             it != vars_i_dim_.end(); ++it)
          names.push_back((*it).first);
      }

      bool remove(const std::string& name) {
        return (vars_i_dim_.erase(name) > 0)
          || (vars_r_dim_.erase(name) > 0);
      }

    };

  }
}

#endif
