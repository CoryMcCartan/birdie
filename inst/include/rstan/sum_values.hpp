#ifndef RSTAN_SUM_VALUES_HPP
#define RSTAN_SUM_VALUES_HPP

#include <stan/callbacks/writer.hpp>
#include <stdexcept>
#include <string>
#include <vector>

namespace rstan {
  class sum_values : public stan::callbacks::writer {
  public:
    explicit sum_values(const size_t N)
      : N_(N), m_(0), skip_(0), sum_(N_, 0.0) { }

    sum_values(const size_t N, const size_t skip)
      : N_(N), m_(0), skip_(skip), sum_(N_, 0.0) { }

    // To deal with C++ name hiding
    using stan::callbacks::writer::operator();

    /**
     * Add values to cumulative sum
     *
     * @param x vector of type T
     */
    void operator()(const std::vector<double>& state) {
      if (N_ != state.size())
        throw std::length_error("vector provided does not "
                                "match the parameter length");
      if (m_ >= skip_) {
        for (size_t n = 0; n < N_; n++)
          sum_[n] += state[n];
      }
      m_++;
    }

    const std::vector<double>& sum() const {
      return sum_;
    }

    const size_t called() const {
      return m_;
    }

    const size_t recorded() const {
      if (m_ >= skip_)
        return m_ - skip_;
      else
        return 0;
    }

  private:
    size_t N_;
    size_t m_;
    size_t skip_;
    std::vector<double> sum_;
  };

}
#endif
