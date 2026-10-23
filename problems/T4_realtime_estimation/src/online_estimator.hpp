#pragma once

#include <array>
#include <cstdint>
#include <deque>
#include <vector>

namespace signal_estimator {

// Independent of ROS: the driver and node use exactly the same causal core.
class OnlineEstimator {
 public:
  explicit OnlineEstimator(int mode = 0) : mode_(mode) {}
  void addObservation(std::int64_t stamp_ns, double value);
  double estimate(std::int64_t stamp_ns) const;
  std::vector<double> predict(std::int64_t start_ns, std::size_t count) const;

 private:
  struct Sample { std::int64_t stamp; double value; };
  double kernel(double delta, double length) const;
  double fit(double length, bool save);
  int mode_;
  std::deque<Sample> history_;
  std::vector<double> times_, values_, alpha_;
  std::int64_t origin_ = 0;
  std::size_t updates_ = 0;
  double mean_ = 0.0;
  double amplitude2_ = 4.0;
  double noise2_ = 1.0;
  double length_ = 0.4;
};
}  // namespace signal_estimator
