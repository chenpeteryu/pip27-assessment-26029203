#pragma once
#include "online_estimator.hpp"
#include <cstdint>
#include <deque>
#include <stdexcept>
#include <vector>

namespace signal_estimator {
constexpr std::int64_t kStepNs = 5'000'000;
struct Prediction { std::int64_t stamp; std::vector<double> values; };

class PredictionStream {
 public:
  explicit PredictionStream(int mode) : mode_(mode), estimator_(mode) {
    if (mode < 0 || mode > 3) throw std::invalid_argument("Invalid mode");
  }
  double observeBasic(std::int64_t stamp, double value) {
    estimator_.addObservation(stamp, value);
    return estimator_.estimate(stamp);
  }
  void enqueue(std::int64_t stamp, double value) {
    if (initialized_ && stamp <= last_input_) throw std::invalid_argument("Unordered input");
    if (pending_.size() == 256) throw std::runtime_error("Observation queue overflow");
    if (!initialized_) { origin_ = stamp; initialized_ = true; }
    last_input_ = stamp;
    pending_.push_back({stamp, value});
  }
  std::int64_t tick() const { return tick_; }
  Prediction next() {
    if (!initialized_) throw std::logic_error("No observation received");
    const auto stamp = origin_ + tick_ * kStepNs;
    while (!pending_.empty() && pending_.front().stamp <= stamp) {
      const auto sample = pending_.front();
      pending_.pop_front();
      estimator_.addObservation(sample.stamp, sample.value);
    }
    const std::size_t count = mode_ == 1 ? 21 : (mode_ == 2 ? 41 : 101);
    Prediction result{stamp, estimator_.predict(stamp, count)};
    ++tick_;
    return result;
  }
 private:
  struct Sample { std::int64_t stamp; double value; };
  int mode_;
  OnlineEstimator estimator_;
  std::deque<Sample> pending_;
  bool initialized_ = false;
  std::int64_t origin_ = 0, last_input_ = 0, tick_ = 0;
};
}  // namespace signal_estimator
