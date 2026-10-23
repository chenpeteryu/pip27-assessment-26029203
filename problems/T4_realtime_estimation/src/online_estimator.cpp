#include "online_estimator.hpp"

#include <algorithm>
#include <cmath>
#include <limits>
#include <stdexcept>

namespace signal_estimator {

double OnlineEstimator::kernel(double delta, double length) const {
  const double r = delta / length;
  // Smooth local dynamics plus an uncertain constant offset.
  return amplitude2_ * (std::exp(-0.5 * r * r) + 10.0);
}

double OnlineEstimator::fit(double length, bool save) {
  const std::size_t n = times_.size();
  std::vector<double> lower(n * n, 0.0), z(n), alpha(n);
  double logdet = 0.0;
  for (std::size_t i = 0; i < n; ++i) {
    for (std::size_t j = 0; j <= i; ++j) {
      double s = kernel(times_[i] - times_[j], length);
      if (i == j) s += noise2_ + amplitude2_ * 1e-9;
      for (std::size_t k = 0; k < j; ++k)
        s -= lower[i * n + k] * lower[j * n + k];
      if (i == j) {
        if (!(s > 0.0)) return std::numeric_limits<double>::infinity();
        lower[i * n + j] = std::sqrt(s);
        logdet += std::log(s);
      } else {
        lower[i * n + j] = s / lower[j * n + j];
      }
    }
    double s = values_[i] - mean_;
    for (std::size_t k = 0; k < i; ++k) s -= lower[i * n + k] * z[k];
    z[i] = s / lower[i * n + i];
  }
  double quadratic = 0.0;
  for (std::size_t ri = n; ri > 0; --ri) {
    const std::size_t i = ri - 1;
    quadratic += z[i] * z[i];
    double s = z[i];
    for (std::size_t k = i + 1; k < n; ++k) s -= lower[k * n + i] * alpha[k];
    alpha[i] = s / lower[i * n + i];
  }
  if (save) alpha_ = std::move(alpha);
  return quadratic + logdet;
}

void OnlineEstimator::addObservation(std::int64_t stamp_ns, double value) {
  if (!std::isfinite(value) || (!history_.empty() && stamp_ns <= history_.back().stamp))
    throw std::invalid_argument("Observations must be finite and strictly ordered");
  if (history_.empty()) origin_ = stamp_ns;
  history_.push_back({stamp_ns, value});
  while (history_.size() > 512) history_.pop_front();
  ++updates_;

  // A nonuniform second difference removes local linear motion. Divide by the
  // sum of squared coefficients, rather than assuming perfectly uniform input.
  double residual2 = 0.0, normalization = 0.0;
  const std::size_t begin = history_.size() > 256 ? history_.size() - 256 : 0;
  for (std::size_t i = begin + 2; i < history_.size(); ++i) {
    const auto& a = history_[i - 2];
    const auto& b = history_[i - 1];
    const auto& c = history_[i];
    const double w = static_cast<double>(b.stamp - a.stamp) /
                     static_cast<double>(c.stamp - a.stamp);
    const double r = b.value - (1.0 - w) * a.value - w * c.value;
    residual2 += r * r;
    normalization += 1.0 + w * w + (1.0 - w) * (1.0 - w);
  }
  if (normalization > 0.0) noise2_ = std::max(1e-10, residual2 / normalization);

  // Bound the dense system. Retain the newest samples and spread older ones
  // across a mode-specific physical-time window, independent of timestamp origin.
  const double window = mode_ == 0 ? 1.5 : 3.0;
  std::vector<Sample> selected;
  for (const auto& sample : history_)
    if (static_cast<double>(stamp_ns - sample.stamp) * 1e-9 <= window)
      selected.push_back(sample);
  times_.clear(); values_.clear();
  const std::size_t n = std::min<std::size_t>(192, selected.size());
  mean_ = 0.0;
  for (std::size_t i = 0; i < n; ++i) {
    const auto& s = selected[n == 1 ? 0 : i * (selected.size() - 1) / (n - 1)];
    times_.push_back(static_cast<double>(s.stamp - origin_) * 1e-9);
    values_.push_back(s.value);
    mean_ += s.value / static_cast<double>(n);
  }
  double variance = 0.0;
  for (double y : values_) variance += (y - mean_) * (y - mean_) / n;
  amplitude2_ = std::max(1.0, 2.0 * std::max(0.0, variance - noise2_));
  if (updates_ == 1) length_ = mode_ == 0 ? 0.35 : 0.75;
  if (updates_ % 10 == 0 && n >= 12) {
    const double base = mode_ == 0 ? 0.26 : (mode_ == 3 ? 0.60 : 0.50);
    double best = std::numeric_limits<double>::infinity();
    for (double factor : {1.0, 1.3, 1.7, 2.2, 2.9}) {
      const double length = base * factor;
      const double score = fit(length, false);
      if (score < best) { best = score; length_ = length; }
    }
  }
  fit(length_, true);
}

double OnlineEstimator::estimate(std::int64_t stamp_ns) const {
  if (history_.empty()) return 0.0;
  if (stamp_ns < history_.back().stamp)
    throw std::invalid_argument("Prediction would use future observations");
  const double t = static_cast<double>(stamp_ns - origin_) * 1e-9;
  double result = mean_;
  for (std::size_t i = 0; i < alpha_.size(); ++i)
    result += kernel(t - times_[i], length_) * alpha_[i];
  return result;
}

std::vector<double> OnlineEstimator::predict(std::int64_t start_ns, std::size_t count) const {
  std::vector<double> result(count);
  for (std::size_t i = 0; i < count; ++i)
    result[i] = estimate(start_ns + static_cast<std::int64_t>(i) * 5'000'000);
  return result;
}
}  // namespace signal_estimator
