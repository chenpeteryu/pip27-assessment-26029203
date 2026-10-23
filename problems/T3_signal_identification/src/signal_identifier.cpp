#include "signal_identifier.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <limits>
#include <utility>
#include <vector>

using namespace std;

namespace {
constexpr double kPi = 3.141592653589793238462643383279502884;
using Coefficients = array<double, 5>;

struct Fit {
  Coefficients coefficients{};
  double error = numeric_limits<double>::infinity();
};


Fit LeastSquares(const vector<double>& times,
                 const vector<double>& observations, int degree,
                 double frequency, double center, double scale) {
  const bool periodic = frequency > 0;
  const int columns = degree + 1 + (periodic ? 2 : 0);
  const size_t count = times.size();
  if (count < static_cast<size_t>(columns)) return {};
  array<vector<double>, 5> q;
  array<array<double, 5>, 5> r{};
  array<double, 5> projected{};
  for (int j = 0; j < columns; ++j) {
    q[j].resize(count);
    for (size_t i = 0; i < count; ++i) {
      const double x = (times[i] - center) / scale;
      if (j <= degree) {
        q[j][i] = j == 0 ? 1.0 : (j == 1 ? x : x * x);
      } else {
        const double angle = 2 * kPi * frequency * times[i];
        q[j][i] = j == degree + 1 ? sin(angle) : cos(angle);
      }
    }
    for (int pass = 0; pass < 2; ++pass) {
      for (int k = 0; k < j; ++k) {
        double dot = 0;
        for (size_t i = 0; i < count; ++i) dot += q[k][i] * q[j][i];
        r[k][j] += dot;
        for (size_t i = 0; i < count; ++i) q[j][i] -= dot * q[k][i];
      }
    }
    double norm = 0;
    for (double value : q[j]) norm += value * value;
    norm = sqrt(norm);
    if (norm < 1e-12) return {};
    r[j][j] = norm;
    for (size_t i = 0; i < count; ++i) {
      q[j][i] /= norm;
      projected[j] += q[j][i] * observations[i];
    }
  }
  Fit result;
  for (int j = columns - 1; j >= 0; --j) {
    double value = projected[j];
    for (int k = j + 1; k < columns; ++k) value -= r[j][k] * result.coefficients[k];
    result.coefficients[j] = value / r[j][j];
  }

  result.error = 0;
  for (size_t i = 0; i < count; ++i) {
    const double x = (times[i] - center) / scale;
    double prediction = result.coefficients[0];
    if (degree >= 1) prediction += result.coefficients[1] * x;
    if (degree >= 2) prediction += result.coefficients[2] * x * x;
    if (periodic) {
      const double angle = 2 * kPi * frequency * times[i];
      prediction += result.coefficients[degree + 1] * sin(angle) +
                    result.coefficients[degree + 2] * cos(angle);
    }
    const double residual = observations[i] - prediction;
    result.error += residual * residual;
  }
  return result;
}

double NormalizePhase(double phase) {
  return phase >= kPi ? phase - 2 * kPi : phase;
}
}

struct SignalIdentifier::Impl {
  explicit Impl(SignalMode signal_mode) : mode(signal_mode) {}
  SignalMode mode;
  double origin = 0;
  vector<double> times;
  vector<double> observations;
};

SignalIdentifier::SignalIdentifier(SignalMode mode, double noise_stddev)
    : impl_(make_unique<Impl>(mode)) {
  (void)noise_stddev;
}

SignalIdentifier::~SignalIdentifier() = default;

void SignalIdentifier::update(double timestamp, double observation) {
  if (impl_->times.empty()) impl_->origin = timestamp;
  impl_->times.push_back(timestamp - impl_->origin);
  impl_->observations.push_back(observation);
}

ParameterEstimate SignalIdentifier::estimate() const {
  ParameterEstimate result{};
  const auto& times = impl_->times;
  const auto& observations = impl_->observations;
  if (times.empty()) return result;
  const bool periodic = impl_->mode == SignalMode::Sine ||
                        impl_->mode == SignalMode::QuadraticSine;
  int degree = 0;
  if (impl_->mode == SignalMode::Linear) degree = 1;
  if (impl_->mode == SignalMode::Quadratic ||
      impl_->mode == SignalMode::QuadraticSine) degree = 2;
  const double duration = times.back();
  const double center = duration / 2;
  const double scale = duration > 0 ? duration / 2 : 1;
  auto fit = [&](double frequency) {
    return LeastSquares(times, observations, degree, frequency, center, scale);
  };
  double frequency = 0;
  Fit best;
  if (!periodic) {
    degree = min(degree, static_cast<int>(times.size() - 1));
    best = fit(0);
  } else if (duration > 0 && times.size() >= static_cast<size_t>(degree + 3)) {
    const double lower = max(0.01, 0.2 - 0.5 / duration);
    const double upper = 2.0 + 0.5 / duration;
    const int intervals = static_cast<int>(ceil((upper - lower) * 8 * duration));
    const double step = (upper - lower) / intervals;
    vector<double> errors(static_cast<size_t>(intervals) + 1);
    for (int i = 0; i <= intervals; ++i) {
      const double f = lower + i * step;
      const Fit candidate = fit(f);
      errors[i] = candidate.error;
      if (candidate.error < best.error) {
        best = candidate;
        frequency = f;
      }
    }
    vector<pair<double, int>> minima;
    for (int i = 0; i <= intervals; ++i) {
      if ((i == 0 || errors[i] <= errors[i - 1]) &&
          (i == intervals || errors[i] <= errors[i + 1])) {
        minima.emplace_back(errors[i], i);
      }
    }
    sort(minima.begin(), minima.end());
    const size_t candidates = min<size_t>(8, minima.size());
    for (size_t index = 0; index < candidates; ++index) {
      const int i = minima[index].second;
      double left = lower + max(0, i - 1) * step;
      double right = lower + min(intervals, i + 1) * step;
      constexpr double kRatio = 0.6180339887498948482;
      double f1 = right - kRatio * (right - left);
      double f2 = left + kRatio * (right - left);
      Fit fit1 = fit(f1), fit2 = fit(f2);
      for (int iteration = 0; iteration < 64; ++iteration) {
        if (fit1.error <= fit2.error) {
          right = f2;
          f2 = f1;
          fit2 = fit1;
          f1 = right - kRatio * (right - left);
          fit1 = fit(f1);
        } else {
          left = f1;
          f1 = f2;
          fit1 = fit2;
          f2 = left + kRatio * (right - left);
          fit2 = fit(f2);
        }
      }
      for (const auto& candidate : {make_pair(f1, fit1), make_pair(f2, fit2)}) {
        if (candidate.second.error < best.error) {
          frequency = candidate.first;
          best = candidate.second;
        }
      }
    }
  } else {
    best.coefficients[0] = observations.front();
    frequency = 0.2;
  }
  const auto& b = best.coefficients;
  auto& p = result.parameters;
  p[0] = b[0];
  if (degree >= 1) p[0] -= b[1] * center / scale;
  if (degree >= 2) p[0] += b[2] * center * center / (scale * scale);
  if (degree >= 1) p[1] = b[1] / scale - 2 * b[2] * center / (scale * scale);
  if (degree >= 2) p[2] = 2 * b[2] / (scale * scale);
  if (periodic) {
    const int offset = degree + 1;
    p[offset] = max(hypot(b[offset], b[offset + 1]),
                         numeric_limits<double>::min());
    p[offset + 1] = frequency;
    p[offset + 2] = NormalizePhase(atan2(b[offset + 1], b[offset]));
  }
  return result;
}
