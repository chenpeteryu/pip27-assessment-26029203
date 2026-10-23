#include "t4_realtime_estimation/msg/signal_packet.hpp"
#include "prediction_stream.hpp"

#include <rcl_interfaces/msg/set_parameters_result.hpp>
#include <rclcpp/rclcpp.hpp>
#include <chrono>
#include <cmath>
#include <cstdint>
#include <memory>
#include <stdexcept>
#include <vector>

namespace signal_estimator {
using Packet = t4_realtime_estimation::msg::SignalPacket;
using SteadyClock = std::chrono::steady_clock;

class EstimatorNode final : public rclcpp::Node {
 public:
  EstimatorNode() : Node("signal_estimator") {
    mode_ = declare_parameter<std::int64_t>("assessment_mode", 0);
    if (mode_ < 0 || mode_ > 3) throw std::invalid_argument("assessment_mode must be 0..3");
    stream_ = std::make_unique<PredictionStream>(static_cast<int>(mode_));
    const auto qos = rclcpp::QoS(rclcpp::KeepLast(256)).reliable().durability_volatile();
    estimate_publisher_ = create_publisher<Packet>("/signal/estimate", qos);
    prediction_publisher_ = create_publisher<Packet>("/signal/prediction", qos);
    observation_subscription_ = create_subscription<Packet>(
        "/signal/observation", qos,
        [this](Packet::ConstSharedPtr message, const rclcpp::MessageInfo& info) {
          onObservation(*message, info);
        });
    timer_ = create_wall_timer(std::chrono::milliseconds(1), [this] { publishDue(); });
    parameter_callback_ = add_on_set_parameters_callback(
        [this](const std::vector<rclcpp::Parameter>& parameters) {
          rcl_interfaces::msg::SetParametersResult result;
          result.successful = true;
          std::int64_t requested = mode_;
          for (const auto& parameter : parameters) {
            if (parameter.get_name() != "assessment_mode") continue;
            if (parameter.get_type() != rclcpp::ParameterType::PARAMETER_INTEGER ||
                parameter.as_int() < 0 || parameter.as_int() > 3 || started_) {
              result.successful = false;
              result.reason = "assessment_mode requires integer 0..3 before the first observation";
              return result;
            }
            requested = parameter.as_int();
          }
          if (requested != mode_) {
            mode_ = requested;
            stream_ = std::make_unique<PredictionStream>(static_cast<int>(mode_));
          }
          return result;
        });
  }

 private:
  void onObservation(const Packet& message, const rclcpp::MessageInfo& info) {
    if (message.values.size() != 1 || !std::isfinite(message.values[0]) ||
        message.header.stamp.nanosec >= 1'000'000'000U) {
      RCLCPP_ERROR(get_logger(), "Invalid observation packet");
      return;
    }
    const std::int64_t stamp = static_cast<std::int64_t>(message.header.stamp.sec) *
                                   1'000'000'000LL + message.header.stamp.nanosec;
    if (!started_) {
      const auto source_ns = info.get_rmw_message_info().source_timestamp;
      if (mode_ != 0 && source_ns <= 0)
        throw std::runtime_error("DDS source_timestamp is required for deadline alignment");
      const auto system_now = std::chrono::system_clock::now().time_since_epoch();
      const auto system_ns = std::chrono::duration_cast<std::chrono::nanoseconds>(system_now).count();
      // DDS wall time only schedules output; Header time remains the signal clock.
      real_origin_ = SteadyClock::now() - std::chrono::nanoseconds(system_ns - source_ns);
      started_ = true;
    }
    if (mode_ == 0) {
      const double y = stream_->observeBasic(stamp, message.values[0]);
      Packet response;
      response.header = message.header;
      response.values = {y};
      estimate_publisher_->publish(response);
    } else {
      stream_->enqueue(stamp, message.values[0]);
      publishDue();
    }
  }

  void publishDue() {
    if (!started_ || mode_ == 0) return;
    for (int emitted = 0; emitted < 8; ++emitted) {
      const auto deadline = real_origin_ + std::chrono::nanoseconds(stream_->tick() * kStepNs);
      if (SteadyClock::now() < deadline) break;
      const auto prediction = stream_->next();
      Packet response;
      std::int64_t sec = prediction.stamp / 1'000'000'000LL;
      std::int64_t ns = prediction.stamp % 1'000'000'000LL;
      if (ns < 0) { --sec; ns += 1'000'000'000LL; }
      response.header.stamp.sec = static_cast<std::int32_t>(sec);
      response.header.stamp.nanosec = static_cast<std::uint32_t>(ns);
      response.values = prediction.values;
      prediction_publisher_->publish(response);
    }
  }

  std::int64_t mode_ = 0;
  bool started_ = false;
  SteadyClock::time_point real_origin_;
  std::unique_ptr<PredictionStream> stream_;
  rclcpp::Subscription<Packet>::SharedPtr observation_subscription_;
  rclcpp::Publisher<Packet>::SharedPtr estimate_publisher_, prediction_publisher_;
  rclcpp::TimerBase::SharedPtr timer_;
  rclcpp::node_interfaces::OnSetParametersCallbackHandle::SharedPtr parameter_callback_;
};
}  // namespace signal_estimator

int main(int argc, char** argv) {
  rclcpp::init(argc, argv);
  int status = 0;
  try {
    rclcpp::spin(std::make_shared<signal_estimator::EstimatorNode>());
  } catch (const std::exception& error) {
    RCLCPP_ERROR(rclcpp::get_logger("signal_estimator"), "%s", error.what());
    status = 1;
  }
  rclcpp::shutdown();
  return status;
}
