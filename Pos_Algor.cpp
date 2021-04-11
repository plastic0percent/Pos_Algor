#include <cstdint>

#include "Pos_Algor.h"

// Mark the start of rest test
template <typename ValueType> void Motion_State<ValueType>::start_rest_test()
{
    _operation_state = Motion_State_State::RestTest;
}

// Mark the end of rest test, state is reset to operational
template <typename ValueType> void Motion_State<ValueType>::end_rest_test()
{
    _g =
        sqrt((_tmp_x_accel) * (_tmp_x_accel) + (_tmp_y_accel) * (_tmp_y_accel) +
             (_tmp_z_accel) * (_tmp_z_accel));
    x_tilt = acos(_tmp_z_accel / _g);
    y_tilt = acos(_tmp_x_accel / _g);
    z_tilt = acos(_tmp_y_accel / _g);
    // Reset states
    _tmp_data_idx = 0;
    _tmp_x_accel = .0;
    _tmp_y_accel = .0;
    _tmp_z_accel = .0;
    _tmp_x_ang_velo = .0;
    _tmp_y_ang_velo = .0;
    _tmp_z_ang_velo = .0;
    _operation_state = Motion_State_State::Operational;
}

// Mark the start of x-direction test
template <typename ValueType>
void Motion_State<ValueType>::start_x_direction_test()
{
    _operation_state = Motion_State_State::XTest;
}

// Mark the end of x-direction test, state is reset to operational
template <typename ValueType>
void Motion_State<ValueType>::end_x_direction_test()
{
    _operation_state = Motion_State_State::Operational;
}

// Add MPU6050 readings according to the state
template <typename ValueType>
void Motion_State<ValueType>::add_data(ValueType phy_ax, ValueType phy_ay,
                                       ValueType phy_az, ValueType phy_gyx,
                                       ValueType phy_gyy, ValueType phy_gyz)
{
    switch (_operation_state)
    {
        case Motion_State_State::Operational:
            _add_data_operational(phy_ax, phy_ay, phy_az, phy_gyx, phy_gyy,
                                  phy_gyz);
            break;
        case Motion_State_State::RestTest:
            _add_data_rest_test(phy_ax, phy_ay, phy_az, phy_gyx, phy_gyy,
                                phy_gyz);
            break;
        case Motion_State_State::XTest:
            _add_data_x_test(phy_ax, phy_ay, phy_az, phy_gyx, phy_gyy, phy_gyz);
            break;
    }
}

// Add MPU6050 readings to be corrected and integrated
template <typename ValueType>
void Motion_State<ValueType>::_do_add_data_operational(
    ValueType phy_ax, ValueType phy_ay, ValueType phy_az, ValueType phy_gyx,
    ValueType phy_gyy, ValueType phy_gyz)

{
    float ax = abs(phy_ax) * cos(y_tilt) * cos(z_tilt),
          ay = abs(phy_ay) * cos(x_tilt) * cos(z_tilt),
          az = abs(phy_az) * cos(x_tilt) * cos(z_tilt);
    // TODO here: Kalman filter
    if (_phase == 2)
    {
        _integr_update(ax, ay, az, phy_gyx, phy_gyy, phy_gyz);
        _phase = 0;
    }
    _prev_x_accel[_phase] = ax;
    _prev_y_accel[_phase] = ay;
    _prev_z_accel[_phase] = az;
    _prev_x_ang_velo[_phase] = phy_gyx;
    _prev_y_ang_velo[_phase] = phy_gyy;
    _prev_z_ang_velo[_phase] = phy_gyz;
    ++_phase;
}

// Add MPU6050 readings for rest correction
template <typename ValueType>
void Motion_State<ValueType>::_do_add_data_rest_test(
    ValueType phy_ax, ValueType phy_ay, ValueType phy_az, ValueType phy_gyx,
    ValueType phy_gyy, ValueType phy_gyz)

{
// Take the average. new_average = (old_average * (n-1) + new_data)/n
// This approach is used instead of storing in an array is due to two reasons:
// 1. Saving memory
// 2. Distributing latency of all operations
#define update_accel(dir)                                                      \
    _tmp_##dir##_accel =                                                       \
        (_tmp_data_idx * (_tmp_data_idx - 1) + phy_a##dir) / _tmp_data_idx
#define update_gyros(dir)                                                      \
    _tmp_##dir##_ang_velo =                                                    \
        (_tmp_data_idx * (_tmp_data_idx - 1) + phy_gy##dir) / _tmp_data_idx
    ++_tmp_data_idx;
    update_accel(x);
    update_accel(y);
    update_accel(z);
    update_gyros(x);
    update_gyros(y);
    update_gyros(z);
}

// Add MPU6050 readings for x-direction test
template <typename ValueType>
void Motion_State<ValueType>::_do_add_data_x_test(
    ValueType phy_ax, ValueType phy_ay, ValueType phy_az, ValueType phy_gyx,
    ValueType phy_gyy, ValueType phy_gyz)

{
}

template <typename XType, typename YType, typename ResultType>
// Do calculate an integration after smoothing to a parabola
// if order_is_two, a second-order integration is also calculated.
// *last_last is used to deduce the constant term. It records the previous
// y-coordinate of the last first-order integral, and will be updated with
// that of the newer function. If order_is_two is false, it is ignored.
// results are the definite integral between x0 and x3, and the length of
// the buffer should be at least equal to the order
static void integral(XType x0, XType x1, XType x2, YType y0, YType y1, YType y2,
                     bool order_is_two, YType *last_last, ResultType *results)
{
    // Smoothed functions
    auto coeff_a = (y1 - y2) / ((x1 - x2) * (x2 - x0)) -
                   (y1 - y0) / ((x1 - x0) * (x2 - x0)),
         coeff_b = (y1 - y0) / (x1 - x0) - (x1 + x0) * coeff_a,
         coeff_c = y2 - coeff_a * x2 * x2 - coeff_b * x2;
    auto int_first = [](XType x) -> YType {
        return 1 / 3 * coeff_a * x * x * x + 1 / 2 * coeff_b * x * x +
               coeff_c * x;
    };
    results[0] = static_cast<ResultType>(int_first(x2) - int_first(x0));
    if (order_is_two)
    {
        auto c_first = *last_last - int_first(x0);
        *last_last = int_first(x2) + c_first;
        auto int_second = [](XType x) -> YType {
            return 1 / 12 * coeff_a * x * x * x * x +
                   1 / 6 * coeff_b * x * x * x + 1 / 2 * coeff_c * x * x +
                   c_first * x;
        };
        results[1] = int_second(x2) - int_second(x0);
    }
}
// Do integration
template <typename ValueType>
void Motion_State<ValueType>::_integr_update(
    ValueType phy_ax, ValueType phy_ay, ValueType phy_az, ValueType phy_gyx,
    ValueType phy_gyy, ValueType phy_gyz, unsigned long time_cur)

{
    ValueType delta_x[2], delta_y[2], delta_z[2], delta_angx, delta_angy,
        delta_angz;

    integral(_time[0], _time[1], time_cur, _prev_x_accel[0], _prev_x_accel[1],
             phy_ax, true, _prev_x_accel + 2, delta_x);
    integral(_time[0], _time[1], time_cur, _prev_y_accel[0], _prev_y_accel[1],
             phy_ay, true, _prev_y_accel + 2, delta_y);
    integral(_time[0], _time[1], time_cur, _prev_z_accel[0], _prev_z_accel[1],
             phy_az, true, _prev_z_accel + 2, delta_z);
    integral(_time[0], _time[1], time_cur, _prev_x_ang_velo[0],
             _prev_x_ang_velo[1], phy_gyx, false, nullptr, &delta_angx);
    integral(_time[0], _time[1], time_cur, _prev_y_ang_velo[0],
             _prev_y_ang_velo[1], phy_gyy, false, nullptr, &delta_angy);
    integral(_time[0], _time[1], time_cur, _prev_z_ang_velo[0],
             _prev_z_ang_velo[1], phy_gyz, false, nullptr, &delta_angz);
    x_velo += delta_x[0];
    y_velo += delta_y[0];
    z_velo += delta_z[0];
    x_disp += delta_x[1];
    y_disp += delta_y[1];
    z_disp += delta_z[1];
    x_tilt += delta_angx;
    y_tilt += delta_angy;
    z_tilt += delta_angz;
}