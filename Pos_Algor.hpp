// Positioning algorithms.
// C++11 extensions required.
#ifndef Pos_Algor_hpp
#define Pos_Algor_hpp

#include <stddef.h>
#include <stdint.h>

/// This class contains routines to maintain the position with the reading from
/// a MPU6050 accelerometer/gyroscope and other sources of position.
///
///
/// The first set of functions correct the error caused by a tilt of the sensor:
///
/// Suppose the reading is unbiased, we first detect the downward position
/// by calculating the resultant direction of acceleration at rest,
/// then the correct components can be calculated.
/// The magnitude of the acceleration due to gravity will also be recorded to
/// assist in deducing the net vertical acceleration.
///
/// However, the above method can only determine the direction of the corrected
/// z axis, while the directions of the corrected x and y axes remain unknown.
/// To solve this, we can only safely assume the physical x axis has a nearly
/// correct direction, and the corrected x axis will have the same direction as
/// the physical x axis when projected onto the physical xy plate.
/// If it is possible to generate a x-only thrust, another routine can also be
/// used to determine a more meaningful corrected x axis.
///
///
/// The kalman_* set of functions try to decrease the random error by utilizing
/// a Kalman filter (Kalman, 1960).
///
/// The integr_* set of functions are meant to do discrete integration by
/// smoothing the data to a parabola.
//
// TODO: Kalman filter, x_direction tests, GPS data source
template <typename ValueType> class Motion_State
{
  public:
    // In the rest test, all added data will go into _tmp_* with average.
    void start_rest_test();
    void end_rest_test();
    void start_x_direction_test();
    void end_x_direction_test();
    // Update MPU6050 reading
    void add_data(ValueType phy_ax, ValueType phy_ay, ValueType phy_az,
                  ValueType phy_gyx, ValueType phy_gyy, ValueType phy_gyz,
                  unsigned long time_now);
    // Update displacement from GPS
    void add_displacement(ValueType dx, ValueType dy);

    // Getters for direct data
    ValueType get_last_ax();
    ValueType get_last_ay();
    ValueType get_last_az();
    ValueType get_last_gyx();
    ValueType get_last_gyy();
    ValueType get_last_gyz();

    // --- EXPORTED DATA ---
    // Tilt with respect to x,y,z axes, in radians.
    // These are affected by neither turns nor pitching. i.e.,
    // The coordinate system is fixed one the calibration is done.
    ValueType x_tilt = .0;
    ValueType y_tilt = .0;
    ValueType z_tilt = .0;
    // x,y,z velocity, their direction never changes and are relative to
    // the predefined start point.
    ValueType x_velo = .0;
    ValueType y_velo = .0;
    ValueType z_velo = .0;
    // x,y,z displacement, their direction never changes and are relative to
    // the predefined start point.
    ValueType x_disp = .0;
    ValueType y_disp = .0;
    ValueType z_disp = .0;
    // Direction of acceleration along the corrected z axis, in radians
    ValueType _direction = .0;

  private:
    void _do_add_data_operational(ValueType phy_ax, ValueType phy_ay,
                                  ValueType phy_az, ValueType phy_gyx,
                                  ValueType phy_gyy, ValueType phy_gyz,
                                  unsigned long time_now);
    void _do_add_data_rest_test(ValueType phy_ax, ValueType phy_ay,
                                ValueType phy_az, ValueType phy_gyx,
                                ValueType phy_gyy, ValueType phy_gyz);
    void _do_add_data_x_test(ValueType phy_ax, ValueType phy_ay,
                             ValueType phy_az, ValueType phy_gyx,
                             ValueType phy_gyy, ValueType phy_gyz);
    // The current state is passed in to save memory
    void _integr_update(ValueType phy_ax, ValueType phy_ay, ValueType phy_az,
                        ValueType phy_gyx, ValueType phy_gyy, ValueType phy_gyz,
                        unsigned long time_cur);
    // This state determines how the data entered will be used
    enum class Motion_State_State
    {
        // Doing regular calculation
        Operational = 1,
        // Doing rest test
        RestTest,
        // Doing x-only test
        XTest,
    } _operation_state = Motion_State_State::Operational;
    // --- INTERNAL DATA ---
    // Calculated acceleration due to gravity
    ValueType _g = .0;
    // Phase of the current data, can be 0, 1 or 2.
    // integration happens whenever phase is 2. Then,
    // phase is reset to 1 because the last (third) data
    // automatically goes into [0].
    uint8_t _phase = 0;
    // Previous states for smoothing.
    // Time in milliseconds
    unsigned long _time[2];
    // The third field is used to store the previous C (constant
    // term of the antiderivative) of the acceleration in order to calculate the
    // displacement
    ValueType _prev_x_accel[3] = {.0};
    ValueType _prev_y_accel[3] = {.0};
    ValueType _prev_z_accel[3] = {.0};
    ValueType _prev_x_ang_velo[2] = {.0};
    ValueType _prev_y_ang_velo[2] = {.0};
    ValueType _prev_z_ang_velo[2] = {.0};
    // Temporary values for corrections
    size_t _tmp_data_idx = 0;
    ValueType _tmp_x_accel = .0;
    ValueType _tmp_y_accel = .0;
    ValueType _tmp_z_accel = .0;
    ValueType _tmp_x_ang_velo = .0;
    ValueType _tmp_y_ang_velo = .0;
    ValueType _tmp_z_ang_velo = .0;
};

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

#define defun_getter(name, valname)                                            \
    template <typename ValueType>                                              \
    ValueType Motion_State<ValueType>::get_last_##name()                       \
    {                                                                          \
        return (_prev_##valname)[_phase ^ 1];                                  \
    }
defun_getter(ax, x_accel);
defun_getter(ay, y_accel);
defun_getter(az, z_accel);
defun_getter(gyx, x_ang_velo);
defun_getter(gyy, y_ang_velo);
defun_getter(gyz, z_ang_velo);

// Add MPU6050 readings according to the state
template <typename ValueType>
void Motion_State<ValueType>::add_data(ValueType phy_ax, ValueType phy_ay,
                                       ValueType phy_az, ValueType phy_gyx,
                                       ValueType phy_gyy, ValueType phy_gyz,
                                       unsigned long time_now)
{
    switch (_operation_state)
    {
        case Motion_State_State::Operational:
            _do_add_data_operational(phy_ax, phy_ay, phy_az, phy_gyx, phy_gyy,
                                     phy_gyz, time_now);
            break;
        case Motion_State_State::RestTest:
            _do_add_data_rest_test(phy_ax, phy_ay, phy_az, phy_gyx, phy_gyy,
                                   phy_gyz);
            break;
        case Motion_State_State::XTest:
            _do_add_data_x_test(phy_ax, phy_ay, phy_az, phy_gyx, phy_gyy,
                                phy_gyz);
            break;
    }
}

// Add MPU6050 readings to be corrected and integrated
template <typename ValueType>
void Motion_State<ValueType>::_do_add_data_operational(
    ValueType phy_ax, ValueType phy_ay, ValueType phy_az, ValueType phy_gyx,
    ValueType phy_gyy, ValueType phy_gyz, unsigned long time_now)

{
    /* First, remove the gravity. Then, project onto the constructed axes */
    float ax = (phy_ax - cos(y_tilt) * _g) * cos(y_tilt) * cos(z_tilt),
          ay = (phy_ay - cos(z_tilt) * _g) * cos(x_tilt) * cos(z_tilt),
          az = (phy_az - cos(x_tilt) * _g) * cos(x_tilt) * cos(z_tilt);
    // TODO here: Kalman filter
    if (_phase == 2)
    {
        _integr_update(ax, ay, az, phy_gyx, phy_gyy, phy_gyz, time_now);
        _phase = 0;
    }
    _prev_x_accel[_phase] = ax;
    _prev_y_accel[_phase] = ay;
    _prev_z_accel[_phase] = az;
    _prev_x_ang_velo[_phase] = phy_gyx;
    _prev_y_ang_velo[_phase] = phy_gyy;
    _prev_z_ang_velo[_phase] = phy_gyz;
    _time[_phase] = time_now;
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
        (_tmp_##dir##_accel * (_tmp_data_idx - 1) + phy_a##dir) /              \
        _tmp_data_idx
#define update_gyros(dir)                                                      \
    _tmp_##dir##_ang_velo =                                                    \
        (_tmp_##dir##_ang_velo * (_tmp_data_idx - 1) + phy_gy##dir) /          \
        _tmp_data_idx
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
    /* TODO */
    (void)phy_ax;
    (void)phy_ay;
    (void)phy_az;
    (void)phy_gyx;
    (void)phy_gyy;
    (void)phy_gyz;
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
    auto int_first = [=](XType x) -> YType {
        return 1 / 3 * coeff_a * x * x * x + 1 / 2 * coeff_b * x * x +
               coeff_c * x;
    };
    results[0] = static_cast<ResultType>(int_first(x2) - int_first(x0));
    if (order_is_two)
    {
        auto c_first = *last_last - int_first(x0);
        *last_last = int_first(x2) + c_first;
        auto int_second = [=](XType x) -> YType {
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

    integral(_time[0] / 1000.0, _time[1] / 1000.0, time_cur / 1000.0,
             _prev_x_accel[0], _prev_x_accel[1], phy_ax, true,
             _prev_x_accel + 2, delta_x);
    integral(_time[0] / 1000.0, _time[1] / 1000.0, time_cur / 1000.0,
             _prev_y_accel[0], _prev_y_accel[1], phy_ay, true,
             _prev_y_accel + 2, delta_y);
    integral(_time[0] / 1000.0, _time[1] / 1000.0, time_cur / 1000.0,
             _prev_z_accel[0], _prev_z_accel[1], phy_az, true,
             _prev_z_accel + 2, delta_z);
    integral(_time[0] / 1000.0, _time[1] / 1000.0, time_cur / 1000.0,
             _prev_x_ang_velo[0], _prev_x_ang_velo[1], phy_gyx, false,
             (float *)nullptr, &delta_angx);
    integral(_time[0] / 1000.0, _time[1] / 1000.0, time_cur / 1000.0,
             _prev_y_ang_velo[0], _prev_y_ang_velo[1], phy_gyy, false,
             (float *)nullptr, &delta_angy);
    integral(_time[0] / 1000.0, _time[1] / 1000.0, time_cur / 1000.0,
             _prev_z_ang_velo[0], _prev_z_ang_velo[1], phy_gyz, false,
             (float *)nullptr, &delta_angz);
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
#endif
