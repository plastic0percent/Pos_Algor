// Positioning algorithms.
// C++11 extensions required.
#ifndef Pos_Algor_h
#define Pos_Algor_h

#include <cstddef>

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
                  ValueType phy_gyx, ValueType phy_gyy, ValueType phy_gyz);
    // Update displacement from GPS
    void add_displacement(ValueType dx, ValueType dy);

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
                                  ValueType phy_gyy, ValueType phy_gyz);
    void _do_add_data_rest_test(ValueType phy_ax, ValueType phy_ay,
                                ValueType phy_az, ValueType phy_gyx,
                                ValueType phy_gyy, ValueType phy_gyz);
    void _do_add_data_x_test(ValueType phy_ax, ValueType phy_ay,
                             ValueType phy_az, ValueType phy_gyx,
                             ValueType phy_gyy, ValueType phy_gyz);
    // The current state is passed in to save memory
    void _integr_update(ValueType phy_ax, ValueType phy_ay, ValueType phy_az,
                        ValueType phy_gyx, ValueType phy_gyy,
                        ValueType phy_gyz, unsigned long time_cur);
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

#endif
