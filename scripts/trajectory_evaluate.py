import json
import numpy as np
import pandas as pd
import os
from scipy.spatial import KDTree
from datetime import datetime


def read_json(json_path):
    with open(json_path, 'r') as f:
        data = json.load(f)
    points = data["Trajectory"]["PointList"]["Point"]
    tra = []
    for point in points:
        ts = point["Time"]
        td = datetime.strptime(ts, '%Y-%m-%d %H:%M:%S')
        td = int(td.timestamp())
        position = point["Location"]
        ps = position.split(",")
        x = float(ps[0])
        y = float(ps[1])
        z = float(ps[2])
        p = (td,x,y,z)
        tra.append(p)
    return tra



def find_closest_spatial_point(ground_truth, measured):
    """
    Find the closest point in ground truth (possibly interpolated) for each measured point.

    Args:
        ground_truth: List of (t, x, y, z) tuples, sorted by time.
        measured: List of (t, x, y, z) tuples.

    Returns:
        List of (t_measure, x_closest, y_closest, z_closest) tuples.
    """
    # Extract ground truth data
    gt_t = np.array([p[0] for p in ground_truth])
    gt_x = np.array([p[1] for p in ground_truth])
    gt_y = np.array([p[2] for p in ground_truth])
    gt_z = np.array([p[3] for p in ground_truth])

    # Interpolate ground truth to higher density (e.g., 0.1s steps)
    t_min, t_max = gt_t.min(), gt_t.max()
    t_interp = np.arange(t_min, t_max, 0.1)  # Adjust step size as needed
    x_interp = np.interp(t_interp, gt_t, gt_x)
    y_interp = np.interp(t_interp, gt_t, gt_y)
    z_interp = np.interp(t_interp, gt_t, gt_z)

    # Compute closest point for each measured point
    closest_points = []
    for t_m, x_m, y_m, z_m in measured:
        # Compute Euclidean distances to all interpolated GT points
        distances = np.sqrt(
            (x_interp - x_m) ** 2 +
            (y_interp - y_m) ** 2 +
            (z_interp - z_m) ** 2
        )
        # Find index of closest point
        idx_closest = np.argmin(distances)
        closest_points.append((
            t_m,
            x_interp[idx_closest],
            y_interp[idx_closest],
            z_interp[idx_closest]
        ))

    return closest_points


def find_closest_points(ground_truth, measured):
    """
    Find the closest ground truth point for each measured point using KDTree,
    calculate Euclidean distance errors, and compute the mean error.

    Args:
        ground_truth: List of (t, x, y, z) tuples representing ground truth trajectory
        measured: List of (t, x, y, z) tuples representing measured trajectory

    Returns:
        closest_points: List of closest ground truth points for each measured point
        errors: List of Euclidean distance errors for each point
        mean_error: Average error across entire trajectory
    """

    # Extract ground truth data into separate arrays
    gt_t = np.array([p[0] for p in ground_truth])  # Time stamps
    gt_x = np.array([p[1] for p in ground_truth])  # X coordinates
    gt_y = np.array([p[2] for p in ground_truth])  # Y coordinates
    gt_z = np.array([p[3] for p in ground_truth])  # Z coordinates

    # Create higher-density ground truth trajectory through interpolation
    # Using 0.1s intervals between original points
    t_interp = np.arange(gt_t.min(), gt_t.max(), 0.1)
    x_interp = np.interp(t_interp, gt_t, gt_x)  # Interpolated X coordinates
    y_interp = np.interp(t_interp, gt_t, gt_y)  # Interpolated Y coordinates
    z_interp = np.interp(t_interp, gt_t, gt_z)  # Interpolated Z coordinates

    # Build KDTree for efficient nearest neighbor searches
    # Stack interpolated coordinates into Nx3 array
    gt_points = np.column_stack((x_interp, y_interp, z_interp))
    tree = KDTree(gt_points)  # Create KDTree structure

    # Initialize result containers
    closest_points = []  # Stores closest ground truth points
    errors = []  # Stores error distances

    # Process each measured point
    for t_m, x_m, y_m, z_m in measured:
        # Query KDTree for nearest neighbor
        # Returns (distance, index) of closest point
        distance, idx = tree.query([x_m, y_m, z_m])

        # Store closest point information
        closest_points.append((
            t_m,  # Keep original timestamp
            x_interp[idx],  # X coordinate of closest point
            y_interp[idx],  # Y coordinate of closest point
            z_interp[idx]  # Z coordinate of closest point
        ))

        # Store error distance
        errors.append(distance)

    # Calculate mean error across all points
    mean_error = np.mean(errors)

    return closest_points, errors, mean_error


def calculate_speed_errors(points, closest_points):
    """
    Calculate speed errors between measured trajectory and ground truth trajectory.
    The first point's speed error is set to 0 instead of NaN.

    Args:
        points: List of measured points as (t, x, y, z) tuples
        closest_points: List of corresponding ground truth points as (t, x, y, z) tuples

    Returns:
        speed_errors: List of speed errors at each point (first point set to 0)
        mean_speed_error: Average speed error including the first point
    """

    # Convert to numpy arrays for vector operations
    # Extract timestamps and position coordinates separately
    t_meas = np.array([p[0] for p in points])  # Measurement timestamps
    pos_meas = np.array([p[1:] for p in points])  # Measurement positions (x,y,z)

    t_gt = np.array([p[0] for p in closest_points])  # Ground truth timestamps
    pos_gt = np.array([p[1:] for p in closest_points])  # Ground truth positions (x,y,z)

    # Calculate position differences between consecutive points
    # delta_pos_meas[i] = pos_meas[i+1] - pos_meas[i]
    delta_pos_meas = np.diff(pos_meas, axis=0)  # Position deltas for measurements
    delta_t_meas = np.diff(t_meas)  # Time deltas for measurements

    delta_pos_gt = np.diff(pos_gt, axis=0)  # Position deltas for ground truth
    delta_t_gt = np.diff(t_gt)  # Time deltas for ground truth

    # Calculate speed magnitude (norm of velocity vector)
    # speed = ||Δposition|| / Δt
    speed_meas = np.linalg.norm(delta_pos_meas, axis=1) / delta_t_meas  # Measured speeds
    speed_gt = np.linalg.norm(delta_pos_gt, axis=1) / delta_t_gt  # Ground truth speeds

    # Calculate absolute speed errors
    speed_errors = np.abs(speed_meas - speed_gt)

    # Insert 0 as error for first point (instead of 0)
    speed_errors = np.insert(speed_errors, 0, 0)  # Prepend 0 for first point

    # Calculate mean error (includes the 0 for first point)
    mean_speed_error = np.mean(speed_errors)

    return speed_errors, mean_speed_error


def calculate_heading_angle_errors(points, closest_points):
    """
    Calculate heading angle errors between measured and ground truth trajectories.
    Heading angle is computed as atan2(dy, dx) in the XY plane (ignoring Z-axis).

    Args:
        points: List of measured points as (t, x, y, z) tuples
        closest_points: List of corresponding ground truth points as (t, x, y, z) tuples

    Returns:
        heading_errors: List of heading angle errors in radians (-π to π)
        mean_heading_error: Mean absolute heading error
    """
    # Convert to numpy arrays
    t_meas = np.array([p[0] for p in points])
    xy_meas = np.array([p[1:3] for p in points])  # Extract only XY coordinates

    t_gt = np.array([p[0] for p in closest_points])
    xy_gt = np.array([p[1:3] for p in closest_points])

    # Calculate XY displacements
    dxy_meas = np.diff(xy_meas, axis=0)
    dt_meas = np.diff(t_meas)

    dxy_gt = np.diff(xy_gt, axis=0)
    dt_gt = np.diff(t_gt)

    # Normalize by time to get XY velocities
    vxy_meas = dxy_meas / dt_meas[:, np.newaxis]
    vxy_gt = dxy_gt / dt_gt[:, np.newaxis]

    # Calculate heading angles using atan2(dy, dx)
    heading_meas = np.arctan2(vxy_meas[:, 1], vxy_meas[:, 0])  # [-π, π]
    heading_gt = np.arctan2(vxy_gt[:, 1], vxy_gt[:, 0])

    # Calculate angular errors (handling circular nature)
    heading_errors = np.arctan2(np.sin(heading_meas - heading_gt),
                                np.cos(heading_meas - heading_gt))

    # Insert 0 for first point and convert to degrees
    heading_errors = np.insert(np.abs(heading_errors), 0, 0)
    heading_errors_deg = np.rad2deg(heading_errors)

    mean_error = np.mean(heading_errors_deg)

    return heading_errors_deg, mean_error

def FJ_target(input_trajectory,real_trajectory):
    closest_points,errors,point_error = find_closest_points(real_trajectory, input_trajectory)
    heading_errors_deg, mean_deg_error = calculate_heading_angle_errors(closest_points, input_trajectory)
    speed_errors, mean_speed_error = calculate_speed_errors(closest_points, input_trajectory)
    return point_error, mean_speed_error,mean_deg_error

def DD_target(input_trajectory,real_trajectory):
    closest_points, errors, position_mean_error = find_closest_points(real_trajectory, input_trajectory)
    speed_errors, mean_speed_error = calculate_speed_errors(closest_points, input_trajectory)
    return position_mean_error, mean_speed_error


def map_values(input_str, mapping):

    # 清理输入
    input_str = input_str.strip()

    # 分割输入
    if ',' in input_str:
        numbers = [num.strip() for num in input_str.split(',')]
    else:
        numbers = [input_str]

    # 映射值
    results = []
    for num in numbers:
        if num in mapping:
            results.append(mapping[num])
        else:
            results.append(f"未知值{num}")

    return results



def run(flag,choice):
    truth = read_json(r"D:\pythoncode\803\test\target6.json")
    measured = read_json(r"D:\pythoncode\803\test\target6.json")
    if flag == 0:
        point_error, mean_speed_error, mean_deg_error = FJ_target(measured, truth)
        mapping = {
            '1': point_error,
            '2': mean_speed_error,
            '3': mean_deg_error
        }
        results = map_values(choice, mapping)
        print(results)
    if flag == 1:
        position_mean_error, mean_speed_error = DD_target(measured, truth)
        mapping = {
            '1': position_mean_error,
            '2': mean_speed_error,
        }
        results = map_values(choice, mapping)
        print(results)




run(0,"1")


    