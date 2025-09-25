import numpy as np
from scipy.signal import correlate
from datetime import datetime, timedelta 
from hsp.modules.geometric_locating import transform_geocentric_to_geographic, transform_geographic_to_geocentric
import json

TRAJECTORY_PREDICTION_DIRNAME = 'prediction'

# 地球相关常数
GM = 3.986004418e14  # 地球引力参数 m³/s²
g = 9.81  # 重力加速度 m/s²

def detect_acceleration_jump(accelerations, threshold_ratio=0.3, min_jump=10.0):
    """
    检测加速度曲线中的突变，判断是否为二级助推

    Args:
        accelerations: 观测到的加速度数组
        threshold_ratio: 突变阈值比例（相对于平均加速度）
        min_jump: 最小突变值 (m/s²)

    Returns:
        tuple: (is_two_stage, jump_index) - 是否二级助推，突变位置索引
    """
    if len(accelerations) < 5:  # 数据点太少无法判断
        return False, -1

    # 计算加速度差分
    acc_diff = np.diff(accelerations)

    # 计算动态阈值
    avg_acc = np.mean(accelerations)
    dynamic_threshold = max(avg_acc * threshold_ratio, min_jump)

    # 寻找显著的正向跳跃（二级点火）
    jump_candidates = []
    for i in range(1, len(acc_diff) - 1):
        # 检查是否有显著的正向跳跃
        if acc_diff[i] > dynamic_threshold:
            # 验证跳跃的持续性（跳跃后保持较高水平）
            if i + 2 < len(accelerations):
                post_jump_avg = np.mean(accelerations[i + 1:min(i + 4, len(accelerations))])
                pre_jump_avg = np.mean(accelerations[max(0, i - 2):i + 1])

                if post_jump_avg - pre_jump_avg > dynamic_threshold * 0.5:
                    jump_candidates.append((i, acc_diff[i], post_jump_avg - pre_jump_avg))

    # 如果找到显著跳跃，选择最大的一个
    if jump_candidates:
        # 选择跳跃幅度最大的位置
        best_jump = max(jump_candidates, key=lambda x: x[2])
        return True, best_jump[0] + 1  # +1 因为diff索引比原数组小1

    return False, -1

def boost_acceleration_model_two_stage(t, params):
    """
    二级助推段加速度模型

    Args:
        t: 时间
        params: 参数字典 {
            'a0_1': 一级起始加速度,
            'a1_1': 一级结束加速度,
            'T1': 一级助推时间,
            'a0_2': 二级起始加速度,
            'a1_2': 二级结束加速度,
            'T2': 二级助推时间,
            'T_total': 总助推时间
        }

    Returns:
        float: 加速度值
    """
    T1 = params['T1']
    T_total = params['T_total']

    if t > T_total:
        return 0.0
    elif t <= T1:
        # 一级助推段
        a0_1 = params['a0_1']
        a1_1 = params['a1_1']
        return a0_1 + (a1_1 - a0_1) * (t / T1) ** 2
    else:
        # 二级助推段
        t_2 = t - T1
        T2 = params['T2']
        a0_2 = params['a0_2']
        a1_2 = params['a1_2']
        return a0_2 + (a1_2 - a0_2) * (t_2 / T2) ** 2

def boost_acceleration_model(t, a0=35.0, a1=70.0, T=20.0, stage_type='single', params=None):
    """
    统一的助推段加速度模型

    Args:
        t: 时间
        a0, a1, T: 单级助推参数（向后兼容）
        stage_type: 'single' 或 'two_stage'
        params: 二级助推参数字典

    Returns:
        float: 加速度值
    """
    if stage_type == 'two_stage' and params is not None:
        return boost_acceleration_model_two_stage(t, params)
    else:
        # 原有单级助推模型
        if t > T:
            return 0.0
        return a0 + (a1 - a0) * (t / T) ** 2

def up_direction_vector(x, y, z):
    lon, lat, alt = transform_geocentric_to_geographic(x, y, z)
    xu, yu, zu = transform_geographic_to_geocentric(lon, lat, alt + 100.0)
    up_vec = np.array([xu - x, yu - y, zu - z])
    return up_vec / np.linalg.norm(up_vec)

def find_boost_phase_position(observed_accelerations, model_params, stage_type='single', two_stage_params=None):
    """
    通过观测到的加速度曲线匹配模型，确定观测数据在助推段中的位置

    Args:
        observed_accelerations: 观测加速度
        model_params: 单级助推参数 (a0, a1, T)
        stage_type: 'single' 或 'two_stage'
        two_stage_params: 二级助推参数字典

    Returns:
        tuple: (best_match_start, best_correlation, best_mse, stage_info)
    """
    best_match_start = 0
    best_correlation = -1
    best_mse = float('inf')

    if stage_type == 'two_stage' and two_stage_params is not None:
        # 二级助推模型
        T_total = two_stage_params['T_total']
        model_times = np.linspace(0, T_total, int(T_total))
        model_acc = [boost_acceleration_model(t, stage_type='two_stage', params=two_stage_params)
                     for t in model_times]
        stage_info = {
            'type': 'two_stage',
            'T1': two_stage_params['T1'],
            'T2': two_stage_params['T2'],
            'T_total': T_total
        }
    else:
        # 单级助推模型
        a0, a1, T = model_params
        model_times = np.linspace(0, T, int(T))
        model_acc = [boost_acceleration_model(t, a0, a1, T) for t in model_times]
        stage_info = {
            'type': 'single',
            'T_total': T
        }

    obs_length = len(observed_accelerations)

    for start_idx in range(max(1, len(model_acc) - obs_length + 1)):
        end_idx = start_idx + obs_length
        if end_idx > len(model_acc):
            continue

        model_segment = model_acc[start_idx:end_idx]

        if len(model_segment) == obs_length:
            if np.std(model_segment) > 0 and np.std(observed_accelerations) > 0:
                observed_accelerations_norm = (observed_accelerations - np.mean(observed_accelerations)) / np.std(
                    observed_accelerations)
                model_segment_norm = (model_segment - np.mean(model_segment)) / np.std(model_segment)
                correlation = np.max(correlate(observed_accelerations_norm, model_segment_norm, mode='valid')) / len(
                    observed_accelerations)

                mse = np.mean((np.array(observed_accelerations) - np.array(model_segment)) ** 2)

                if not np.isnan(correlation) and mse < best_mse:
                    best_correlation = correlation
                    best_match_start = start_idx
                    best_mse = mse

    return best_match_start, best_correlation, best_mse, stage_info

def estimate_launch_point_kinematic(
        observed_positions_ecef, observed_velocities_ecef, match_start_time, stop_altitude, model_params, stage_type
        ):
    """
    增强版发射点估计，支持二级助推
    """
    dt = 1.0
    current_pos = observed_positions_ecef[0]
    current_vel = observed_velocities_ecef[0]

    pos = []

    for t_step in range(int(match_start_time)):
        t = match_start_time - t_step - 1

        if stage_type == 'two_stage':
            acc_magnitude = boost_acceleration_model(t, stage_type='two_stage', params=model_params)
        else:
            a0, a1, T = model_params
            acc_magnitude = boost_acceleration_model(t, a0, a1, T)

        if np.linalg.norm(current_vel) > 1.0:
            acc_direction = current_vel / np.linalg.norm(current_vel)
        else:
            acc_direction = up_direction_vector(*current_pos)

        acc_vector = acc_magnitude * acc_direction
        pos.append(current_pos-(current_vel * dt + 0.5 * acc_vector * dt ** 2))
        current_vel -=  acc_vector * dt
        current_pos = pos[-1]

        temp_blh = transform_geocentric_to_geographic(*current_pos)
        if temp_blh[2] <= stop_altitude:
            break

    return pos

def estimate_shutdown_point_kinematic(observed_positions_ecef, observed_velocities_ecef,
                                               observed_accelerations, match_start_time,
                                               model_params, stage_type):
    """增强版关机点推算，支持二级助推"""
    dt = 1.0

    obs_end_time = match_start_time + len(observed_accelerations) - 1

    current_pos = observed_positions_ecef[-1]
    current_vel = observed_velocities_ecef[-1]

    # 确定总助推时间
    if stage_type == 'two_stage':
        T_total = model_params['T_total']
    else:
        a0, a1, T = model_params
        T_total = T

    remaining_time = int(T_total - obs_end_time)
    if remaining_time <= 0:
        return [], current_vel

    pos = []
    for t_step in range(remaining_time):
        t = obs_end_time + t_step + 1

        if t >= T_total:
            break

        # 根据助推类型计算加速度
        if stage_type == 'two_stage':
            acc_magnitude = boost_acceleration_model(t, stage_type='two_stage', params=model_params)
        else:
            a0, a1, T = model_params
            acc_magnitude = boost_acceleration_model(t, a0, a1, T)

        if np.linalg.norm(current_vel) > 1.0:
            acc_direction = current_vel / np.linalg.norm(current_vel)
        else:
            acc_direction = up_direction_vector(*current_pos)

        acc_vector = acc_magnitude * acc_direction

        pos.append(current_pos+current_vel * dt + 0.5 * acc_vector * dt ** 2)
        current_pos = pos[-1] 
        current_vel += acc_vector * dt

    return pos, current_vel

def missile_impact_prediction( start_ecef, velocity_ecef, stop_altitude, time_step=0.1):
    """使用ECEF坐标系进行导弹落点预测"""

    predicted_points = []
    current_pos_ecef = np.array(start_ecef)
    current_vel_ecef = velocity_ecef.copy()

    point_counter = 1
    max_iterations = int(3600 / time_step)

    while point_counter < max_iterations:

        r = current_pos_ecef
        r_norm = np.linalg.norm(r)

        gravity_acc_ecef = -GM * r / (r_norm ** 3)

        current_vel_ecef += gravity_acc_ecef * time_step
        predicted_points.append(current_pos_ecef+current_vel_ecef * time_step)
        current_pos_ecef = predicted_points[-1]

        current_blh = transform_geocentric_to_geographic(current_pos_ecef[0], current_pos_ecef[1], current_pos_ecef[2])

        if current_blh[2] <= stop_altitude:
            break

        point_counter += 1

    return predicted_points
