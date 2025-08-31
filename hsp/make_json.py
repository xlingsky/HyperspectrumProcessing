import os
import pandas as pd
import numpy as np
import json
from datetime import datetime
from rpcm import rpc_from_geotiff, RPCModel, compute_epsg
from osgeo import osr

def make_json(dirpath, csv_path):
    current_time = datetime.now()
    # 定义文件夹路径
    folder_path = csv_path

    # 读取real_targets的csv
    real_targets_file_path = os.path.join(folder_path, 'real_targets.csv')

    real_targets_df = pd.read_csv(real_targets_file_path)

    # 读取real_targets.csv中的trajectory_id项
    trajectory_ids = real_targets_df['trajectory_id'].tolist()
    trajectory_distance = real_targets_df['total_distance'].tolist()

    # 定义时间原点
    time_origin = pd.Timestamp('2025-01-02 12:24:00')

    src = osr.SpatialReference()
    src.ImportFromEPSG(4326)
    dst = osr.SpatialReference()
    dst.ImportFromEPSG(4479)
    cvt = None

    # 遍历每个trajectory_id
    for idx,trajectory_id in enumerate(trajectory_ids):
        trajectory_file_path = os.path.join(folder_path, f'{trajectory_id}.csv')
        if os.path.exists(trajectory_file_path):
            # 读取轨迹信息
            trajectory_df = pd.read_csv(trajectory_file_path)

            for index, row in trajectory_df.iterrows():
                rpc = rpc_from_geotiff(os.path.join(dirpath, row['frame']))
                row['longitude'], row['latitude'] = rpc.localization(row['x'], row['y'], rpc.alt_offset)

            if cvt is None:
                dst.ImportFromEPSG(compute_epsg(trajectory_df['longitude'][0], trajectory_df['latitude'][0]))
                cvt = osr.CoordinateTransformation(src, dst)

            for index, row in trajectory_df.iterrows():
                row['x'], row['y'], _ = cvt.TransformPoint(row['latitude'], row['longitude'], rpc.alt_offset)

            # 移除frame列中的.tif后缀（如果存在）
            trajectory_df['frame'] = trajectory_df['frame'].str.replace('.tif', '', regex=False)

            # 将frame列转换为数值类型
            trajectory_df['frame'] = pd.to_numeric(trajectory_df['frame'], errors='coerce')

            # 过滤掉frame为NaN的行
            trajectory_df = trajectory_df.dropna(subset=['frame'])

            # 计算时间（秒）
            trajectory_df['time_seconds'] = (trajectory_df['frame'] - 10001).astype(float)

            # 计算位移（米）
            trajectory_df['dx'] = trajectory_df['x'].diff().fillna(0) 
            trajectory_df['dy'] = trajectory_df['y'].diff().fillna(0) 

            # 计算时间间隔（秒）
            trajectory_df['dt'] = trajectory_df['time_seconds'].diff().fillna(0)

            # 计算x和y方向的速度分量（m/s），处理dt为0的情况
            trajectory_df['vx'] = np.where(
                trajectory_df['dt'] > 0,
                trajectory_df['dx'] / trajectory_df['dt'],
                0
            )

            trajectory_df['vy'] = np.where(
                trajectory_df['dt'] > 0,
                trajectory_df['dy'] / trajectory_df['dt'],
                0
            )

            # 计算时间
            trajectory_df['time'] = time_origin + pd.to_timedelta(trajectory_df['frame'] - 10001, unit='s')

            # 初始化单条轨迹的JSON数据结构
            trajectory_data = {
                "Trajectory": {
                    "trajectory_id": trajectory_id,
                    "datatime":f"{current_time}",
                    "SatelliteList": {
                        "Satellite": [
                            {
                                "ID": "1"
                            }
                        ]
                    },
                    "total_distance":f"{trajectory_distance[idx]:.2f}",
                    "PointList": {
                        "Point": []
                    }
                }
            }

            # 将轨迹点添加到当前轨迹的数据中
            for index, row in trajectory_df.iterrows():
                if index == 0:
                    warning_status = "01H"
                elif index == len(trajectory_df) - 1:
                    warning_status = "03H"
                else:
                    warning_status = "02H"

                point = {
                    "Time": row['time'].strftime('%Y-%m-%d %H:%M:%S'),
                    "Location": f"{row['latitude']:.2f},{row['longitude']:.2f},{0.00:.2f}",
                    "Row and Column Numbers":f"{row['y']:.2f},{row['x']:.2f}",
                    "Velocity": f"{row['vx']:.2f},{row['vy']:.2f},{0.00:.2f}",
                    "gray_value":f"{row['gray_value']:.2f}",
                    "WarningStatus": warning_status,
                    "DigitalNumber": "100",
                    "Energy": "0.01",
                    "Coordinate": f"{row['x']:.2f},{row['y']:.2f},{0.00:.2f}"  # 添加(x,y,z)坐标信息
                }
                trajectory_data["Trajectory"]["PointList"]["Point"].append(point)

            # 保存单条轨迹数据到独立的JSON文件
            trajectory_json_path = os.path.join(folder_path, f'trajectory_{trajectory_id}.json')
            with open(trajectory_json_path, 'w') as f:
                json.dump(trajectory_data, f, indent=4)

            print(f"轨迹 {trajectory_id} 的JSON文件已生成")

if __name__ == '__main__':
    csv_path = r'D:\pythoncode\803\output\order\file\output_0'
    make_json(csv_path)