import numpy as np
from datetime import datetime
from bisect import bisect_right
import math
try:
    import hsp.sofa.ctypes.PyMsOfa as sofa
    SOFA_CTYPES = True
except:
    import hsp.sofa.python.PyMsOfa as sofa
    SOFA_CTYPES = False

# 闰秒表 (MJD, TAI-UTC)
LEAP_SECONDS = [
    (44239, 19),  # 1980-01-01
    (44786, 20),  # 1981-07-01
    (45151, 21),  # 1982-07-01
    (45516, 22),  # 1983-07-01
    (46247, 23),  # 1985-07-01
    (47161, 24),  # 1988-01-01
    (47892, 25),  # 1990-01-01
    (48257, 26),  # 1991-01-01
    (48804, 27),  # 1992-07-01
    (49169, 28),  # 1993-07-01
    (49534, 29),  # 1994-07-01
    (50083, 30),  # 1996-01-01
    (50630, 31),  # 1997-07-01
    (51179, 32),  # 1999-01-01
    (53736, 33),  # 2006-01-01
    (54832, 34),  # 2009-01-01
    (56109, 35),  # 2012-07-01
    (57204, 36),  # 2015-07-01
    (57754, 37),  # 2017-01-01
]

def utc2jdTime(year, month, day, hour, minute, second):
    """
    Convert UTC time to Julian Date.
    
    Args:
        year: Integer year
        month: Integer month (1-12)
        day: Integer day
        hour: Integer hour (0-23)
        minute: Integer minute (0-59)
        second: Double/float second
    
    Returns:
        Julian Date as float
    """
    # Handle January and February as months 13 and 14 of previous year
    y = year
    m = month
    if m <= 2:
        y -= 1
        m += 12

    A = y // 100
    B = 2 - A + (A // 4)

    # Calculate Julian Date
    jd = (math.floor(365.25 * (y + 4716)) + 
          math.floor(30.6001 * (m + 1)) + 
          day + B - 1524.5)

    # Add time of day
    timeFraction = (hour + minute / 60.0 + second / 3600.0) / 24.0
    jd += timeFraction

    return jd

def get_tai_utc_diff(mjd):
    """获取TAI-UTC差值"""
    if mjd < LEAP_SECONDS[0][0]:
        return 19

    mjds = [ls[0] for ls in LEAP_SECONDS]
    idx = bisect_right(mjds, mjd)

    if idx == 0:
        return LEAP_SECONDS[0][1]
    return LEAP_SECONDS[idx - 1][1]

def lagrange_interp(t, times, values):
    """拉格朗日插值"""
    result = 0.0
    n = len(times)
    for i in range(n):
        term = values[i]
        for j in range(n):
            if i != j:
                term *= (t - times[j]) / (times[i] - times[j])
        result += term
    return result


def interpolate_eop(target_mjd, eop_list):
    """插值EOP参数"""
    if not eop_list:
        # 如果没有EOP数据，返回默认值
        return {
            'xp': 0.0,
            'yp': 0.0,
            'dut': 0.0,
            'dat': get_tai_utc_diff(int(target_mjd))
        }

    # 找到合适的插值点
    mjds = [eop['mjd'] for eop in eop_list]
    k = bisect_right(mjds, target_mjd)
    k = max(1, min(k, len(eop_list) - 3))

    # 4点拉格朗日插值
    times = [float(eop_list[k - 1 + i]['mjd']) for i in range(4)]

    xp_values = [eop_list[k - 1 + i]['xp'] for i in range(4)]
    yp_values = [eop_list[k - 1 + i]['yp'] for i in range(4)]
    dut_values = [eop_list[k - 1 + i]['dut'] for i in range(4)]

    interp_eop = {
        'xp': lagrange_interp(target_mjd, times, xp_values),
        'yp': lagrange_interp(target_mjd, times, yp_values),
        'dut': lagrange_interp(target_mjd, times, dut_values),
        'dat': get_tai_utc_diff(int(target_mjd))
    }

    return interp_eop

class Transformer:
    def __init__(self, eop_list=None):
        self.eop_list = eop_list if eop_list is not None else []
    def load_eop_data(self, filename):
        eop_list = []
        try:
            with open(filename, 'r') as f:
                for line in f:
                    if line.strip() and 'I' in line:
                        # 去掉前7个字符
                        data = line[7:].split()
                        if len(data) >= 8:
                            eop = {
                                'mjd': float(data[0]),
                                'xp': float(data[2]),  # 极移X (角秒)
                                'yp': float(data[4]),  # 极移Y (角秒)
                                'dut': float(data[7])  # UT1-UTC (秒)
                            }
                            eop_list.append(eop)
        except FileNotFoundError:
            print(f"[WARNING]: EOP file '{filename}' not found. Using default values.")
            return False

        self.eop_list = eop_list
        return True

    def j2000_to_cgcs2000_matrix(self, time):
        if isinstance(time, str):
            time = datetime.strptime(time, "%Y-%m-%d %H:%M:%S")

        year, month, day, hour, minute, second = time.year, time.month, time.day, time.hour, time.minute, time.second

        eop_list = self.eop_list
        mjd = utc2jdTime(year, month, day, hour, minute, second) - 2400000.5

        # 插值EOP参数
        eop = interpolate_eop(mjd, eop_list)

        # 计算djmjd0和date
        if SOFA_CTYPES:
            djmjd0, date = sofa.pymCal2jd(year, month, day)
        else:
            djmjd0, date, _ = sofa.pymCal2jd(year, month, day)

        # 计算TT (地球时)
        # TT = UTC + DAT + 32.184
        TTMTAI = 32.184  # TT - TAI constant
        tt = date + (hour * 3600 + minute * 60 + second) / 86400.0 + (eop['dat'] + TTMTAI) / 86400.0

        # 计算TUT (UT1时间)
        tut = (hour * 3600 + minute * 60 + second + eop['dut']) / 86400.0

        # 1. 计算J2000到GCRS的转换矩阵
        # 使用pym 2006偏差-岁差模型
        rb, rp, rbp = sofa.pymBp06(djmjd0, tt)
        # 转置rb得到J2000到GCRS的矩阵
        rJ2000_to_GCRS = sofa.pymTr(rb)

        # 2. 计算GCRS到CIRS的转换矩阵（考虑岁差章动）
        x, y, s = sofa.pymXys06a(djmjd0, tt)
        rc2it = sofa.pymC2ixys(x, y, s)

        # 3. 计算CIRS到TIRS的转换矩阵（考虑地球自转）
        era = sofa.pymEra00(djmjd0 + date, tut)
        rc2ti = sofa.pymCr(rc2it)  # 复制矩阵
        rc2ti = sofa.pymRz(era, rc2ti)  # 绕Z轴旋转ERA角

        # 4. 计算TIRS到ITRS的转换矩阵（考虑极移）
        # 将角秒转换为弧度
        DAS2R = 4.848136811095359935899141e-6  # 角秒到弧度
        xp_rad = eop['xp'] * DAS2R
        yp_rad = eop['yp'] * DAS2R
        sp = sofa.pymSp00(djmjd0, tt)
        rpom = sofa.pymPom00(xp_rad, yp_rad, sp)

        # 组合CIRS到ITRS的转换
        rc2it_final = sofa.pymRxr(rpom, rc2ti)

        # 5. 计算从J2000到ITRS的完整转换矩阵
        rj2it = sofa.pymRxr(rc2it_final, rJ2000_to_GCRS)

        # CGCS2000采用ITRS框架，因此J2000到CGCS2000的转换矩阵就是rj2it
        return np.array(rj2it)
