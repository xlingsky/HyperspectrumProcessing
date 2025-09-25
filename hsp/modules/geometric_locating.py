import os
from hsp.rpcm import rpc_from_geotiff, compute_epsg
from hsp.utils import common
import hsp.sofa as sofa
import pyproj
import numpy as np

GEOLOCATION_DIRNAME = 'geolocation'

geographic_to_geocentric = pyproj.Transformer.from_crs(
            "EPSG:4979",  # Geographic 3D CRS (WGS84 with height)
            "EPSG:4978",  # Geocentric CRS (WGS84)
            always_xy=True  # Ensure lon, lat order
        )

def transform_geographic_to_geocentric(lon, lat, alt):
    return geographic_to_geocentric.transform(lon, lat, alt)

def transform_geocentric_to_geographic(x, y, z): 
    return geographic_to_geocentric.transform(x, y, z, direction= pyproj.enums.TransformDirection.INVERSE)

def init(config : dict):
    eop = os.path.join(common.parent_dir, 'data', config['sofa'])
    if not os.path.exists(eop):
        config['sofa'] = None
        return False

    config['sofa'] = sofa.Transformer()

    if not config['sofa'].load_eop_data(eop):
        config['sofa'] = None
        return False

    return True

class Locating:
    def __init__(self, rpcpath):
        self.rpc = None
        if not self.load_rpc(rpcpath):
            return
        lon, lat = self.rpc.lon_offset, self.rpc.lat_offset
        self.geographic_to_projection = pyproj.Transformer.from_crs(4326, compute_epsg(lon, lat), always_xy=True)
        self.geographic_to_geocentric = geographic_to_geocentric 

    @property
    def good(self) -> bool:
        return self.rpc is not None

    def load_rpc(self, geotiff):
        try:
            self.rpc = rpc_from_geotiff(geotiff)
        except Exception as e:
            return False

        return self.rpc is not None

    def height_off(self):
        return self.rpc.alt_offset

    def height_scale(self):
        return self.rpc.alt_scale

    def transform(self, col, row, z=0):
        try:
            lon, lat = self.rpc.localization( col, row, z)
            x, y = self.geographic_to_projection.transform(lon, lat)
            X, Y, Z = self.geographic_to_geocentric.transform(lon, lat, z)
            return [lon, lat], [x, y], [X, Y, Z]
        except Exception as e:
            return [0,0], [0,0], [0,0,0]

def generate_trajectory_from_rpcs(rpcs : list, sample_number:int = 100, height:float = None, image_point:tuple = None):

    if len(rpcs) == 0:
        return {}

    if len(rpcs) < sample_number:
        idx = list(range(len(rpcs)))
    else:
        idx = np.linspace(0, len(rpcs)-1, sample_number, dtype=int)

    if height is None:
        height = rpcs[0].alt_offset

    if image_point is None:
        image_point = (rpcs[0].col_offset, rpcs[0].row_offset)

    return {i:rpcs[i].localization(image_point[0], image_point[1], height) for i in idx}

def generate_trajectory_from_geotiffs(geotiffs : list, sample_number:int = 100, height:float = None, image_point:tuple = None):

    if len(geotiffs) == 0:
        return {}

    rpcs = []
    for tiff in geotiffs:
        try:
            rpc = rpc_from_geotiff(tiff)
            if rpc is not None:
                rpcs.append(rpc)
        except Exception as e:
            continue

    lonlat = generate_trajectory_from_rpcs(rpcs, sample_number, height, image_point)

    return { geotiffs[k]:v for k,v in lonlat.items()}

def filter_2d_points(points:np.ndarray, window_size:int = 31, poly_order:int = 3):
    from scipy.signal import savgol_filter
    if window_size % 2 == 0:
        window_size -= 1
    x_filtered = savgol_filter( points[:, 0], window_size, poly_order)
    y_filtered = savgol_filter( points[:, 1], window_size, poly_order)
    return np.column_stack((x_filtered, y_filtered))

if __name__ == '__main__':
    import matplotlib.pyplot as plt
    import cv2

    dirpath = '/Users/xlingsky/Desktop/satellite_data_tiles/mountain/GF7'
    time_interval = 90
    max_speed = None
    show_original = False

    files = []
    for root, dirs, filenames in os.walk(dirpath):
        for filename in filenames:
            if filename.lower().endswith('.tif') or filename.lower().endswith('.tiff'):
                files.append(os.path.join(root, filename))
    files.sort()
    img_lonlat = generate_trajectory_from_geotiffs(files, len(files)//time_interval)

    lonlat = np.array( list(img_lonlat.values()) )

    geographic_to_projection = pyproj.Transformer.from_crs(4326, compute_epsg(lonlat[0,0], lonlat[0,1]), always_xy=True)

    proj = np.array( geographic_to_projection.transform(lonlat[:,0], lonlat[:,1]) ).T

    proj_offset = proj[0,:]
    proj -= proj_offset

    proj_filtered = filter_2d_points( proj, window_size=min(31,proj.shape[0]) )

    v = np.gradient(proj_filtered, axis=0)/time_interval
    vn = np.linalg.norm(v, axis=1)

    if max_speed is not None:
        flags = vn > max_speed
    else:
        thresh, flags  = cv2.threshold(vn.astype(np.uint16), 0, 1, cv2.THRESH_BINARY+cv2.THRESH_OTSU)
        flags = flags.astype(np.bool).flatten()

    plt.figure()
    if show_original:
        plt.plot(proj[:, 0], proj[:, 1], 'b-', alpha=0.7, label='Clean')
    plt.plot(proj_filtered[:, 0], proj_filtered[:, 1], 'g-', linewidth=2, label='SG Filtered')
    plt.scatter(proj_filtered[:, 0][flags], proj_filtered[:, 1][flags], color='red', s=20)
    plt.scatter(proj_filtered[:, 0][~flags], proj_filtered[:, 1][~flags], color='green', s=20)
    plt.title('Comparison: Clean vs Filtered')
    plt.xlabel('X')
    plt.ylabel('Y')
    plt.axis('equal')
    plt.grid(True, alpha=0.3)
    plt.legend()
    
    plt.tight_layout()
    plt.show()

    exit(0)