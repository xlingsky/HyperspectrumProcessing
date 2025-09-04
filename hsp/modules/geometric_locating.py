import os
from hsp.rpcm import rpc_from_geotiff, compute_epsg
from hsp.utils import common
import hsp.sofa as sofa
import pyproj

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
        self.load_rpc(rpcpath)
        lon, lat = self.rpc.lon_offset, self.rpc.lat_offset
        self.geographic_to_projection = pyproj.Transformer.from_crs(4326, compute_epsg(lon, lat), always_xy=True)
        self.geographic_to_geocentric = geographic_to_geocentric 

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