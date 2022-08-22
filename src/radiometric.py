import numpy as np
from osgeo import gdal
import os
import sys

def getFilePath(dirpath, filename):
    if os.path.exists(filename):
        return filename
    else:
        path = os.path.join(dirpath, filename)
        if os.path.exists(path):
            return path
    return None

def detectBadpixels(light_mean, light_std, dark_mean, dark_std):
    factors = [1.5, 0.5, 2, 1.5, 0.5]
    badpixels = list()
    module = [[0,512], [512,1024], [1024,1536], [1536,2048]]
    for m in module:
        lm = light_mean[:,m[0]:m[1]]
        ls = light_std[:,m[0]:m[1]]
        dm = dark_mean[:,m[0]:m[1]]
        ds = dark_std[:,m[0]:m[1]]
        dd = lm-dm
        ddm = np.mean(dd, axis=1)
        lsm = np.mean(ls, axis=1)
        dmm = np.mean(dm, axis=1)
        rows, cols = dd.shape
        for r in range(rows):
            for c in range(cols):
                if dd[r,c]>factors[0]*ddm[r] or dd[r,c]<factors[1]*ddm[r] or ls[r,c]>factors[2]*lsm[r] or dm[r,c]>factors[3]*dmm[r] or dm[r,c]<factors[4]*dmm[r]:
                    badpixels.append([r+1, c+m[0]+1])
    return badpixels

def computeNonuniform(high_mean, low_mean, high_dn = None, low_dn = None):
    rows, cols = high_mean.shape
    a = np.zeros_like(high_mean)
    b = np.zeros_like(high_mean)
    if high_dn is None:
        high_dn = np.median(high_mean, axis=1)
    if low_dn is None:
        low_dn = np.median(low_mean, axis=1)

    diff1 = high_dn-low_dn
    diff2 = high_mean-low_mean

    diff2[diff2==0] = 0.0001

    for r in range(rows):
        try:
            a[r,:] = (diff1[r]/diff2[r,:]) # .reshape((rows,1))
            b[r,:] = low_dn[r]-a[r,:]*low_mean[r,:]
        except:
            print("{}".format(r))

    return a, b

def moduleShift(imagepath, shift=577, yr=[0,2230]):
    src = gdal.Open(imagepath)
    if src is None:
        return None

    margin = 0
    yr[0] = 0
    yr[1] = src.RasterYSize-shift

    outpath = os.path.splitext(imagepath)[0]+"_shift.tif"
    driver = gdal.GetDriverByName('GTiff')
    dst = driver.Create(outpath, src.RasterXSize-3*margin, yr[1], src.RasterCount, src.GetRasterBand(1).DataType)
    for i in range(1, src.RasterCount+1):
        data = src.GetRasterBand(i).ReadAsArray()
        ye = shift+yr[0]+yr[1]
        if ye > data.shape[0]:
            ye = data.shape[0]
        dst.GetRasterBand(i).WriteRaster(  0,  0, 512,yr[1], data[yr[0]:yr[1], 0:512])
        dst.GetRasterBand(i).WriteRaster(512,  0, 512-margin,yr[1], data[shift+yr[0]:ye,512+margin:1024])
        dst.GetRasterBand(i).WriteRaster(1024-margin, 0, 512-margin,yr[1], data[yr[0]:yr[1], 1024+margin:1536])
        # dst.GetRasterBand(i).WriteRaster(1536-2*margin, 0, 512,yr[1], data[yr[0]:yr[1],1536+margin:])

    dst = None
    src = None

def vnir_badpixel(filepath):
    xrng = range(1461, 1464)
    brng = range(156, 206)
    badpixels = list()
    for b in brng:
        for x in xrng:
            badpixels.append([b, x])
    np.savetxt(filepath, np.asarray(badpixels), fmt='%d')

def simulateGeo(imagepath):
    wkt_wgs84 = "GEOGCS[\"WGS 84\",DATUM[\"WGS_1984\",SPHEROID[\"WGS 84\",6378137,298.257223563,AUTHORITY[\"EPSG\",\"7030\"]],AUTHORITY[\"EPSG\",\"6326\"]],PRIMEM[\"Greenwich\",0,AUTHORITY[\"EPSG\",\"8901\"]],UNIT[\"degree\",0.0174532925199433,AUTHORITY[\"EPSG\",\"9122\"]],AXIS[\"Latitude\",NORTH],AXIS[\"Longitude\",EAST],AUTHORITY[\"EPSG\",\"4326\"]]"
    xmin = 120
    ymax = 34
    gsd = 1/110000
    img = gdal.Open(imagepath, gdal.GA_Update)
    img.SetGeoTransform([xmin, gsd, 0, ymax, 0, -gsd])
    img.SetProjection(wkt_wgs84)

    geopath = os.path.splitext(imagepath)[0]+"_xyz.dat"
    geo = gdal.GetDriverByName('ENVI').Create(geopath, img.RasterXSize, img.RasterYSize, 3, gdal.GDT_Float64)
    geo.SetGeoTransform([xmin, gsd, 0, ymax, 0, -gsd])
    geo.SetProjection(wkt_wgs84)

    x = np.zeros((img.RasterYSize, img.RasterXSize))
    geo.GetRasterBand(3).WriteArray(x)
    for i in range(img.RasterXSize):
        x[:, i] = xmin+i*gsd
    geo.GetRasterBand(1).WriteArray(x)
    for i in range(img.RasterYSize):
        x[i, :] = ymax-i*gsd
    geo.GetRasterBand(2).WriteArray(x)
    geo = None


    img = None

def main():
    # light_mean = np.loadtxt('/Volumes/data/863/lab/swir/E4/8.8/budark/20170525-172007-SWIR推扫-直通-全谱段-存储缓存99.3/mean.txt')
    # light_std  = np.loadtxt('/Volumes/data/863/lab/swir/E4/8.8/budark/20170525-172007-SWIR推扫-直通-全谱段-存储缓存99.3/mean_std.txt')
    # dark_mean  = np.loadtxt('/Volumes/data/863/lab/swir/E4/8.8/dark/20170525-172154-SWIR直通-全谱段-丢行帧数5,帧号0xa7c1e-3084丢1次-图像帧丢4次-多数4次共5848字节-存储缓存100.0/mean.txt')
    # dark_std   = np.loadtxt('/Volumes/data/863/lab/swir/E4/8.8/dark/20170525-172154-SWIR直通-全谱段-丢行帧数5,帧号0xa7c1e-3084丢1次-图像帧丢4次-多数4次共5848字节-存储缓存100.0/mean_std.txt')

    # badpixels = detectBadpixels(light_mean, light_std, dark_mean, dark_std)
    # np.savetxt('/Volumes/data/863/lab/swir/badpixels.txt', np.asarray(badpixels))

    # return None
    highpath = '/Volumes/data/20220525/863FEIXIN/E7/8.8/budark/20170525-173647-SWIR推扫-直通-全谱段/00000000_mod_mean.txt'
    lowpath = '/Volumes/data/20220525/863FEIXIN/E2/8.8/budark/20170525-170227-SWIR推扫-直通-全谱段/00000000_mod_mean.txt'

    high_mean = np.loadtxt(highpath)
    low_mean = np.loadtxt(lowpath)

    a,b = computeNonuniform(high_mean, low_mean)

    np.savetxt('/Volumes/data/20220525/863FEIXIN/nonuniform_a.txt', a)
    np.savetxt('/Volumes/data/20220525/863FEIXIN/nonuniform_b.txt', b)
    
if __name__ == "__main__":
   moduleShift(sys.argv[1]);
