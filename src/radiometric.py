import numpy as np
from osgeo import gdal
import os
from pykml.factory import KML_ElementMaker as KML
from pykml.factory import GX_ElementMaker as GX
import lxml.etree as etree

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

def moduleShift(imagepath, shift=548, yr=[0,2230]):
    src = gdal.Open(imagepath)
    if src is None:
        return None

    margin = 10

    outpath = os.path.splitext(imagepath)[0]+"_shift.tif"
    driver = gdal.GetDriverByName('GTiff')
    dst = driver.Create(outpath, src.RasterXSize-3*margin, yr[1], src.RasterCount, src.GetRasterBand(1).DataType)
    for i in range(1, src.RasterCount+1):
        data = src.GetRasterBand(i).ReadAsArray()
        ye = shift+yr[0]+yr[1]
        if ye > data.shape[0]:
            ye = data.shape[0]
        dst.GetRasterBand(i).WriteRaster(  0,  0, 512,yr[1], data[shift+yr[0]:ye, 0:512])
        dst.GetRasterBand(i).WriteRaster(512,  0, 512-margin,yr[1], data[yr[0]:yr[1],512+margin:1024])
        dst.GetRasterBand(i).WriteRaster(1024-margin, 0, 512-margin,yr[1], data[shift+yr[0]:ye, 1024+margin:1536])
        dst.GetRasterBand(i).WriteRaster(1536-2*margin, 0, 512,yr[1], data[yr[0]:yr[1],1536+margin:])

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

def dms2dd(s):
    d = s.split(':')
    return float(d[0])+float(d[1])/60+float(d[2])/3600

def loadGCP(gcppath):
    ret = list()
    with open(gcppath, 'r') as f:
        lines = f.readlines()
        for l in lines:
            d = l.split()
            name = d[0]
            lon = dms2dd(d[5])
            lat = dms2dd(d[4])
            height = float(d[6])
            ret.append([name, lon, lat, height])
    return ret

def saveKML(kmlpath, gcps):
    fold = KML.Folder(KML.Placemark())

    for g in gcps:
        fold.append(KML.Placemark(
            KML.name(g[0]),
            KML.Point(KML.coordinates(str(g[1]) +','+ str(g[2]) +','+str(g[3]))))
                    )
    content = etree.tostring(etree.ElementTree(fold), pretty_print=True)

    with open(kmlpath, 'w') as fp:
        fp.write(content.decode('utf-8'))

def MARKPOSA(s):
    pos = s.find('SINGLE')
    if pos == -1:
        return 0,0,0
    t = s[pos:].split(',')
    return float(t[2]),float(t[1]),float(t[3])

def MARK1PVAA(s):
    pos = s.find(';')
    if pos == -1:
        return 0,0,0
    t = s[pos:].split(',')
    return float(t[3]),float(t[2]),float(t[4])

class POS:
    def __init__(self, filepath):
        self.load(filepath)
    def load(self, filepath):
        pos = list()
        keyword = 'MARK1PVAA'
        with open(filepath,'r') as f:
            self._name = os.path.splitext(os.path.basename(filepath))[0]
            lines = f.readlines()
            last = ''
            for i, l in enumerate(lines):
                p = l.find(keyword)
                if p==-1:
                    continue
                t = l[p:]
                if last != t:
                    last = t
                    lon,lat,hei = MARK1PVAA(t)
                    pos.append([i, lon,lat,hei])
        self._pos = pos
        return len(pos)>0
    def toPoints(self, num=2):
        fold = KML.Folder()
        if num == 1:
            interval = len(self._pos)
        else:
            interval = (len(self._pos)-1)//(num-1)
        for i in range(0,len(self._pos),interval):
            p = self._pos[i]
            fold.append(KML.Placemark(
                KML.name('[{}]{}'.format(self._name,p[0])),
                KML.styleUrl("#point"),
                KML.Point(
                    GX.altitudeMode("relativeToSeaFloor"),
                    KML.coordinates(str(p[1]) +','+ str(p[2]) +','+str(p[3])))
            ))
        return fold
    def toLineString(self):
        coord = ''
        for p in self._pos:
            coord += ' '+str(p[1]) +','+ str(p[2]) +','+str(p[3])
        fold = KML.Folder(KML.Placemark(
            KML.name(self._name),
            KML.styleUrl("#line"),
            KML.LineString(
                GX.altitudeMode("relativeToSeaFloor"),
                KML.coordinates(coord)
                    )
        ))
        return fold
    def toKML(self, kmlpath, type='line'):
        if type == 'line':
            fold = self.toLineString()
            style = KML.Style(
                KML.LineStyle(
                    KML.color('7f0000ff'),
                    KML.width('3'),
                ),
                id='line'
            )
        else:
            fold = self.toPoints()
            style = KML.Style(
                KML.IconStyle(
                    KML.color('7f0000ff'),
                    KML.scale(1.0),
                ),
                id='point'
            )
        doc = KML.Document()
        doc.append(fold)
        doc.append(style)

        content = etree.tostring(etree.ElementTree(doc), pretty_print=True)

        with open(kmlpath, 'w') as fp:
            fp.write(content.decode('utf-8'))

def PosFiles2KML(poslist, kmlpath, name='real flight'):
    doc = KML.Document()
    doc.append(KML.name(name))
    style = KML.Style(
        KML.IconStyle(
            KML.color('7f0000ff'),
            KML.scale(1.0),
        ),
        id='point'
    )
    doc.append(style)
    for p in poslist:
        pos = POS(p)
        info = pos.toPoints()
        doc.append(info)

    content = etree.tostring(etree.ElementTree(doc), pretty_print=True)

    with open(kmlpath, 'w') as fp:
        fp.write(content.decode('utf-8'))

def PosDir2KML(dirpath, kmlpath):
    name = os.path.basename(dirpath)
    poslist = list()
    for root, _, files in os.walk(os.path.abspath(dirpath)):
        for f in files:
            if f.find(".pos")!=-1:
                poslist.append(os.path.join(root, f))
    PosFiles2KML(poslist, kmlpath)

def test(stripfile, kmlpath):
    h = "3800"
    strips = list()
    with open(stripfile, 'r') as f:
        lines = f.readlines()
        for t in lines:
            s = t.split()
            if len(s)<5:
                continue
            info = dict()
            info['name'] = s[0]
            info['coord'] = s[1]+','+s[2]+','+h+' '+s[3]+','+s[4]+','+h
            info['xyz'] = [float(s[1]),float(s[2]),float(h),float(s[3]),float(s[4]),float(h)]
            strips.append(info)
    if len(strips)<1:
        return False

    doc = KML.Document()
    style = KML.Style(
        KML.LineStyle(
            KML.color('007f00ff'),
            KML.width('3'),
        ),
        id='line'
    )
    doc.append(KML.name("designed flight"))
    fold = KML.Folder()
    for s in strips:
        fold.append(
            KML.Placemark(
            KML.name(s['name']),
            KML.styleUrl("#line"),
            KML.MultiGeometry(
            KML.LineString(
                GX.altitudeMode("relativeToSeaFloor"),
                KML.coordinates(s['coord'])
                    ),
            KML.Point(
                GX.altitudeMode("relativeToSeaFloor"),
                KML.coordinates(str((s['xyz'][0]+s['xyz'][3])/2)+','+str((s['xyz'][1]+s['xyz'][4])/2)+','+str((s['xyz'][2]+s['xyz'][5])/2))
            )
            )
        ))
    doc.append(fold)
    content = etree.tostring(etree.ElementTree(doc), pretty_print=True)

    with open(kmlpath, 'w') as fp:
        fp.write(content.decode('utf-8'))

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
