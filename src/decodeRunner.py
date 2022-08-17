import os
from osgeo import gdal
import numpy as np

def getFilePath(dirpath, filename):
    if os.path.exists(filename):
        return filename
    else:
        path = os.path.join(dirpath, filename)
        if os.path.exists(path):
            return path
    return None

class decodeRunner:
    def __init__(self):
        self._decode = r'E:\xlingsky\repo\HyperspectrumProcessing\build\Release\decode.exe'
        self._rar = "\"C:\\Program Files\\WinRAR\\Rar.exe\""
        self._rm = 'rm '

    def load(self, filepath):
        with open(filepath) as f:
            dirpath = os.path.dirname(filepath)
            lines = f.readlines()
            data = list()
            for t in lines:
                filepath = getFilePath(dirpath, t.strip())
                if filepath:
                    data.append(filepath)
            if len(data) > 0:
                self._dirpath = dirpath
                self._data = data
                return True
        return False

    def splicing(self, srclist, taskfile, dstdir, ext='.dat'):
        with open(taskfile,'w') as f:
            for s in srclist:
                f.write('{} {} -o {} -v=2\n'.format(self._decode, s, os.path.join(dstdir, os.path.splitext(os.path.basename(s))[0]+ext)))
    
    def run(self, srclist, poslist, dstdir):
        xrange = [300,1650]
        cnt = len(srclist)
        for i in range(cnt):
            img = gdal.Open(srclist[i])
            rows = img.RasterYSize
            img = None
            cutpath = os.path.join(dstdir,os.path.splitext(os.path.basename(srclist[i]))[0]+'_cut.tif')
            cmd = ('gdal_translate -srcwin {} {} {} {} {} {}'.format(xrange[0],0,xrange[1]-xrange[0],rows, srclist[i], cutpath))
            print(cmd)
            os.system(cmd)
            cmd = ('{} {}'.format(self._decode, poslist[i]))
            print(cmd)
            os.system(cmd)
            gcpfile = os.path.splitext(poslist[i])[0]+".gcp"
            prjfile = os.path.splitext(poslist[i])[0]+".prj"
            arr = np.loadtxt(gcpfile)
            gcpimg = os.path.join(dstdir,os.path.splitext(os.path.basename(srclist[i]))[0]+'_gcp.tif')
            cmd = ('gdal_translate {} {} '.format(cutpath,gcpimg))
            for r in arr:
                cmd += ('-gcp {} {} {} {} '.format(r[0],r[1],r[2],r[3]))
            print(cmd)
            os.system(cmd)
            dstfile = os.path.join(dstdir,os.path.splitext(os.path.basename(srclist[i]))[0]+'_geo.tif')
            cmd = ('gdalwarp {} {}'.format(gcpimg, dstfile))
            print(cmd)
            os.system(cmd)

    # [rar] 0:none 1:rar 2:rar+rm
    def generateIndividualTask(self, taskfile, dstdir, rar=2, rarpath=''):
        with open(taskfile, 'w') as f:
            if rar>0 and len(rarpath)<3:
                rarpath = self._dirpath
            for d in self._data:
                if not os.path.exists(d):
                    continue
                filename = os.path.basename(d)
                f.write('{} {} -v=2 -o {}\n'.format(self._decode, d, os.path.join(dstdir, filename)))
                if rar>0:
                    rarfile = os.path.join(rarpath, os.path.splitext(filename)[0]+'.rar')
                    f.write('{} a {} {}\n'.format(self._rar, rarfile, d))
                    if rar==2:
                        f.write('{} {}\n'.format(self._rm, d))

def main():
    with open(r'L:\AOI\vnir\data.txt', 'r') as f:
        srclist = list()
        poslist = list()
        for line in f:
            d = line.split()
            srclist.append(d[0].rstrip())
            poslist.append(d[1].rstrip())
        runner = decodeRunner()
        runner.run(srclist, poslist, r'L:\AOI\vnir\geo')
