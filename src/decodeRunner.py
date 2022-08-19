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

    def vnir(self, strips, dstdir = r'H:\jinchang\ws\vnir'):
        config_radiometric = r'D:\jinchang\config\vnir\radiometric.xml'
        dir_radiometric = os.path.join(dstdir, '1_radio')
        xrange= [300, 1649]
        if not os.path.exists(dir_radiometric):
            os.mkdir(dir_radiometric)
        tsk_radiometric = os.path.join(dstdir, '1_radio.bat')
        with open(tsk_radiometric,'w') as f:
            for strip in strips:
                radiolist = list()
                for pos in strip[1]:
                    img = os.path.splitext(pos)[0]+'.dat'
                    name = strip[0]+'_'+os.path.splitext(os.path.basename(img))[0]
                    dstpath = os.path.join(dir_radiometric, name)+'.tif'
                    dataset = gdal.Open(img)
                    if not dataset:
                        print('NOT FOUND {}\n'.format(img))
                        continue
                    f.write('{} -task {} {} -o {}\n'.format(self._decode, config_radiometric, img, dstpath))
                    radiolist.append([dstpath, dataset.RasterXSize, dataset.RasterYSize])
                strip.append(radiolist)
        dir_cut = os.path.join(dstdir, '2_cut')
        if not os.path.exists(dir_cut):
            os.mkdir(dir_cut)
        tsk_cut = os.path.join(dstdir,'2_cut.bat')
        with open(tsk_cut,'w') as f:
            hdr = 'gdal_translate -srcwin {} {} {} '.format(xrange[0], 0, xrange[1]-xrange[0]+1)
            for strip in strips:
                cutlist = list()
                for img in strip[2]:
                    dstpath = os.path.join(dir_cut, os.path.basename(img[0]))
                    f.write('{} {} {} {}\n'.format(hdr, img[2], img[0], dstpath))
                    cutlist.append(dstpath)
                strip.append(cutlist)
        with open(os.path.join(dstdir,'2_poscopy.bat'), 'w') as f:
            for strip in strips:
                for pos in strip[1]:
                    name = strip[0]+'_'+os.path.basename(pos)
                    f.write('cp {} {}\n'.format(pos, os.path.join(dir_cut, name)))
        dir_strip = os.path.join(dstdir, '3_strip')
        if not os.path.exists(dir_strip):
            os.mkdir(dir_strip)
        dir_strip_tsk = os.path.join(dir_strip, 'tsk')
        if not os.path.exists(dir_strip_tsk):
            os.mkdir(dir_strip_tsk)
        dir_strip_prod = os.path.join(dir_strip,'product')
        if not os.path.exists(dir_strip_prod):
            os.mkdir(dir_strip_prod)
        tsk_strip = os.path.join(dstdir, '3_strip.bat')
        with open(tsk_strip,'w') as f:
            for strip in strips:
                imglist = os.path.join(dir_strip_tsk,'{}.txt'.format(strip[0]))
                with open(imglist, 'w') as l:
                    for img in strip[3]:
                        l.write('{}\n'.format(img))
                fname = os.path.splitext(os.path.basename(strip[1][0]))[0]
                lname = os.path.splitext(os.path.basename(strip[1][-1]))[0]
                f.write('{} {} -o {}.tif\n'.format(self._decode,imglist, os.path.join(dir_strip_prod, strip[0]+'_'+fname+'_'+lname)))


    def swir(self, strips, dstdir = r'D:\jinchang\ws\swir'):
        config_radiometric = r'D:\jinchang\config\swir\radiometric.xml'
        dir_radiometric = os.path.join(dstdir, '1_radio')
        xrange= [0, 1535]
        if not os.path.exists(dir_radiometric):
            os.mkdir(dir_radiometric)
        tsk_radiometric = os.path.join(dstdir, '1_radio.bat')
        with open(tsk_radiometric,'w') as f:
            for strip in strips:
                radiolist = list()
                for pos in strip[1]:
                    img = os.path.splitext(pos)[0]+'.dat'
                    name = strip[0]+'_'+os.path.splitext(os.path.basename(img))[0]
                    dstpath = os.path.join(dir_radiometric, name)+'.tif'
                    dataset = gdal.Open(img)
                    if not dataset:
                        print('NOT FOUND {}\n'.format(img))
                        continue
                    f.write('{} -task {} {} -o {}\n'.format(self._decode, config_radiometric, img, dstpath))
                    radiolist.append([dstpath, dataset.RasterXSize, dataset.RasterYSize])
                strip.append(radiolist)
        dir_cut = os.path.join(dstdir, '2_cut')
        if not os.path.exists(dir_cut):
            os.mkdir(dir_cut)
        tsk_cut = os.path.join(dstdir,'2_cut.bat')
        with open(tsk_cut,'w') as f:
            hdr = 'gdal_translate -srcwin {} {} {} '.format(xrange[0], 0, xrange[1]-xrange[0]+1)
            for strip in strips:
                cutlist = list()
                for img in strip[2]:
                    dstpath = os.path.join(dir_cut, os.path.basename(img[0]))
                    f.write('{} {} {} {}\n'.format(hdr, img[2], img[0], dstpath))
                    cutlist.append(dstpath)
                strip.append(cutlist)
        with open(os.path.join(dstdir,'2_poscopy.bat'), 'w') as f:
            for strip in strips:
                for pos in strip[1]:
                    name = strip[0]+'_'+os.path.basename(pos)
                    f.write('cp {} {}\n'.format(pos, os.path.join(dir_cut, name)))
        # dir_strip = os.path.join(dstdir, '3_strip')
        # if not os.path.exists(dir_strip):
        #     os.mkdir(dir_strip)
        # dir_strip_tsk = os.path.join(dir_strip, 'tsk')
        # if not os.path.exists(dir_strip_tsk):
        #     os.mkdir(dir_strip_tsk)
        # dir_strip_prod = os.path.join(dir_strip,'product')
        # if not os.path.exists(dir_strip_prod):
        #     os.mkdir(dir_strip_prod)
        # tsk_strip = os.path.join(dstdir, '3_strip.bat')
        # with open(tsk_strip,'w') as f:
        #     for strip in strips:
        #         imglist = os.path.join(dir_strip_tsk,'{}.txt'.format(strip[0]))
        #         with open(imglist, 'w') as l:
        #             for img in strip[3]:
        #                 l.write('{}\n'.format(img))
        #         fname = os.path.splitext(os.path.basename(strip[1][0]))[0]
        #         lname = os.path.splitext(os.path.basename(strip[1][-1]))[0]
        #         f.write('{} {} -o {}\n'.format(self._decode,imglist, os.path.join(dir_strip_prod, strip[0]+'_'+fname+'_'+lname)))

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
