import os
from pykml.factory import KML_ElementMaker as KML
from pykml.factory import GX_ElementMaker as GX
import lxml.etree as etree
from shapely.geometry import Point, Polygon
import numpy as np

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
        try:
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
        except:
            self._pos = []
            return False
    def search(self, poly):
        ret = list()
        for p in self._pos:
            pt = Point(p[1], p[2])
            if poly.contains(pt):
                ret.append(p)
        return ret

    def geotransform(self, col = 1024):
        p1 = self._pos[1]
        p2 = self._pos[-1]
        a = np.array([[col, -p1[0], 1, 0],[col, p1[0], 0, 1],[col, -p2[0], 1, 0],[col, p2[0], 0, 1] ])
        b = np.array([p1[1],p1[2],p2[1],p2[2]])
        x = np.linalg.solve(a,b)
        return x

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

def PosFiles2CSV(poslist, csvpath):
    with open(csvpath,'w') as f:
        for p in poslist:
            pos = POS(p)
            if len(pos._pos) < 1:
                print('{} NOT FOUND'.format(p))
                continue
            for i in pos._pos:
                f.write('{},{},{},{}\n'.format(p,i[0],i[1],i[2]))

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
