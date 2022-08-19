from pos import POS
import os
from shapely.geometry import Point, Polygon

def load(filepath):
    ret = list()
    with open(filepath,'r') as f:
        for line in f:
            pospath = line.rstrip()
            if len(pospath)<3:
                continue
            if not os.path.exists(pospath):
                print("warning: NOT FOUND {}".format(pospath))
            p = POS(pospath)
            ret.append([pospath,p])
    return ret

def search(poslist, poly):
    ret = list()
    for p in poslist:
        record = p[1].search(poly)
        if len(record)>0:
            ret.append([p[0],record])
    return ret

def searchRect(poslist, xmin, xmax, ymin, ymax):
    poly = Polygon([(xmin, ymin),(xmin, ymax),(xmax, ymax),(xmax, ymin)])
    return search(poslist, poly)

def rerange(searchret):
    ret = list()
    last_id = -1
    group = list()
    for sr in searchret:
        name = os.path.splitext(os.path.basename(sr[0]))[0]
        id = int(name,16)
        if abs(id - last_id) != 1 and len(group)>0:
            ret.append(group)
            group = list()
        group.append(sr)
        last_id = id
    return ret

def createTask(rr, dirpath):
    ext = '.dat'
    ret = list()
    for r in rr:
        n1 = os.path.splitext(os.path.basename(r[0][0]))[0]
        n2 = os.path.splitext(os.path.basename(r[-1][0]))[0]
        taskpath = os.path.join(dirpath,'{}_{}.txt'.format(n1,n2))
        with open(taskpath,'w') as f:
            for i in r:
                f.write('{}\n'.format(os.path.splitext(i[0])[0]+ext))
            ret.append(taskpath)
    return ret

def loadCSV(csvpath):
    with open(csvpath,'r') as f:
        ret = list()
        for l in f:
            items = l.rstrip().split(',')
            items[2] = float(items[2])
            items[3] = float(items[3])
            ret.append(items)
        return ret

def check(items, pt):
    dist = 0.006**2
    ret = list()
    for i in items:
        err = (i[2]-pt[0])**2+(i[3]-pt[1])**2
        if err < dist:
            ret.append(i)
    return ret

def main(filepath):
    xmin = 101.748
    xmax = 101.902
    ymin = 38.551
    ymax = 38.589

    margin = 0.006
    poslist = load(filepath)
    sr = searchRect(poslist, xmin-margin, xmax+margin, ymin-margin, ymax+margin)
    rr = rerange(sr)

    outfile = os.path.splitext(filepath)[0]+'_aoi.txt'
    with open(outfile,'w') as f:
        f.write('{}\n'.format(len(rr)))
        id = 1
        for g in rr:
            f.write('{} {}\n'.format(id, len(g)))
            id+=1
            for item in g:
                f.write('{}\n'.format(item[0]))

def loadStrips(filepath):
    with open(filepath,'r') as f:
        ret = list()
        ls = f.readlines()
        num_strip = int(ls[0].rstrip())
        lid = 1
        for i in range(num_strip):
            strip = list()
            info = ls[lid].split()
            sid = info[0]
            num_img = int(info[1])
            lid+=1
            strip.append(sid)
            strip.append([ls[lid+x].rstrip() for x in range(num_img)])
            lid+=num_img
            ret.append(strip)
        return ret
