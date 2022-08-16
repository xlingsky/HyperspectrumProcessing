import os

def getFilePath(dirpath, filename):
    if os.path.exists(filename):
        return filename
    else:
        path = os.path.join(dirpath, filename)
        if os.path.exists(path):
            return path
    return None

class decodeRunner:

    def __init__(self, filepath):
        load(filepath)
        self._decode = r''
        self._rar = r''
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

    # [rar] 0:none 1:rar 2:rar+rm
    def generateIndividualTask(self, taskfile, dstdir, rar=2, rarpath=''):
        with open(taskfile, 'w') as f:
            if rar>0 and len(rarpath)<3:
                rarpath = self._dirpath
            for d in self._data:
                filename = os.path.basename(d)
                f.write('{} {} -v=2 -o {}\n'.format(self._decode, d, os.path.join(dstdir, filename)))
                if rar>0:
                    f.write('{} {} {}\n'.format(self._rar, os.path.join(rarpath, os.path.splitext(filename)[0]+'.rar'), d))
                    if rar==2:
                        f.write('{} {}\n'.format(self._rm, d))
