from xml.etree import ElementTree as et
import os

class task:
    def __init__(self) -> None:
        self.io_ = {'in':dict(), 'out':dict()}
        self.params_ = dict()
        self.flags_ = dict()
    def successed(self):
        for o in self.io_['out']:
            if 'path' in o and not os.path.exists(o['path']):
                return False
        return True
    def __str__(self):
        return ' '.join('{} \"{}\"'.format(k,v) for k,v in self.flags_.items())+' '.join(str(x['path']) for x in self.io_['in'])
    def input(self):
        return self.io_['in']
    def output(self):
        return self.io_['out']

    def export(self, tags, cfg):
        for t in tags:
            cfg[t] = self.io_['out'][t]
    def exportall(self, cfg):
        cfg.update(self.io_['out'])

    def initialize_statistic(self, cfg):
        self.io_['in']['img'] = cfg['img']
        self.params_['dim_order'] = cfg['dim_order']
        self.io_['out']['v_mean'] = {'path':, 'order':}
        self.io_['out']['v_std'] = {'path':, 'order':}
        self.flags_['-v'] = cfg['verbose']
        io to jsonfile

