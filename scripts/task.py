from enum import Enum
import multiprocessing 
import time
from multiprocessing import Event
import websockets
import json

class Task:
    class Status(Enum):
        WAIT = 0
        START = 1
        RUNNING = 2
        PAUSE = 3
        FINISH = 4
        TERMINATED = 5
        DIE = 6

    def __init__(self):
       self._process = None
       self._status = Task.Status.WAIT
       self._pause_event = Event()
       self._stop_event = Event()
    def resume(self):
        if self._status != Task.Status.PAUSE:
            return False
        self._pause_event.clear()
        return True
    def start(self, target, jobid, taskid):
        if self._status == Task.Status.RUNNING:
            return False
        self._process = multiprocessing.Process(target=self._worker, args=(jobid, taskid, self._stop_event, self._pause_event, self._status), kwargs={})
        self._process.start()
        time.sleep(1)
        return self._status == Task.Status.RUNNING
    def pause(self):
        if self._status != Task.Status.RUNNING:
            return False
        self._pause_event.set()
        return True
    def terminate(self):
        if self._status != Task.Status.RUNNING:
            return False
        self._stop_event.set()
        return True
    def finished(self):
        return self._status == Task.Status.FINISH
    def check(self):
        if self._process.is_alive():
            self._status = Task.Status.RUNNING
        else:
            self._status = Task.Status.FINISH
    def status(self):
        return self._status
    def _worker(self, jobid, taskid, stop_event, pause_event, status):
        websocket = websockets.connect("")
        status = Task.Status.START
        message = {
            "jobID": jobid,
            "taskID": taskid,
            "total":100,
            "progress": 0
        }
        for i in range(100):
            status = Task.Status.RUNNING

            if stop_event.is_set():
                status = Task.Status.TERMINATED
                return 1
            while pause_event.is_set():
                time.sleep(0.5)
                status = Task.Status.PAUSE
                if stop_event.is_set():
                    status = Task.Status.TERMINATED
                    return 1
            time.sleep(10)

            message['progress'] = i+1
            websocket.send(json.dumps(message))
        status = Task.Status.FINISH
        

    