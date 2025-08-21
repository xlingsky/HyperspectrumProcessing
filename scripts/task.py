from enum import Enum
import multiprocessing 
from multiprocessing import Event
import time
from websocketserver import WebSocketServer
import json

class Status(Enum):
    SUCCESS=200
    ERROR=500

class Task:
    def __init__(self):
       self._process = None
       self._pause_event = Event()
       self._stop_event = Event()
    def resume(self):
        self._pause_event.clear()
        return True
    def start(self, target, jobid, taskid, server):
        self._process = multiprocessing.Process(target=self._worker, args=(jobid, taskid, self._stop_event, self._pause_event, server), kwargs={})
        self._process.start()
        time.sleep(1)
        return True
    def pause(self):
        self._pause_event.set()
        return True
    def terminate(self):
        self._stop_event.set()
        return True

    def is_paused(self):
        return self._pause_event.is_set()
    def is_terminated(self):
        return self._stop_event.is_set()
    def is_alive(self):
        return self._process.is_alive() if self._process else False
    def is_running(self):
        return not (self.is_paused() or self.is_terminated()) if self.is_alive() else False
    def is_finished(self):
        return self._process.exitcode is not None if self._process else False
    def _worker(self, jobid, taskid, stop_event, pause_event, server: WebSocketServer):
        from flask import jsonify
        message = {
            "jobID": jobid,
            "taskID": taskid,
            "total":100,
            "progress": 0
        }
        for i in range(100):

            if stop_event.is_set():
                return 1
            while pause_event.is_set():
                time.sleep(0.5)
                if stop_event.is_set():
                    return 1
            time.sleep(10)

            message['progress'] = i+1
            server.add_message_to_queue(json.dumps(message))
        
        return jsonify({"status": Status.SUCCESS.value, "message": "节点运行完成", "result": {}})

    