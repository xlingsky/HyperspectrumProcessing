from enum import Enum

class Event:
    class Status(Enum):
        WAIT = 0
        RUNNING = 1
        PAUSED = 2
        FINISHED = 3
        KILLED = 4
        
    def __init__(self):
       self._status = Event.Status.RUNNING
    def resume(self):
        self._status = Event.Status.RUNNING
        return True
    def start(self):
        self._status = Event.Status.RUNNING
        return True #self._worker(target, argv, jobid, taskid, server)
    def pause(self):
        self._status = Event.Status.PAUSED
        return True
    def terminate(self):
        self._status = Event.Status.KILLED
        return True

    def is_paused(self):
        return self._status == Event.Status.PAUSED
    def is_terminated(self):
        return self._status == Event.Status.KILLED
    def is_running(self):
        return self._status == Event.Status.RUNNING
    def is_finished(self):
        return self._status == Event.Status.FINISHED
            

    