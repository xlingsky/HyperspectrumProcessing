import time
from watchdog.observers import Observer
from watchdog.events import FileSystemEventHandler
from threading import Timer

class TimeoutFileWatcher(FileSystemEventHandler):
    def __init__(self, newfile_callback, finished_callback, timeout_seconds=30):
        self.callback_newfile = newfile_callback
        self.callback_finished = finished_callback
        self.timeout_seconds = timeout_seconds
        self.last_activity_time = time.time()
        self.observer = None
        self.timeout_timer = None
        
    def start_timeout_timer(self):
        """Start or restart the timeout timer"""
        if self.timeout_timer:
            self.timeout_timer.cancel()
        
        self.timeout_timer = Timer(self.timeout_seconds, self.check_timeout)
        self.timeout_timer.daemon = True
        self.timeout_timer.start()
    
    def check_timeout(self):
        """Check if timeout has been reached"""
        current_time = time.time()
        if current_time - self.last_activity_time >= self.timeout_seconds:
            print(f"No new files for {self.timeout_seconds} seconds. Stopping watch...")
            if self.observer:
                self.observer.stop()
            if self.callback_finished is not None:
                self.callback_finished()
    
    def on_created(self, event):
        """Handle file creation events"""
        import os
        if not event.is_directory:
            file_path = event.src_path
            if self.callback_newfile(os.path.dirname(file_path), os.path.basename(file_path)):
                self.last_activity_time = time.time()
                self.start_timeout_timer()
    
    def  start(self, directory):          
        observer = Observer()
        self.observer = observer

        observer.schedule(self, directory, recursive=True)
        observer.start()
    
        # Start the initial timeout timer
        self.start_timeout_timer()
    
        print(f"Watching {directory} for available images. Will stop after {self.timeout_seconds}s of inactivity.")
    
        # try:
        #     while observer.is_alive():
        #         observer.join(1)
        # except KeyboardInterrupt:
        #     observer.stop()
        # observer.join()
    
    def is_alive(self):
        return self.observer.is_alive()

    def join(self, timeout = None):
        self.observer.join(timeout)