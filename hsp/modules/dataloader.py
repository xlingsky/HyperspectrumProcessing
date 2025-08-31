import os
from hsp.utils.filewatcher import FileWatcher

class DataLoader:
    def __init__(self, directory, extensions, callback):
        self.extensions = [ext.lower() for ext in extensions]
        self.callback = callback
        self.directory = directory
        self.known_files = list()
        
        # Initialize with existing files
        self.scan_existing_files(directory)

    def scan_existing_files(self, directory):
        for file in os.listdir(directory):
            if any(file.lower().endswith(ext) for ext in self.extensions):
                self.known_files.append(file)

    # def scan_existing_files_recursive(self, directory):
    #     """Scan for existing files on startup"""
    #     for root, _, files in os.walk(directory):
    #         for file in files:
    #             if any(file.lower().endswith(ext) for ext in self.extensions):
    #                 self.known_files.append(os.path.join(root, file))

    def load_data(self):
        ...
