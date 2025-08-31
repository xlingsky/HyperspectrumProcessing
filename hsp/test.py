import subprocess
import sys
import os
import threading
import time

from hsp import service
from hsp import workflow

def run_script_with_realtime_output(script_path, tag = None):
    if tag is None:
        tag = os.path.splitext(os.path.basename(script_path))[0]
    try:
        process = subprocess.Popen(
            [sys.executable, script_path],
            stdout = subprocess.PIPE,
            stderr = subprocess.PIPE,
            text=True,
            bufsize=1,
            universal_newlines=True
        )

        while True:
            output = process.stdout.readline()
            if output == '' and process.poll() is not None:
                break
            if output:
                print(f'[{tag}:INFO]: {output.strip()}')

        stderr = process.stderr.read()
        if stderr:
            print(f"[{tag}:ERR]: {stderr}")

        return process.returncode
    except Exception as e:
        print(f"[{tag}:ERR]: {e}")
        return -1
            

if __name__ == "__main__":
    dirpath = os.path.dirname(__file__)

    share = dict()
    workflow.tracking_processing(os.path.join(os.path.dirname(
        dirpath), "template/Order.json"), service.Event(), service.Logger(), share)

    exit(0)

    server = threading.Thread(target=run_script_with_realtime_output, args=(os.path.join(dirpath, "service.py"),), daemon=True)
    client = threading.Thread(target=run_script_with_realtime_output, args=(os.path.join(dirpath, "client.py"),), daemon=True)

    server.start()
    time.sleep(5)
    
    client.start()

    client.join()
    server.join()