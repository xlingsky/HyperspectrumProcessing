from flask import Flask, request, jsonify
from enum import Enum
from hsp.modules.event import Event
import workflow 
from hsp.modules.websocketserver import WebSocketServer
import asyncio
import threading
import os
import subprocess
import ast
import input_check
import json


app = Flask("HSP")
server = WebSocketServer()

class Logger:
    def __init__(self):
        self.message = {
            "step": 0,
            "total":0,
            "progress": 0
        }
    def progress_update(self, progress):
        self.message["progress"] = progress
        print(json.dumps(self.message))

    def set_total(self, total):
        self.message["total"] = total

    def set_step(self, step):
        self.message["step"] = step

    def set_taskid(self, taskid):
        self.message["taskID"] = taskid

class ServerLogger(Logger):
    def __init__(self, jobid, taskid, server):
        super().__init__()
        self.message['jobID'] = jobid
        self.message['taskID'] = taskid
        self.server = server
    def progress_update(self, progress):
        self.message["progress"] = progress
        server.add_message_to_queue(json.dumps(self.message))

class Status(Enum):
    SUCCESS=200
    ERROR=500

class Action(Enum):
    START = "run"
    PAUSE = "mount"
    KILL = "terminate"
    FETCH_IMAGES = "showImage"
    FETCH_VIDEO = "showVideo"

algorithm = {
    "ID_WF_MBJC": workflow.detection_processing,
    "ID_WF_MBGZ": workflow.tracking_processing,
    "ID_WF_GJSC": workflow.automatic_processing,
    "ID_WF_GJRH": workflow.automatic_processing,
    "ID_WF_MBSB": workflow.automatic_processing,
    "ID_WF_GJYC": workflow.automatic_processing,
    "ID_WF_YWPG": workflow.automatic_processing,
    "ID_WF_SSCL": workflow.automatic_processing
}

tasklist = dict()

@app.route('/hsp/service', methods=['POST'])
def message_processing():
    global tasklist
    data = request.json
    hdr = [data.get('JobID'), data.get('TaskID'), data.get('Action')]

    if not hdr[0] or not hdr[1] or not hdr[2]:
        return jsonify({"status":Status.ERROR.value, "message": "缺少信息头信息", "result":{}})

    action = Action(hdr[2])
    job = tasklist.get(hdr[0])
    taskid = hdr[1]
    if job is None:
        if action != Action.START:
            return jsonify({"status":Status.ERROR.value, "message": "节点未运行", "result":{}})
        tasklist[hdr[0]] = {'logger': ServerLogger(data.get('JobID'), taskid, server), 'share': dict()}
        job = tasklist[hdr[0]]
        task = None
    else:
        task = next(reversed(job.items()))
        if task[0] != taskid and not task[1].is_finished():
            return jsonify({"status":Status.ERROR.value, "message": "上个节点未完成", "result":{}})

    if action == Action.START:
        if task is not None and task[0] == taskid:
            if not task[1].resume():
                return jsonify({"status": Status.ERROR.value, "message": "节点重启失败", "result": {}})
            else:
                return jsonify({"status": Status.SUCCESS.value, "message": "节点重启成功", "result": {}})
        job[taskid] = Event()
        flag, ret = algorithm[taskid](data.get('OrderPath'), job[taskid], job['logger'], job['share'])
        if flag :
            return jsonify({"status": Status.SUCCESS.value, "message": "节点运行完成", "result": ret})
        else:
            return jsonify({"status": Status.ERROR.value, "message": ret, "result": {}})
    elif action == Action.PAUSE:
        assert(task is not None)
        if not task[1].pause():
            return jsonify({"status":Status.ERROR.value, "message": "节点暂停失败", "result":{}})
        return jsonify({"status":Status.SUCCESS.value, "message": "节点暂停成功", "result":{}})
    elif action == Action.KILL:
        assert(task is not None)
        if not task[1].terminate():
            return jsonify({"status":Status.ERROR.value, "message": "节点中止失败", "result":{}})
        return jsonify({"status":Status.SUCCESS.value, "message": "节点中止成功", "result":{}})
    else:
        return jsonify({"status":Status.ERROR.value, "message": "无效操作", "result":{}})

@app.route('/Adsb', methods=['POST'])
def analyze_adsb():
    dat_paths = request.data.decode('utf-8').strip()
    dat_paths = ast.literal_eval(dat_paths)
    results = []
    print(dat_paths)

    if not dat_paths:
        return jsonify({"status": 500, 'message': 'No path provided'})


    try:
        for dat_path in dat_paths:
            if not os.path.exists(dat_path):
                return jsonify({"status": 500, 'message': 'file does not exist'})
            result = subprocess.run(
                ['/home/nuaa/adsb_json.sh', dat_path],
                stdout=subprocess.PIPE,
                stderr=subprocess.PIPE,
                text=True
            )
            print(result.stdout)
            results.append(result.stdout.rstrip())

        result = {"result": results}
        return jsonify({"status": 200, 'message': "Success",'result':result})

    except Exception as e:
        return jsonify({"status": 500, 'message': str(e)})

@app.route('/inputCheck', methods=['POST'])
def analyze_images_endpoint():
    folder_path = request.data.decode('utf-8').strip()
    print(folder_path)

    if not folder_path:
        return jsonify({"status": "500", 'error': 'No folder path provided'})

    if not os.path.exists(folder_path):
        return jsonify({"status": "500", 'error': 'Folder does not exist'})

    try:
        result = input_check.analyze_images(folder_path)
        if result is None:
            return jsonify({"status": "500", 'error': 'No valid images found in folder'})

        result_json = {"status": 200, "message": "Success", "result": {}}

        result.update({
            'folder_path': folder_path,
        })
        result_json["result"] = result
        return jsonify(result_json)
    except Exception as e:
        return jsonify({"status": "500", 'error': str(e)})

def run_websocket_server(host, port):
    asyncio.run(server.run_server(host, port))
def run_flask(host, port):
    app.run(host=host, port=port, debug=False)

if __name__ == '__main__':
    import sys
    ws_host = '0.0.0.0'
    ws_port = 9876
    flask_host = '127.0.0.1'
    flask_port = 7000
    if len(sys.argv) >= 3:
        ws_host = sys.argv[1]
        ws_port = int(sys.argv[2])
    if len(sys.argv) >= 5:
        flask_host = sys.argv[3]
        flask_port = int(sys.argv[4])

    ws_thread = threading.Thread(target=run_websocket_server, args=(ws_host, ws_port), daemon=True)
    ws_thread.start()

    run_flask(flask_host, flask_port)
    ws_thread.join()
