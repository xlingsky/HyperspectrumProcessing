from flask import Flask, request, jsonify
from enum import Enum
from task import Task, Status
from workflow import object_detection, object_tracking, trajectory_generation, trajectory_fusion, object_recognition, trajectory_prediction, evaluation, realtime_processing
from websocketserver import WebSocketServer
import asyncio

app = Flask("HSP")
server = WebSocketServer()

class Action(Enum):
    START = "run"
    PAUSE = "mount"
    KILL = "terminate"
    FETCH_IMAGES = "showImage"
    FETCH_VIDEO = "showVideo"

algorithm = {
    "ID_WF_MBJC": object_detection,
    "ID_WF_MBGZ": object_tracking,
    "ID_WF_GJSC": trajectory_generation,
    "ID_WF_GJRH": trajectory_fusion,
    "ID_WF_MBSB": object_recognition,
    "ID_WF_GJYC": trajectory_prediction,
    "ID_WF_YWPG": evaluation,
    "ID_WF_SSCL": realtime_processing
}

tasklist = dict()

async def main():
    await server.run_server('192.168.2.200', 90)
    app.run(host='0.0.0.0', port=5000, debug=True)


@app.route('/hsp/service', methods=['POST'])
def message_processing():
    global tasklist
    data = request.json
    hdr = {data.get('JobID'), data.get('TaskID'), data.get('Action')}

    if not hdr[0] or not hdr[1] or not hdr[2]:
        return jsonify({"status":Status.ERROR.value, "message": "缺少信息头信息", "result":{}})

    action = Action(hdr[2])
    job = tasklist.get(hdr[0])
    if job is None:
        tasklist[hdr[0]] = dict()
        job = tasklist[hdr[0]]
        task = None
        if action != Action.START:
            return jsonify({"status":Status.ERROR.value, "message": "节点未运行", "result":{}})
    else:
        task = next(reversed(job.items()))
        if task[0] != taskid and not task[1].finished():
            return jsonify({"status":Status.ERROR.value, "message": "上个节点未完成", "result":{}})

    taskid = hdr[1]

    if action == Action.START:
        if task is not None and task[0] == taskid:
            if not task[1].resume():
                return jsonify({"status": Status.ERROR.value, "message": "节点重启失败", "result": {}})
        task = Task()
        if task.start(algorithm[taskid], hdr[0], hdr[1], server):
            job[taskid] = task
            return True
        else:
            return jsonify({"status":Status.ERROR.value, "message": "节点启动失败", "result":{}})
    elif action == Action.PAUSE:
        assert(task is not None)
        if not task[1].pause():
            return jsonify({"status":Status.ERROR.value, "message": "节点暂停失败", "result":{}})
    elif action == Action.KILL:
        assert(task is not None)
        if not task[1].terminate():
            return jsonify({"status":Status.ERROR.value, "message": "节点暂停失败", "result":{}})
    else:
        return jsonify({"status":Status.ERROR.value, "message": "无效操作", "result":{}})
        

if __name__ == '__main__':
    asyncio.run(main())
