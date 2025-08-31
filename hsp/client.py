import asyncio
import websockets
import requests
import threading
import os
from hsp.utils import common

def run_pipeline(url):
    algorithm = [
        "ID_WF_SSCL",
        # "ID_WF_MBJC",
        # "ID_WF_MBGZ",
        # "ID_WF_GJSC",
        # "ID_WF_GJRH",
        # "ID_WF_MBSB",
        # "ID_WF_GJYC",
        # "ID_WF_YWPG",
        # "ID_WF_SSCL",
    ]
    for taskid in algorithm:
        jsondata = {
            "JobID": "1",
            "TaskID": taskid,
            "Action": "run",
            "OrderPath": os.path.join(common.parent_dir, "template/Order.json"),
            "Parameters": {
            }
        }
        print(jsondata)
        try:
            # 方法1：使用json参数（推荐）
            response = requests.post(
                url,
                json=jsondata,  # 使用json参数，会自动设置正确的Content-Type
            )
            if response is not None:
                print("状态码:", response.status_code)
                print("响应内容:", response.json())

        except requests.exceptions.Timeout:
            print("请求超时，请检查网络连接或服务器状态")
        except requests.exceptions.ConnectionError:
            print("连接错误，请检查URL和网络连接")
        except requests.exceptions.RequestException as e:
            print(f"请求发生错误: {e}")
        except ValueError:
            print("响应不是有效的JSON格式")
            print("原始响应:", response.text)

async def test_websocket_client(server_host, server_port):
    # 配置服务器地址 - 根据你的实际情况修改
    # server_host = "localhost"  # 如果是本机测试
    uri = f"ws://{server_host}:{server_port}"

    try:
        print(f"尝试连接到: {uri}")
        print("按 Ctrl+C 停止客户端\n")

        async with websockets.connect(uri) as websocket:
            print("连接成功!")
            print("输入消息发送到服务器，输入 'quit' 退出")

            # 测试发送一些初始消息
            test_messages = [
                # "Hello Server!",
                # "这是一条测试消息",
                "{\"type\": \"test\", \"data\": \"JSON数据\"}"
            ]

            while True:

                # 接收服务器响应
                response = await websocket.recv()
                print(f"收到: {response}")

                print("-" * 50)
                await asyncio.sleep(0.1)

    except websockets.exceptions.ConnectionClosed:
        print("连接已关闭")
    except ConnectionRefusedError:
        print(f"连接被拒绝，请检查:")
        print(f"1. 服务器是否运行在 {uri}")
        print(f"2. 防火墙是否开放端口 {server_port}")
        print(f"3. 服务器host设置是否正确")
    except Exception as e:
        print(f"连接错误: {e}")

def run_websocket_client(host, port):
    asyncio.run(test_websocket_client(host, port))

if __name__ == "__main__":
    import sys
    ws_host = '127.0.0.1'
    ws_port = 9876
    flask_host = '127.0.0.1'
    flask_port = 7000
    if len(sys.argv) >= 3:
        ws_host = sys.argv[1]
        ws_port = int(sys.argv[2])
    if len(sys.argv) >= 5:
        flask_host = sys.argv[3]
        flask_port = int(sys.argv[4])

    ws = threading.Thread(target=run_websocket_client, args=(ws_host, ws_port), daemon=True)
    ws.start()

    run_pipeline(f'http://{flask_host}:{flask_port}/hsp/service')
    ws.join()