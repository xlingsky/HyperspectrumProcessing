import asyncio
import websockets
import queue

class WebSocketServer:
    def __init__(self):
        self.clients = set()
        self.message_queue = queue.Queue()  # Message queue

    async def register_client(self, websocket):
        self.clients.add(websocket)
        print(f"New client connected: {websocket.remote_address}")

    async def unregister_client(self, websocket):
        self.clients.remove(websocket)
        print(f"Client disconnected: {websocket.remote_address}")

    async def notify_clients(self, message):
        try:
            if self.clients:
                for client in self.clients:
                    await client.send(message)
        except Exception as e:
            self.clients.clear()
            self.clients = set()
            print(f"[Websocket ERROR]: {e}. Clear all clients.")

    async def send_progress(self, progress):
        message = f"Progress: {progress}%"
        await self.notify_clients(message)

    def add_message_to_queue(self, message):
        self.message_queue.put_nowait(message)

    async def send_messages_from_queue(self):
        while True:
            if not self.message_queue.empty():
                message = self.message_queue.get_nowait()
                await self.notify_clients(message)
            await asyncio.sleep(0.001)  # Sleep to prevent busy waiting

    async def sender(self, websocket):
        await self.register_client(websocket)
        try:
            while True:
                if not self.message_queue.empty():
                    message = self.message_queue.get_nowait()
                    await self.notify_clients(message)
                await asyncio.sleep(0.001)  # Sleep to prevent busy waiting
        finally:
            await self.unregister_client(websocket)

    async def run_server(self, host, port):
        server = await websockets.serve(self.sender, host, port)
        await server.wait_closed()

    async def handler(self, websocket):
        await self.register_client(websocket)
        try:
            while True:
                # Wait for a message from the client, if any
                message = await websocket.recv()
                print(f"Received message from client: {message}")
        finally:
            await self.unregister_client(websocket)

# Example usage:
if __name__ == "__main__":
    server = WebSocketServer()
    asyncio.run(server.run_server('0.0.0.0', 8765))