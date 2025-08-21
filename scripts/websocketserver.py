import asyncio
import websockets
import queue

class WebSocketServer:
    def __init__(self, host="localhost", port=8765):
        self.host = host
        self.port = port
        self.clients = set()
        self.message_queue = queue.Queue()  # Message queue

    async def register_client(self, websocket):
        self.clients.add(websocket)
        await self.notify_clients(f"New client connected: {websocket.remote_address}")

    async def unregister_client(self, websocket):
        self.clients.remove(websocket)
        await self.notify_clients(f"Client disconnected: {websocket.remote_address}")

    async def notify_clients(self, message):
        if self.clients:
            await asyncio.wait([client.send(message) for client in self.clients])

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
            await asyncio.sleep(0.1)  # Sleep to prevent busy waiting

    async def run_server(self, host, port):
        async with websockets.serve(self.handler, host, port):
            await self.send_messages_from_queue()  # Run the message sending loop

    async def handler(self, websocket, path):
        await self.register_client(websocket)
        try:
            while True:
                # Wait for a message from the client, if any
                message = await websocket.recv()
                print(f"Received message from client: {message}")
        finally:
            await self.unregister_client(websocket)

    def start(self):
        start_server = websockets.serve(self.handler, self.host, self.port)
        asyncio.get_event_loop().run_until_complete(start_server)
        asyncio.get_event_loop().run_forever()

# Example usage:
if __name__ == "__main__":
    server = WebSocketServer()
    server.start()