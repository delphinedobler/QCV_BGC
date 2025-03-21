from .attribute_tmpl import ReactAttrTmpl
import asyncio


class ReactQueueTmpl(ReactAttrTmpl):
    def __init__(self, sem_nb: int = 1):
        super().__init__(default_value=[], sem_nb=sem_nb)
    
    async def get(self):
        """
        Retrieve an item from the queue and synchronize the reactive attribute.
        """
        async with self.semaphore:
            if not self.value:
                raise asyncio.QueueEmpty("Queue is empty")
            value = self.value.pop(0)
        await self.sync_non_reactive_attrs_with_reactive_attrs()
        return value
    
    async def put(self, item):
        """
        Put an item into the queue and synchronize the reactive attribute.
        """
        queue = await self.get_objects()
        queue.append(item)
        await self.set_objects(value=queue)
        # await self.sync_non_reactive_attrs_with_reactive_attrs()
    
    def is_empty(self):
        queue = self.value
        return len(queue) == 0