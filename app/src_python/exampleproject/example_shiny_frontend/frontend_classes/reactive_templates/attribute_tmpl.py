import asyncio
from shiny import reactive
import copy


class ReactAttrTmpl:
    def __init__(self, default_value, sem_nb: int = 1):
        self.semaphore = asyncio.Semaphore(value = sem_nb)
        self.value = default_value
        self.reactive = reactive.value(self.value)
    
    
    # def __repr__(self) -> str:
    #     return self.value
    
    def __call__(self):
        return self.value
    
    async def get_objects(self):
        await self.sync_non_reactive_attrs_with_reactive_attrs()
        async with self.semaphore:
            return copy.deepcopy(self.value)
    
    async def set_objects(self, value):
        async with self.semaphore:
            self.value = value
        await self.sync_non_reactive_attrs_with_reactive_attrs()
        
    async def sync_non_reactive_attrs_with_reactive_attrs(self):
        """
        Synchronize non-reactive attributes with shiny reactive attributes.
        """
        try:
            async with self.semaphore:
                value = self.value
                self.reactive.set(value)
        except Exception as e:
            print(e)
