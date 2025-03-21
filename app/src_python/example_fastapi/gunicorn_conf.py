import asyncio
import logging
import os

from example_fastapi.logic.cache import reset_cache


logger = logging.getLogger('gunicorn.error')

        
async def initialize_app():
    # await db_manager.initialize_database()
    # await initialize_cache()
    # logger.info("Initialization complete.")
    pass

def on_starting(server):
    logger.info("Gunicorn server is starting up.")
    # loop = asyncio.get_event_loop()
    # loop.run_until_complete(initialize_app())

def when_ready(server):
    logger.info("Gunicorn server is ready.")

def on_exit(server):
    logger.info("Gunicorn server is shutting down.")
    loop = asyncio.get_event_loop()
    loop.run_until_complete(reset_cache())
    logger.info("Cache has been reset")

