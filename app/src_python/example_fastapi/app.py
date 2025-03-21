import sys
import logging
import asyncio

from pydantic import ValidationError
from fastapi import FastAPI, HTTPException, Request, Depends
from fastapi import status
from fastapi.responses import JSONResponse, RedirectResponse
from fastapi.exceptions import RequestValidationError
from fastapi.openapi.utils import get_openapi
from fastapi.openapi.docs import (
    get_redoc_html,
    get_swagger_ui_html,
    get_swagger_ui_oauth2_redirect_html,
)
from starlette.requests import Request as StarletteRequest
from starlette.middleware.sessions import SessionMiddleware
from humps import camelize
from contextlib import asynccontextmanager



from boilerplate_shared_module.schemes_pydantic.basic_schemes import ResponseCapsule
from boilerplate_shared_module.consumers.tools import get_app_route
from boilerplate_shared_module.exceptions.redirect import RedirectNeededException


# from .logic.storage.cache.manager import redis_manager
from .metadata.api import get_api_meta, app_type
from .middlewares.sessions import *
from .routers.routers import *
from .routers.urls import docs_url, open_api_url
from .handlers.authentication import *


logger = logging.getLogger('uvicorn.error')
logger.setLevel(logging.DEBUG)
# StreamHandler for the console
stream_handler = logging.StreamHandler(sys.stdout)
log_formatter = logging.Formatter("%(asctime)s [%(processName)s: %(process)d] [%(threadName)s: %(thread)d] [%(levelname)s] %(name)s: %(message)s")
stream_handler.setFormatter(log_formatter)
logger.addHandler(stream_handler)
logger.info('API is starting up')




# APP INIT VARS

redis_app_lock_var = 'APP_INIT_LOCK'
LOCK_TIMEOUT = 60  # seconds

# APP INIT FUNCTIONS ------------------------

async def sync_with_db_routine(time_s=40):
    lock = None
    redis_lock_key = "cache_sync_lock"
    try:
        lock = await redis_manager.acquire_lock(redis_lock_key)
        if not lock:
            logger.info("Another instance is already running the sync task.")
        else:
            logger.info("Launching cache sync routine ...")
            while True:
                await sync_cache_with_db_periodically(time_s=time_s)
    except asyncio.CancelledError:
        logger.info("Sync task cancelled. Releasing lock...")
    except Exception as e:
        logger.error(f"Error in sync routine: {e}")
        logger.error(traceback.format_exc())
    finally:
        if lock:
            await redis_manager.release_lock(lock)



async def initialize_app(from_db: bool = False):
    try:
        # if not os.path.exists(app_init_lock_path):
        async with redis_manager.async_redis_client as redis:
            lock_acquired = await redis.set(redis_app_lock_var, "locked", ex=LOCK_TIMEOUT, nx=True)
            # with open(redis_manager.init_lock_path, 'w') as lock_file:
            if lock_acquired:
                await db_manager.initialize_database()
                # lock_file.write("Db initialized\n")
                logger.info("Database initialization complete.")
                await initialize_cache(from_db=from_db)
                # lock_file.write("Cache initialized\n")
                logger.info("Cache initialization complete.")
            else:
                logger.info("Another worker is handling initialization.")
    except Exception as e:
        logger.error(str(e))
        logger.error(traceback.format_exc())
        
    
# ---------------------------------------------------------------------------------------------------


# FASTAPI APP


@asynccontextmanager
async def lifespan(app: FastAPI):
    # await initialize_database()
    # await initialize_cache()
    # cache_task = asyncio.create_task(sync_cache_with_db_periodically(time_s=15))
    await initialize_app(
        from_db=True,
        )
    cache_task = asyncio.create_task(sync_with_db_routine(time_s=40))
    yield
    cache_task.cancel()

api_meta: ApiMetadata = get_api_meta()
api_meta_dic = api_meta_to_rc_dic(api_meta=api_meta)
 
app = FastAPI(
    title=api_meta.name,
    description="Python Boilerplate API",
    version=api_meta.version,
    docs_url='',
    openapi_url='',
    # docs_url=docs_url,
    # openapi_url=open_api_url,
    # openapi_prefix='test'
    lifespan=lifespan
)

# Authentification middleWare
# app.add_middleware(AuthMiddleware,
#                 #    allow_credentials=True,
#                 #    allow_methods=["*"],
#                 #    allow_headers=["*"]
#                    )

# Necessary
# Add middleware that converts cookies into session
app.add_middleware(
    AddCookiesToSessionMiddleware
)
# Add session middleware at the first layer
app.add_middleware(
    SessionMiddleware, secret_key="Z1cgzqCQVJCebnEm4K130Ga80iuu1LSLRbEH20Od",
    )

# app.mount("/static", StaticFiles(directory="/home/serviceuser/service/app/src_python/fast_api/static"), name="static")
# app.add_middleware(AuthMiddleware)


@router.get("/docs", include_in_schema=False)
async def custom_swagger_ui_html():
    return get_swagger_ui_html(
        openapi_url=open_api_url,
        title="Custom API Docs",
        swagger_favicon_url="https://fastapi.tiangolo.com/img/favicon.png"
    )

@router.get("/openapi.json", include_in_schema=False)
async def custom_openapi():
    return JSONResponse(get_openapi(
        title=app.title,
        version=app.version,
        description=app.description,
        routes=app.routes,
    ))


app.include_router(router=router)
app.include_router(router=not_auth_router)

logger.info(f"Login page is located at {get_app_route(route='login', api_type=app_type)}")

@app.exception_handler(RedirectNeededException)
def redirect_needed_exception_handler(request: Request, exc: RedirectNeededException):
    logger.info(exc.message)
    return RedirectResponse(url=exc.redirect_url)

@app.exception_handler(HTTPException)
def my_http_exception_handler(request: Request, exc: HTTPException):
    camel_response = camelize(exc.detail) # convert keys to camelCase
    return JSONResponse(status_code=exc.status_code, content=camel_response)

# Rewrite exception for unprocessable entity
@app.exception_handler(RequestValidationError)
def validation_exception_handler(request: Request, exc: RequestValidationError):
    status_code = status.HTTP_422_UNPROCESSABLE_ENTITY
    error = '. '.join([err.get('msg') for err in exc.errors()])
    msg = 'Unprocessable entity'
    detail = ResponseCapsule(message=msg, http_code=status_code, error=error, **api_meta_dic).model_dump()
    return JSONResponse(status_code=status_code, content=detail)

