from fastapi import APIRouter, Depends
from .urls import *
from ..routes import (
    users, 
    authentication, 
    tests, 
    cache as cache_routers, 
    objects, 
    sessions,
    workspaces,
    )
# from ..middlewares.authentication import *
from boilerplate_shared_module.middlewares.authentication import *
from boilerplate_shared_module.schemes_pydantic.api_metadata import SoftwareUsageEnum
from ..logic.storage.cache import manager as cache_manager
from ..logic.authentication import keycloak_manager
# from starlette.middleware import Middleware
from starlette.types import ASGIApp
from typing import Type
from fastapi.routing import APIRoute


app_middleware = AppAuthMiddlewares(
    api_meta=api_meta,
    redis_manager=cache_manager.redis_manager,
    keycloak_manager=keycloak_manager,
    redirect=False,
    software_used=SoftwareUsageEnum.fastapi
    )
authentication_middlewares = app_middleware.get_middlewares()


# Authenticated router
router = APIRouter(
    prefix=common_prefix, 
    # dependencies=[Depends(Middleware(AuthMiddleware))]  
    dependencies=authentication_middlewares
    )

router.include_router(sessions.router, tags=['sessions'])
router.include_router(cache_routers.router, tags=['cache'])
router.include_router(workspaces.router, tags=['workspaces'])
# Not authenticated router
not_auth_router = APIRouter(prefix=common_prefix)
not_auth_router.include_router(users.not_auth_router, tags=['users'])
not_auth_router.include_router(tests.not_auth_router, tags=['tests'])
not_auth_router.include_router(objects.not_auth_router, tags=['objects'])
not_auth_router.include_router(authentication.not_auth_router, tags=['authentication'])