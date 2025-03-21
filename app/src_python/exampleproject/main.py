from starlette.applications import Starlette
from starlette.routing import Mount, Route
from starlette.requests import Request
from starlette.responses import RedirectResponse, HTMLResponse
from starlette.middleware import Middleware
from starlette.middleware.base import BaseHTTPMiddleware
from starlette.middleware.sessions import SessionMiddleware
from example_shiny_frontend.not_authenticated_page import NotAuthenticatedPage
from example_shiny_frontend.authenticated_page import AuthenticatedPage
from example_shiny_frontend.frontend_classes.objects_store import ObjectStore
from boilerplate_shared_module.middlewares.authentication import *
from example_shiny_frontend.login import LoginPage
from example_shiny_frontend.tools.proxy import *
from example_shiny_frontend.tools.authentication import *
from boilerplate_shared_module.consumers.tools import get_app_route
from example_shiny_frontend import app_metadata
from example_shiny_frontend.admin import AdminPortal
import logging
import sys


logger = logging.getLogger('shiny')
logger.setLevel(logging.DEBUG)
# StreamHandler for the console
stream_handler = logging.StreamHandler(sys.stdout)
log_formatter = logging.Formatter("%(asctime)s [%(processName)s: %(process)d] [%(threadName)s: %(thread)d] [%(levelname)s] %(name)s: %(message)s")
stream_handler.setFormatter(log_formatter)
logger.addHandler(stream_handler)
logger.info('Shiny is starting up')
logger.propagate = True


# ███╗   ███╗ █████╗ ██╗███╗   ██╗     █████╗ ██████╗ ██████╗     ███╗   ███╗███████╗████████╗ █████╗ 
# ████╗ ████║██╔══██╗██║████╗  ██║    ██╔══██╗██╔══██╗██╔══██╗    ████╗ ████║██╔════╝╚══██╔══╝██╔══██╗
# ██╔████╔██║███████║██║██╔██╗ ██║    ███████║██████╔╝██████╔╝    ██╔████╔██║█████╗     ██║   ███████║
# ██║╚██╔╝██║██╔══██║██║██║╚██╗██║    ██╔══██║██╔═══╝ ██╔═══╝     ██║╚██╔╝██║██╔══╝     ██║   ██╔══██║
# ██║ ╚═╝ ██║██║  ██║██║██║ ╚████║    ██║  ██║██║     ██║         ██║ ╚═╝ ██║███████╗   ██║   ██║  ██║
# ╚═╝     ╚═╝╚═╝  ╚═╝╚═╝╚═╝  ╚═══╝    ╚═╝  ╚═╝╚═╝     ╚═╝         ╚═╝     ╚═╝╚══════╝   ╚═╝   ╚═╝  ╚═╝

# Custom these in `app_metadata.py`                                                                                                    
software_used = app_metadata.software_used
app_meta = app_metadata.app_meta


# ███████╗███╗   ██╗██████╗ ██████╗  ██████╗ ██╗███╗   ██╗████████╗███████╗
# ██╔════╝████╗  ██║██╔══██╗██╔══██╗██╔═══██╗██║████╗  ██║╚══██╔══╝██╔════╝
# █████╗  ██╔██╗ ██║██║  ██║██████╔╝██║   ██║██║██╔██╗ ██║   ██║   ███████╗
# ██╔══╝  ██║╚██╗██║██║  ██║██╔═══╝ ██║   ██║██║██║╚██╗██║   ██║   ╚════██║
# ███████╗██║ ╚████║██████╔╝██║     ╚██████╔╝██║██║ ╚████║   ██║   ███████║
# ╚══════╝╚═╝  ╚═══╝╚═════╝ ╚═╝      ╚═════╝ ╚═╝╚═╝  ╚═══╝   ╚═╝   ╚══════╝
                                                                         
ROUTE_NON_ROUTE_AUTH_HOMEPAGE = '/anonymous'
ROUTE_AUTH_HOMEPAGE = '/home'
ROUTE_DEFAULT_PAGE = '/'


#  ██████╗ ██████╗ ███╗   ███╗██████╗  ██████╗ ███╗   ██╗███████╗███╗   ██╗████████╗███████╗
# ██╔════╝██╔═══██╗████╗ ████║██╔══██╗██╔═══██╗████╗  ██║██╔════╝████╗  ██║╚══██╔══╝██╔════╝
# ██║     ██║   ██║██╔████╔██║██████╔╝██║   ██║██╔██╗ ██║█████╗  ██╔██╗ ██║   ██║   ███████╗
# ██║     ██║   ██║██║╚██╔╝██║██╔═══╝ ██║   ██║██║╚██╗██║██╔══╝  ██║╚██╗██║   ██║   ╚════██║
# ╚██████╗╚██████╔╝██║ ╚═╝ ██║██║     ╚██████╔╝██║ ╚████║███████╗██║ ╚████║   ██║   ███████║
#  ╚═════╝ ╚═════╝ ╚═╝     ╚═╝╚═╝      ╚═════╝ ╚═╝  ╚═══╝╚══════╝╚═╝  ╚═══╝   ╚═╝   ╚══════╝
     
# Cache manager + authentication backend config
if software_used == SoftwareUsageEnum.simple_shiny:
    #cache manager
    redis_manager = None
    # Authentication backend
    authentication_backend_type = None
else:
    #cache manager
    redis_manager = RedisManager()
    # Authentication backend: Can be changed if authentication is made on other backend
    authentication_backend_type=ApiTypeEnum.main

# Object store                                                                                     
object_store = ObjectStore(
    front_ROUTE_AUTH_HOMEPAGE=get_app_route(
        app_meta=app_meta,
        route=ROUTE_AUTH_HOMEPAGE,
        ),
    front_app_meta=app_meta,
    redis_manager=redis_manager,
    software_used=software_used,
    backend_type=authentication_backend_type,
    )

# Middleware
app_auth_middleware_class = AppAuthMiddlewares(
    api_meta=app_meta,
    redis_manager=object_store.redis_manager,
    keycloak_manager=object_store.keycloak_manager,
    redirect=True, # True because it is a frontend
    software_used=object_store.software_used,
    )

auth_middlewares = app_auth_middleware_class.get_middlewares()


# ██████╗  ██████╗ ██╗   ██╗████████╗███████╗███████╗    ██████╗ ███████╗███████╗
# ██╔══██╗██╔═══██╗██║   ██║╚══██╔══╝██╔════╝██╔════╝    ██╔══██╗██╔════╝██╔════╝
# ██████╔╝██║   ██║██║   ██║   ██║   █████╗  ███████╗    ██║  ██║█████╗  █████╗  
# ██╔══██╗██║   ██║██║   ██║   ██║   ██╔══╝  ╚════██║    ██║  ██║██╔══╝  ██╔══╝  
# ██║  ██║╚██████╔╝╚██████╔╝   ██║   ███████╗███████║    ██████╔╝███████╗██║     
# ╚═╝  ╚═╝ ╚═════╝  ╚═════╝    ╚═╝   ╚══════╝╚══════╝    ╚═════╝ ╚══════╝╚═╝     
                                                                               
async def default_page(request):
    return HTMLResponse('Hello, you')


async def authenticated_ping(request):
    return HTMLResponse('PONG')


async def auth(request):
    redirect_url = request.query_params.get('redirect_url')
    if redirect_url:
        request.session['redirect_url'] = redirect_url
    # Get the keycloak url for simple shiny app; or the backend auth url for advanced shiny app
    keycloak_auth_url = object_store.get_keycloak_auth_url(
        request_state=request.state,
    )
    if redirect_url:
        keycloak_auth_url += f"?redirect_url={redirect_url}"
    return RedirectResponse(keycloak_auth_url)
    
    
async def auth_callback(request):
    # This callback is only called for simple shiny app, elsewhere it is backend callback which is called
    auth_code = request.query_params.get('code')
    authorization = object_store.keycloak_manager.create_authorization_token(
        code=auth_code,
        redirect_url=get_app_route(
            app_meta=object_store.app_meta, 
            route=ROUTE_AUTH_CALLBACK,
            )
    )
    print(f"{authorization=}")
    request.session['authorization'] = authorization
    redirect_url = request.session.get('redirect_url')
    if redirect_url:
        # Remove redirect url from session
        request.session.pop('redirect_url', {})
    # If no redirect url is specified, redirect to default homepage
    else:
        redirect_url = object_store.front_ROUTE_AUTH_HOMEPAGE
    return RedirectResponse(redirect_url)


async def logout(request: Request):
    if object_store.software_used == SoftwareUsageEnum.simple_shiny:
        auth = request.session.get("authorization", None)
        if auth is not None:
            try:
                refresh_token = auth.get('refresh_token')
                object_store.keycloak_manager.open_id.logout(refresh_token)
            except Exception as e:
                print(e)
            # Clear http session and cookies
            request.session.pop("authorization", None)
            request.cookies.pop("authorization", None)
            
            # Redirect to login
            return RedirectResponse(ROUTE_FRONT_LOGIN_PAGE)
    elif object_store.software_used == SoftwareUsageEnum.advanced_shiny:
        # Use backend logout to logout session (and so remove token from the session)
        backend_logout_route = get_app_route(
            api_type=object_store.backend_type,
            route='logout',
        )
        backend_logout_url = f"{backend_logout_route}?redirect_url={ROUTE_FRONT_LOGIN_PAGE}"
        return RedirectResponse(backend_logout_url)
    else:
        raise NotImplementedError(
            f"`logout` isn't implemented yet for {object_store.software_used=}"
        )


auth_middlewares = auth_middlewares

not_authenticated_routes = [
    Route(ROUTE_DEFAULT_PAGE,endpoint=default_page),
    Route(ROUTE_AUTH, endpoint=auth),
    Route(ROUTE_AUTH_CALLBACK, endpoint=auth_callback),
    Route(ROUTE_LOGOUT, endpoint=logout),
    # Route('/login_logical', endpoint=login_logical),
    Mount(
        ROUTE_FRONT_LOGIN_PAGE, 
        app=LoginPage(object_store=object_store).app, 
        name='login'
        ),
    Mount(ROUTE_NON_ROUTE_AUTH_HOMEPAGE, app=NotAuthenticatedPage(object_store=object_store).app),
]

authenticated_routes = [
    Mount(ROUTE_AUTH_HOMEPAGE, app=AuthenticatedPage(object_store=object_store).app, 
        ),
    Route(
        '/authenticated_ping',
        endpoint=authenticated_ping,
        ),
    Mount('/admin', app=AdminPortal(object_store=object_store).app, 
        ),
]

def wrap_with_middleware(route, middlewares):
    if isinstance(route, Mount):
        return Mount(
            route.path,
            app=route.app,
            middleware=middlewares
        )
    elif isinstance(route, Route):
        return Route(
            route.path,
            endpoint=route.endpoint,
            methods=route.methods,
            name=route.name,
            middleware=middlewares
        )
    return route

# Apply authentication middleware to authenticated routes / mounts
authenticated_routes_with_middleware = [
    wrap_with_middleware(route, auth_middlewares) for route in authenticated_routes
]

app_routes = authenticated_routes_with_middleware + \
             not_authenticated_routes


#  █████╗ ██████╗ ██████╗ 
# ██╔══██╗██╔══██╗██╔══██╗
# ███████║██████╔╝██████╔╝
# ██╔══██║██╔═══╝ ██╔═══╝ 
# ██║  ██║██║     ██║     
# ╚═╝  ╚═╝╚═╝     ╚═╝     
                        
app_middlewares = [
    Middleware(SessionMiddleware, secret_key="Z1cgzqCQVJCebnEm4K130Ga80iuu1LSLRbEH20Od"),
    # Middleware(AuthMiddleware)
]

app = Starlette(
    routes=app_routes, 
    middleware=app_middlewares,
    )
