from ..app_metadata import *
from boilerplate_shared_module.consumers.tools import get_app_route
import os
from shiny import Session


keycloak_config_path = '/home/serviceuser/service/dev-secrets/ServiceClient.json'

# Main authentication routes
ROUTE_FRONT_LOGIN_PAGE = '/login'
ROUTE_AUTH = '/auth'
ROUTE_AUTH_CALLBACK = '/auth_callback'
ROUTE_LOGOUT = '/logout'
        

def get_logout_url():
    logout_url = get_app_route(
        app_meta=app_meta,
        route=ROUTE_LOGOUT
        )
    return logout_url


async def logout_session(session: Session):
    # session.http_conn.session.clear()
    # session.http_conn.cookies.clear()
    await session.send_custom_message('redirect', {
                                                    'url': get_logout_url(),
                                                    }
    )
