from fastapi.templating import Jinja2Templates
from starlette.requests import Request
from starlette.responses import RedirectResponse
from fastapi import Depends, Header
from datetime import datetime, timezone
import base64
import json
import traceback


from boilerplate_shared_module.schemes_pydantic.orm import *
from boilerplate_shared_module.schemes_pydantic.basic_schemes import *
from boilerplate_shared_module.tools.models import *
from boilerplate_shared_module.consumers.tools import get_app_route
from boilerplate_shared_module.tools.api_infos import ApiMetadata, api_meta_to_rc_dic, ApiTypeEnum

from ..metadata.api import get_api_meta
from ..logic.authentication import *
from ..routers.urls import *


# Define templates dir
templates_dir = '/home/serviceuser/service/app/src_python/example_fastapi/templates'
templates = Jinja2Templates(directory=templates_dir)


async def handler_login(
    request: Request,
    api_meta: ApiMetadata = Depends(get_api_meta),
    ):
    """
    Handles the login route.

    Parameters
    ----------
    request : Request
        Incoming request object.

    Returns
    -------
    ResponseCapsule
        Response capsule containing the data for rendering the login page or redirection.
    """
    api_meta_dic = api_meta_to_rc_dic(api_meta=api_meta)
    # session_id = request.session.get('session_id')
    # if session_id:
    #     try:
    #         session_validity = await is_session_valid(request=request)
    #         if session_validity:
    #             data = RedirectResponse(docs_url)
    #             return ResponseCapsule(data=data, type='RedirectResponse', **api_meta_dic)
    #         else:
    #             # data = RedirectResponse(not_auth_router.url_path_for('route_logout'))
    #             data = RedirectResponse(get_app_route(
    #                 api_type=ApiTypeEnum.main,
    #                 route='logout'
    #                 )
    #             )
    #             return ResponseCapsule(data=data, type='RedirectResponse', **api_meta_dic)
    #     except Exception as e:
    #         # data = RedirectResponse(not_auth_router.url_path_for('route_logout'))
    #         data = RedirectResponse(get_app_route(
    #                 api_type=ApiTypeEnum.main,
    #                 route='logout'
    #                 )
    #             )
    #         return ResponseCapsule(data=data, type='RedirectResponse', **api_meta_dic)
    redirect_url = request.url_for("custom_swagger_ui_html")
    data = templates.TemplateResponse(
        "login.html", 
        {"request": request, "redirect_url": redirect_url}
    )
    return ResponseCapsule(data=data, type='TemplateResponse', **api_meta_dic)


async def handler_logout(
    request: Request,
    redirect_url: str | None = None,
    api_meta: ApiMetadata = Depends(get_api_meta),
    ):
    """
    Handles the logout route.

    Parameters
    ----------
    request : Request
        Incoming request object.
    redirect_url : str | None
        URL where the user will be redirected after being logged out.

    Returns
    -------
    ResponseCapsule
        Response capsule containing the redirection data.
    """
    api_meta_dic = api_meta_to_rc_dic(api_meta=api_meta)
    access_token = request.session.get('access_token')
    if access_token:
        try:
            await logout(request=request, access_token=access_token)
        except Exception as e:
            msg = f'A problem occured while logging out : {e}'
            error_traceback = str(traceback.format_exc())
            return ResponseCapsule(
                message=msg, 
                http_code=500, 
                error=error_traceback, 
                **api_meta_dic,
                )
    else:
        msg = f'User is already logged out'
        return ResponseCapsule(
            message=msg, 
            data=request.session,
            type='dict',
            **api_meta_dic,
            )
    # data = RedirectResponse(not_auth_router.url_path_for('route_login'))
    if redirect_url:
        data = RedirectResponse(redirect_url)
        return ResponseCapsule(
            data=data, 
            type='RedirectResponse', 
            **api_meta_dic)
    else:
        msg = f'User has been logged out'
        return ResponseCapsule(
            message=msg, 
            data=request.session,
            type='dict',
            **api_meta_dic,
            )


async def handler_refresh_token(
    request: Request,
    api_meta: ApiMetadata = Depends(get_api_meta),
    ):
    """
    Handles the refresh token route

    Parameters
    ----------
    request : Request
        Incoming request object.

    Returns
    -------
    ResponseCapsule
        Response capsule containing the new authorization token or error message.
    """
    api_meta_dic = api_meta_to_rc_dic(api_meta=api_meta)
    
    access_token = request.session.get('access_token')
    if not access_token:
        msg = "User isn't authenticated. Please login before refreshing the token"
        error_traceback = msg
        return ResponseCapsule(http_code=401, error=error_traceback, **api_meta_dic)
    try:
        async with db_manager.AsyncSessionLocal() as db:
            authorization = await refresh_token(request=request, db=db)
        # Return result
        typee = 'authorization_token'
        return ResponseCapsule(data=authorization, type=typee)
    
    except Exception as e:
        msg = f'A problem occured while refreshing the token : {e}'
        error_traceback = str(traceback.format_exc())
        return ResponseCapsule(http_code=401, error=error_traceback, **api_meta_dic)
     

async def handler_auth(
    request: Request,
    api_meta: ApiMetadata = Depends(get_api_meta),
    ):
    """
    Handles the authentication route.

    Parameters
    ----------
    request : Request
        Incoming request object.

    Returns
    -------
    ResponseCapsule
        Response capsule containing the redirection data.
    """
    api_meta_dic = api_meta_to_rc_dic(api_meta=api_meta)
    # request_state = {}
    # if redirect_url:
    #     request_state['redirect_url'] = redirect_url
    # Rediriger l'utilisateur vers l'URL d'authentification de Keycloak
    auth_url = keycloak_manager.auth_url(
        authentication_callback_url=route_auth_callback,
        request_state=request.state
    )
    data = RedirectResponse(url=auth_url)
    return ResponseCapsule(data=data, type='RedirectResponse', **api_meta_dic)


async def handler_auth_callback(
    request: Request, 
    code: str,
    state: str,
    redirect_url: str | None = None,
    api_meta: ApiMetadata = Depends(get_api_meta),
    ):
    """
    Handles the authentication callback route.

    Parameters
    ----------
    request : Request
        Incoming request object.
    code : str
        Authorization code.
    state: str
        State of the incoming request
    redirect_url : str | None
        URL where the user will be redirected after authentication.

    Returns
    -------
    ResponseCapsule
        Response capsule containing the redirection data.
    """
    api_meta_dic = api_meta_to_rc_dic(api_meta=api_meta)
    async with db_manager.AsyncSessionLocal() as db:
        auth = await add_token(request=request, code=code, db=db)
    
    if not redirect_url:
        redirect_url = request.session.get('redirect_url')
    if redirect_url:
        # Remove redirect url from session
        request.session.pop('redirect_url', {})
        data = RedirectResponse(url=redirect_url)
        return ResponseCapsule(data=data, type='RedirectResponse', **api_meta_dic)
    return ResponseCapsule(data=auth, type='AuthorizationToken', **api_meta_dic)


# async def handler_get_current_user(
#     request: Request,
#     api_meta: ApiMetadata = Depends(get_api_meta),
#     ):
#     """
#     Handles the route to get information about the current user.

#     Parameters
#     ----------
#     request : Request
#         Incoming request object.

#     Returns
#     -------
#     ResponseCapsule
#         Response capsule containing the user information.
#     """
#     api_meta_dic = api_meta_to_rc_dic(api_meta=api_meta)
#     # authorization = request.session.get('authorization')
#     # session_id = request.session.get('session_id')
#     access_token = request.session.get('access_token')
#     data = {'message': "User isn't connected"}
#     # if session_id:
#     if access_token:
#         session_info = await redis_manager.get_session_cache(access_token=access_token)
#         user_email = session_info.user_email
#         if user_email:
#             user = await redis_manager.get_user_cache(user_email=user_email)
#         # user_infos : ApiUserInfos = await get_user_infos(authorization=authorization,
#         #                                                 session_id=session_id)
#         # user = user_infos.user
#             data = {"message": f"Hello {user.first_name} {user.last_name}"}
            
#     return ResponseCapsule(data=data, type='dict', **api_meta_dic)