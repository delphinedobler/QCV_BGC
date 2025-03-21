from starlette.requests import Request
from starlette.responses import RedirectResponse
from starlette.middleware.base import BaseHTTPMiddleware
import logging
import json

from boilerplate_shared_module.exceptions.redirect import RedirectNeededException

from ..handlers.authentication import *

        
logging_level = logging.DEBUG
logger = logging.getLogger('AuthMiddleware')
logging.basicConfig(level=logging_level)

#  Configure SysLogHandler
# syslog_handler = logging.handlers.SysLogHandler(address='/dev/log')
# syslog_handler = logging.StreamHandler()
# formatter = logging.Formatter(
#     '%(asctime)s [%(processName)s: %(process)d] [%(threadName)s: %(thread)d] [%(levelname)s] %(name)s: %(message)s',
#     datefmt='%Y-%m-%d %H:%M:%S'
# )
# syslog_handler.setFormatter(formatter)
# logger.addHandler(syslog_handler)

# bind the middleware to an APIRoute subclass


class AppMiddlewares:
    def __init__(
        self, 
        api_meta: ApiMetadata,
        redirect: bool = False,
        is_fastapi: bool = True,
        ):
        self.api_meta: ApiMetadata = api_meta
        self.redirect: bool = redirect
        self.is_fastapi: bool = is_fastapi
        
    def get_middleware_redirection(
        self,
        redirect_url: str, 
        msg: str | None = None,
        ):
        if self.is_fastapi:
            if msg is None:
                msg = f"Middleware: {status.HTTP_401_UNAUTHORIZED} - redirecting to {redirect_url}..."
            raise RedirectNeededException(
                redirect_url=redirect_url, 
                message=msg
                )
        else:
            return RedirectResponse(
                url=redirect_url, 
                status_code=status.HTTP_401_UNAUTHORIZED
                ) 
            
            
    def get_auth_info_in_request(self, request: Request):
        login_url = get_app_route(
                    app_meta=api_meta,
                    route='login'
                    )
        logger.debug("URL called: %s", request.url)
        # Try getting the session_id
        session_id = request.session.get('session_id')
        if not session_id:
            logger.debug("Session id not found in the session")
            encoded_session_id = request.cookies.get('session_id')
            if encoded_session_id:
                logger.debug("Session id found in the cookies")
                # Décoder en base64
                session_id_json_bytes = base64.b64decode(encoded_session_id)

                # Convertir les octets en chaîne JSON
                session_id_json = session_id_json_bytes.decode('utf-8')

                # Charger la chaîne JSON en objet Python (dictionnaire, liste, etc.)
                session_id = json.loads(session_id_json)
                request.session['session_id'] = session_id
            else:
                logger.debug("Session id not found (cookies + session)")
                status_code = status.HTTP_401_UNAUTHORIZED
                msg = 'Session id not found in the cookies and / or the session. Please login.'
                if self.redirect:
                    response = self.get_middleware_redirection(
                        redirect_url=login_url, 
                        msg=msg,
                        )
                    return response
                else:
                    error = msg
                    response = ResponseCapsule(
                        http_code=status_code,
                        error=error,
                        message=msg,
                        **api_meta_to_rc_dic(api_meta=api_meta)
                    )
                    raise HTTPException(status_code=response.http_code, detail=response.model_dump())

        # Try getting the user_email
        user_email = request.session.get('user_email')
        if not user_email:
            logger.debug("User email not found in the session")
            encoded_user_email = request.cookies.get('user_email')
            if encoded_user_email:
                logger.debug("User email found in the cookies")
                # Décoder en base64
                user_email_json_bytes = base64.b64decode(encoded_user_email)

                # Convertir les octets en chaîne JSON
                user_email_json = user_email_json_bytes.decode('utf-8')

                # Charger la chaîne JSON en objet Python (dictionnaire, liste, etc.)
                user_email = json.loads(user_email_json)
                request.session['user_email'] = user_email
            else:
                logger.debug("User email not found (cookies + session)")
                status_code = status.HTTP_401_UNAUTHORIZED
                msg = 'User email not found in the cookies and / or the session. Please login.'
                if self.redirect:
                    response = self.get_middleware_redirection(
                        redirect_url=login_url, 
                        msg=msg,
                        )
                    return response
                else:
                    error = msg
                    response = ResponseCapsule(
                        http_code=status_code,
                        error=error,
                        message=msg,
                        **api_meta_to_rc_dic(api_meta=api_meta)
                    )
                    raise HTTPException(status_code=response.http_code, detail=response.model_dump())
        return


    async def ensure_token_validity(
        self,
        request: Request, 
        ):
        session_id = request.session.get('session_id')
        logger.debug(f"session id is {session_id}")
        session_infos: ApiSessionInfo = await redis_manager.get_session_cache(session_id=session_id)    
        # sessions_cache = await get_sessions_cache()
        # print(f"{sessions_cache=}")
        if not session_infos:
            logger.debug("Any session was found in the cache")
            redirect_url = get_app_route(
                    api_type=ApiTypeEnum.main,
                    route='login'
                    )
            return RedirectResponse(redirect_url)

        # Get the token
        token = session_infos.session.token
        # Redirect to login page if token not found
        if not token:
            logger.debug("Access token not found in the cache session")
            status_code = status.HTTP_401_UNAUTHORIZED
            msg = 'Access token not found in the cache session. Please login'
            error = msg
            response = ResponseCapsule(
                http_code=status_code,
                error=error,
                message=msg,
                **api_meta_to_rc_dic(api_meta=self.api_meta)
            )
            raise HTTPException(status_code=response.http_code, detail=response.model_dump())
        
        else:
            session_validity = await is_session_valid(session_id=str(session_infos.session.id))
            if session_validity:
                # Continue request normally
                # return await call_next(request)
                return
            else:
                # redirect_url = get_app_route(
                #     api_type=ApiTypeEnum.main,
                #     route='logout'
                #     )
                # return RedirectResponse(url=redirect_url)
                status_code = status.HTTP_401_UNAUTHORIZED
                msg = 'Session has expired. Please login'
                error = msg
                response = ResponseCapsule(
                    http_code=status_code,
                    error=error,
                    message=msg,
                    **api_meta_to_rc_dic(api_meta=api_meta)
                )
                raise HTTPException(status_code=response.http_code, detail=response.model_dump())



# class AuthMiddleware(BaseHTTPMiddleware):
#     """
#     Middleware for handling authentication.

#     Checks for authentication on all routes except login, logout, and authentication-related routes.
#     If authentication fails, redirects to the login page.
#     """
    
#     # def __init__(self, app: ASGIApp):
#     #     super().__init__(app)
    
#     async def dispatch(self, request: Request, call_next):
# # async def auth_middleware_dispatch(request: Request, call_next):
#         authorized_urls = [
#             # not_auth_router.url_path_for('route_login'),
#             # not_auth_router.url_path_for('route_logout'),
#             # not_auth_router.url_path_for('route_auth'),
#             # not_auth_router.url_path_for('route_auth_callback'),
#             # not_auth_router.url_path_for('route_read_ping'),            
#         ]
#         logger.debug("URL called: %s", request.url)
#         if request.url in authorized_urls:
#             logger.debug('URL authorized')
#             return await call_next(request)
        
#         # Try getting the session_id
#         session_id = request.session.get('session_id')
#         if not session_id:
#             logger.debug("Session id not found in the session")
#             encoded_session_id = request.cookies.get('session_id')
#             if encoded_session_id:
#                 logger.debug("Session id found in the cookies")
#                 # Décoder en base64
#                 session_id_json_bytes = base64.b64decode(encoded_session_id)

#                 # Convertir les octets en chaîne JSON
#                 session_id_json = session_id_json_bytes.decode('utf-8')

#                 # Charger la chaîne JSON en objet Python (dictionnaire, liste, etc.)
#                 session_id = json.loads(session_id_json)
#                 request.session['session_id'] = session_id
#             else:
#                 logger.debug("Session id not found (cookies + session)")
#                 redirect_url = get_app_route(
#                     api_type=ApiTypeEnum.main,
#                     route='login'
#                     )
#                 return RedirectResponse(redirect_url)
        
#         # Try getting the user_email
#         user_email = request.session.get('user_email')
#         if not user_email:
#             logger.debug("User email not found in the session")
#             encoded_user_email = request.cookies.get('user_email')
#             if encoded_user_email:
#                 logger.debug("User email found in the cookies")
#                 # Décoder en base64
#                 user_email_json_bytes = base64.b64decode(encoded_user_email)

#                 # Convertir les octets en chaîne JSON
#                 user_email_json = user_email_json_bytes.decode('utf-8')

#                 # Charger la chaîne JSON en objet Python (dictionnaire, liste, etc.)
#                 user_email = json.loads(user_email_json)
#                 request.session['user_email'] = user_email
#             else:
#                 logger.debug("User email not found (cookies + session)")
#                 redirect_url = get_app_route(
#                     api_type=ApiTypeEnum.main,
#                     route='login'
#                     )
#                 return RedirectResponse(redirect_url)
#         logger.debug(f"session id is {session_id}")
#         session_infos: ApiSessionInfo = await get_session_cache(session_id=session_id)    
#         # sessions_cache = await get_sessions_cache()
#         # print(f"{sessions_cache=}")
#         if not session_infos:
#             logger.debug("Any session was found in the cache")
#             redirect_url = get_app_route(
#                     api_type=ApiTypeEnum.main,
#                     route='login'
#                     )
#             return RedirectResponse(redirect_url)

#         # Get the token
#         token = session_infos.session.token
#         # Redirect to login page if token not found
#         if not token:
#             logger.debug("Access token not found in the cache session")
#             redirect_url = get_app_route(
#                     api_type=ApiTypeEnum.main,
#                     route='login'
#                     )
#             return RedirectResponse(redirect_url)
        
#         else:
#             # Decode token (include token expiration time checking)
#             # try:
#             #     decoded_token = keycloak_openid.decode_token(token)
#             #     logger.debug("Token successfully decoded")
#             # except Exception as e:
#             #     logger.debug("A problem occured while decoding the token, ensure that the token is valid.")
#             #     return RedirectResponse(url=not_auth_not_auth_router.url_path_for('route_logout'))
            
#             # await refresh_token(request=request)
            
#             # Verify the session validity (deleted and / or token valid)
#             session_validity = await is_session_valid(session_id=str(session_infos.session.id))
#             if session_validity:
#                 # Continue request normally
#                 return await call_next(request)
#             else:
#                 redirect_url = get_app_route(
#                     api_type=ApiTypeEnum.main,
#                     route='logout'
#                     )
#                 return RedirectResponse(url=redirect_url)
    
