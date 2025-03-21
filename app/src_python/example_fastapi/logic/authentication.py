from starlette.requests import Request
from datetime import datetime, timezone
from fastapi import status, HTTPException

from boilerplate_shared_module.schemes_pydantic.orm import ApiUser, ApiSessionLog, ApiUserInfos
from boilerplate_shared_module.managers.authentication import KeycloakManager
from boilerplate_shared_module.managers.redis import RedisManager
from boilerplate_shared_module.managers.cache_auth_tools import *
from boilerplate_shared_module.consumers.tools import get_app_route

from .proxy_tools import *
from ..metadata.api import ApiTypeEnum
from .orm import *
from .storage.cache import manager as redis_management
from .cache import *


route_auth_callback = get_app_route(
    route="/auth_callback",
    app_meta=api_meta,
)

keycloak_manager = KeycloakManager(
    keycloak_config_path='/home/serviceuser/service/dev-secrets/ServiceClient.json',
    )
redis_manager: redis_management.RedisManager = redis_management.redis_manager
        
        
async def logout(
    request: Request, 
    access_token: str,
    ):
    session_infos: ApiSessionInfo = await redis_manager.get_session_cache(access_token=access_token)
    if session_infos:
        try:
            access_token = session_infos.session.token
            api_user_infos: ApiUserInfos = await get_and_update_user_infos(
                keycloak_manager=keycloak_manager,
                redis_manager=redis_manager,
                update_in_cache=False,
                is_user_action=False, 
                access_token=access_token,
                refresh_token=session_infos.session.refresh_token,
                )
            refresh_token = api_user_infos.session.refresh_token
            keycloak_manager.open_id.logout(refresh_token)
            # exp_token_date = datetime.now(tz=timezone.utc)
            
            # api_user_infos.session.token_expiration_date = exp_token_date
            
            # async with db_manager.AsyncSessionLocal() as db:
            #     await add_or_refresh_user_in_cache(api_user=api_user_infos.user, db=db)
                # await add_or_refresh_user_infos_in_cache_and_db(user_infos=api_user_infos,
                #                                             db=db)
            # Soft delete the session + set the update the token expiration date
            await redis_manager.delete_session_in_cache(api_session=api_user_infos.session,
                                        api_user=api_user_infos.user)
        except Exception as e:
            pass
    
    # Clear http session and cookies
    request.session.pop('access_token', None)
    request.cookies.pop('access_token', None)


async def add_token(
    request: Request, 
    code: str, db: AsyncSession | None = None
    ):
    # Utiliser le code d'autorisation pour obtenir le jeton d'accès
    authorization = keycloak_manager.create_authorization_token(
        code=code, 
        redirect_url=route_auth_callback,
        )
    # Traiter la réponse de Keycloak
    api_user_infos = await get_and_update_user_infos(
        keycloak_manager=keycloak_manager,
        redis_manager=redis_manager,
        update_in_cache=False,
        is_user_action=True,
        authorization=authorization
        )
    # request.session['user_email'] = api_user_infos.user.email
    
    # async with db_manager.AsyncSessionLocal() as db:
    new_api_user_infos = await add_or_refresh_user_infos_in_cache_and_db(
        user_infos=api_user_infos,
        db=db,
        )
    # request.session['access_token'] = 
    # Store the access token
    request.session['access_token'] = str(new_api_user_infos.session.token)
    
    # Store the user_email
    # request.session['user_email'] = str(new_api_user_infos.user.email)
    # Store the authorization
    # request.session['authorization'] = authorization
    return authorization
    

async def refresh_token(
    request: Request, 
    db: AsyncSession | None = None
    ):
    """
    Refreshes the token in the request and in the cache.

    Parameters
    ----------
    request : Request
        Incoming request object.
    """
    access_token = request.session.get('access_token')
    try:
        session_infos: ApiSessionInfo = await redis_manager.get_session_cache(access_token=access_token)
        refresh_token = session_infos.session.refresh_token
        authorization = keycloak_manager.refresh_authorization_token(refresh_token=refresh_token)
        
        new_session_fields = {
            "token": authorization.get('access_token'),
            "refresh_token": authorization.get('refresh_token'),
        }
        # Update user token in the cache
        api_user_infos = await get_and_update_user_infos(
            keycloak_manager=keycloak_manager,
            redis_manager=redis_manager,
            update_in_cache=False,
            is_user_action=False,
            access_token=access_token, 
            refresh_token=refresh_token, 
            session_update_fields=new_session_fields,
            )
        # async with db_manager.AsyncSessionLocal() as db:
        new_api_user_infos: ApiUserInfos = await add_or_refresh_user_infos_in_cache_and_db(
            user_infos=api_user_infos,
            db=db,
            )

        # Store the session id
        request.session['access_token'] = str(new_api_user_infos.session.token)

        return authorization
        # Return result
        # return authorization
    
    except Exception as e:
        raise e


async def is_session_valid(
    access_token: str | None = None,
    ) -> bool:
    """
    Check if the session is valid.

    Parameters
    ----------
    access_token : str | None
        Access token
        
    Returns
    -------
    bool
        Boolean that determine if the session is valid (True if valid).
    """
    if access_token:
        try:
            session_infos: ApiSessionInfo = await redis_manager.get_session_cache(access_token=access_token)
            
            # Session isn't valid if it has been deleted
            if session_infos.session.deleted_at != None:
                return False
            # Check the token expiration date
            else:
                token_exp = session_infos.session.token_expiration_date
                if token_exp < datetime.now(tz=timezone.utc):
                    return False
                else:
                    return True        
        except Exception as e:
            return False
    else:
        return False