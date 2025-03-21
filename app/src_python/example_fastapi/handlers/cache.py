from sqlalchemy.ext.asyncio import AsyncSession
from fastapi import HTTPException, Path, Query, Depends, Request
from ..logic.orm import *
from ..logic.cache import *
from boilerplate_shared_module.schemes_pydantic.orm import *
from boilerplate_shared_module.schemes_pydantic.basic_schemes import *
from boilerplate_shared_module.tools.models import *
from ..metadata.api import get_api_meta
from boilerplate_shared_module.tools.api_infos import ApiMetadata, api_meta_to_rc_dic
from pydantic import ValidationError
from pydantic import TypeAdapter
import traceback
from ..logic.storage.cache import manager as redis_management
from ..logic.storage.database import manager as db_management


redis_manager: redis_management.RedisManager = redis_management.redis_manager
db_manager: db_management.DatabaseManager = db_management.db_manager


## USERS ##

async def handler_get_users_cache(
    api_meta: ApiMetadata = Depends(get_api_meta),
                            ):
    """
    Handles the route to get users.

    Parameters
    ----------


    Returns
    -------
    ResponseCapsule
        Response capsule containing the requested users or error message.
    """
    api_meta_dic = api_meta_to_rc_dic(api_meta=api_meta)
    try:
        objects = await redis_manager.get_users_cache()
        api_model = ApiUser
    except Exception as e:
        msg = f'A problem occured while reading users in the cache : {e}'
        error_traceback = str(traceback.format_exc())
        return ResponseCapsule(http_code=500, error=error_traceback, message=msg, **api_meta_dic)
    try:
        ta = TypeAdapter(dict[str, api_model])
        objects = ta.validate_python(objects)
        typee = get_api_model_name(api_model=api_model)
        return ResponseCapsule(data=objects, type=typee)
    except ValidationError as e:
        msg = f'A problem occured during the validation of the model : {e}'
        error_traceback = str(traceback.format_exc())
        return ResponseCapsule(http_code=422, error=error_traceback, message=msg, **api_meta_dic)


async def handler_get_user_cache_by_id(
    user_id: uuid.UUID, 
    api_meta: ApiMetadata = Depends(get_api_meta),
    ) -> ResponseCapsule:
    """
    Handles the route to get a user by ID.

    Parameters
    ----------
    user_id : uuid.UUID
        ID of the user.


    Returns
    -------
    ResponseCapsule
        Response capsule containing the requested user or error message.
    """
    api_meta_dic = api_meta_to_rc_dic(api_meta=api_meta)
    try:
        obj = redis_manager.get_user_cache_by_id(user_id=user_id)
        api_model = ApiUser
    except Exception as e:
        msg = f'A problem occured while reading users cache in the cache : {e}'
        error_traceback = str(traceback.format_exc())
        return ResponseCapsule(http_code=500, error=error_traceback, message=msg, **api_meta_dic)
    if not obj:
        msg = "Missing user with user_id given"
        return ResponseCapsule(http_code=400, error=msg, message=msg, **api_meta_dic)
    try:
        obj = api_model.model_validate(obj)
        typee = get_api_model_name(api_model=api_model)
        return ResponseCapsule(data=obj, type=typee)
    except ValidationError as e:
        msg = f'A problem occured during the validation of the model : {e}'
        error_traceback = str(traceback.format_exc())
        return ResponseCapsule(http_code=422, error=error_traceback, message=msg, **api_meta_dic)


async def handler_get_user_cache_by_email(
    email: str, 
    api_meta: ApiMetadata = Depends(get_api_meta),
    ) -> ResponseCapsule:
    """
    Handles the route to get a user by email.

    Parameters
    ----------
    email : str
        Email of the user.

    Returns
    -------
    ResponseCapsule
        Response capsule containing the requested user or error message.
    """
    api_meta_dic = api_meta_to_rc_dic(api_meta=api_meta)
    try:
        obj = await redis_manager.get_user_cache(user_email=email)
        api_model = ApiUser
    except Exception as e:
        msg = f'A problem occured while reading users in the cache : {e}'
        error_traceback = str(traceback.format_exc())
        return ResponseCapsule(http_code=500, error=error_traceback, message=msg, **api_meta_dic)
    if not obj:
        msg = "Missing user with email given"
        return ResponseCapsule(http_code=400, error=msg, message=msg, **api_meta_dic)
    try:
        obj = api_model.model_validate(obj)
        typee = get_api_model_name(api_model=api_model)
        return ResponseCapsule(data=obj, type=typee)
    except ValidationError as e:
        msg = f'A problem occured during the validation of the model : {e}'
        error_traceback = str(traceback.format_exc())
        return ResponseCapsule(http_code=422, error=error_traceback, message=msg, **api_meta_dic)



async def handler_post_user_cache(user: ApiUser, 
    db: AsyncSession = Depends(db_manager.async_get_db), 
    api_meta: ApiMetadata = Depends(get_api_meta),
    ) -> ResponseCapsule:
    """
    Creates a new user in the cache.

    Parameters
    ----------
    user : ApiUser
        The user data to create.
    db : AsyncSession, optional
        Database session, by default Depends(db_manager.async_get_db)

    Returns
    -------
    ResponseCapsule
        Response capsule containing the created user or error message.
    """
    api_meta_dic = api_meta_to_rc_dic(api_meta=api_meta)
    try:
        user = await redis_manager.get_user_cache(user_email=user.email)
        api_model = ApiUser
    except Exception as e:
        msg = f'A problem occured while reading user in the cache : {e}'
        error_traceback = str(traceback.format_exc())
        return ResponseCapsule(http_code=500, error=error_traceback, message=msg, **api_meta_dic)
    if user:
        msg = "Email is already registered in the cache"
        return ResponseCapsule(http_code=400, error=msg, message=msg, **api_meta_dic)
    try:
        data = await add_or_refresh_user_in_cache(
            db=db, 
            api_user=user,
            )
    except Exception as e:
        msg = f'A problem occured while adding the user in the cache : {e}'
        error_traceback = str(traceback.format_exc())
        return ResponseCapsule(http_code=500, error=error_traceback, message=msg, **api_meta_dic)
    try:
        data = api_model.model_validate(data)
        typee = get_api_model_name(api_model=api_model)
        return ResponseCapsule(data=data, type=typee, **api_meta_dic)
    except ValidationError as e:
        msg = f'A problem occured during the validation of the model : {e}'
        error_traceback = str(traceback.format_exc())
        return ResponseCapsule(http_code=422, error=error_traceback, message=msg, **api_meta_dic)
    
    
async def handler_put_user_cache(user: ApiUser, user_id: uuid.UUID, 
    db: AsyncSession = Depends(db_manager.async_get_db), 
    api_meta: ApiMetadata = Depends(get_api_meta),
    ) -> ResponseCapsule:
    """
    Updates an existing user in the cache.

    Parameters
    ----------
    user : ApiUser
        The updated user data.
    user_id : uuid.UUID
        The ID of the user to update.
    db : AsyncSession, optional
        Database session, by default Depends(db_manager.async_get_db)

    Returns
    -------
    ResponseCapsule
        Response capsule containing the updated user or error message.
    """
    api_meta_dic = api_meta_to_rc_dic(api_meta=api_meta)
    try:
        user = await redis_manager.get_user_cache_by_id(
            user_id=user_id,
            )
        api_model = ApiUser
    except Exception as e:
        msg = f'A problem occured while reading user in the cache : {e}'
        error_traceback = str(traceback.format_exc())
        return ResponseCapsule(http_code=500, error=error_traceback, message=msg, **api_meta_dic)
    if not user:
        msg = "Missing user in cache with user_id given"
        return ResponseCapsule(http_code=400, error=msg, message=msg, **api_meta_dic)
    try:
        data = await add_or_refresh_user_in_cache(api_user=user,
                                                db=db,
                                                )
    except Exception as e:
        msg = f'A problem occured while updating the user in the cache : {e}'
        error_traceback = str(traceback.format_exc())
        return ResponseCapsule(http_code=500, error=error_traceback, message=msg, **api_meta_dic)

    try:
        data = api_model.model_validate(data)
        typee = get_api_model_name(api_model=api_model)
        
        return ResponseCapsule(data=data, type=typee, **api_meta_dic)
    except ValidationError as e:
        msg = f'A problem occured during the validation of the model : {e}'
        error_traceback = str(traceback.format_exc())
        return ResponseCapsule(http_code=422, error=error_traceback, message=msg, **api_meta_dic)
                
    
async def handler_soft_delete_user_cache(
    user_id: uuid.UUID,
    api_meta: ApiMetadata = Depends(get_api_meta),
    ) -> ResponseCapsule:
    """
    Soft deletes a user from the cache.

    Parameters
    ----------
    user_id : uuid.UUID
        The ID of the user to soft delete.
    edis_client: Redis
        Redis Client

    Returns
    -------
    ResponseCapsule
        Response capsule containing the soft deleted user or error message.
    """
    api_meta_dic = api_meta_to_rc_dic(api_meta=api_meta)
    try:
        user = await redis_manager.get_user_cache_by_id(user_id=user_id)
        api_model = ApiUser
    except Exception as e:
        
        msg = f'A problem occured while reading user in the cache : {e}'
        error_traceback = str(traceback.format_exc())
        return ResponseCapsule(http_code=500, error=error_traceback, message=msg, **api_meta_dic)
    if not user:
        msg = "Missing user with user_id given"
        return ResponseCapsule(http_code=400, error=msg, message=msg, **api_meta_dic)

    try:
        data = await redis_manager.delete_user_in_cache(
            api_user=user
            )
    except Exception as e:
        msg = f'A problem occured while deleting the user in the cache : {e}'
        error_traceback = str(traceback.format_exc())
        return ResponseCapsule(http_code=500, error=error_traceback, message=msg, **api_meta_dic)

    try:
        data = api_model.model_validate(data)
        typee = get_api_model_name(api_model=api_model)
        
        return ResponseCapsule(data=data, type=typee, **api_meta_dic)
    except ValidationError as e:
        msg = f'A problem occured during the validation of the model : {e}'
        error_traceback = str(traceback.format_exc())
        return ResponseCapsule(http_code=422, error=error_traceback, message=msg, **api_meta_dic)
