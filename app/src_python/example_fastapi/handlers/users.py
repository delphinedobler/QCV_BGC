from sqlalchemy.ext.asyncio import AsyncSession
from fastapi import HTTPException, Path, Query, Depends, Request
from ..logic.orm import *
from boilerplate_shared_module.schemes_pydantic.orm import *
from boilerplate_shared_module.schemes_pydantic.basic_schemes import *
from boilerplate_shared_module.tools.models import *
from ..metadata.api import get_api_meta
from boilerplate_shared_module.tools.api_infos import ApiMetadata, api_meta_to_rc_dic
from ..logic.cache import *
from pydantic import ValidationError
from pydantic import TypeAdapter
import traceback
from starlette import status
from ..logic.storage.cache import manager as redis_management
from ..logic.storage.database import manager as db_management


redis_manager: redis_management.RedisManager = redis_management.redis_manager
db_manager: db_management.DatabaseManager = db_management.db_manager


## USER MIDDLEWARE ##

async def get_authenticated_user(
    request: Request,
    api_meta: ApiMetadata = Depends(get_api_meta),
    ) -> ApiUser:
    """
    Middleware to retrieve the authenticated user from the session.

    Parameters
    ----------
    request : Request
        Incoming request.

    Returns
    -------
    ApiUser
        Authenticated user.
    """
    access_token = request.session.get('access_token')
    if not access_token:
        status_code = status.HTTP_401_UNAUTHORIZED
        msg = 'No user found in the session. Please login.'
        error = msg
        response = ResponseCapsule(
            http_code=status_code,
            error=error,
            message=msg,
            **api_meta_to_rc_dic(api_meta=api_meta)
        )
        raise HTTPException(status_code=response.http_code, detail=response.model_dump())
    session_infos = await redis_manager.get_session_cache(access_token=access_token)
    api_user = await redis_manager.get_user_cache(user_email=session_infos.user_email)
    if not api_user:
        status_code = status.HTTP_401_UNAUTHORIZED
        msg = 'User not found. Please login'
        error = msg
        response = ResponseCapsule(
            http_code=status_code,
            error=error,
            message=msg,
            **api_meta_to_rc_dic(api_meta=api_meta)
        )
        raise HTTPException(status_code=response.http_code, detail=response.model_dump())
    return api_user


## USERS ##

async def handler_get_users(
    preloaded: bool, 
    limit: int | None = None, 
    db: AsyncSession = Depends(db_manager.async_get_db), 
    api_meta: ApiMetadata = Depends(get_api_meta),
    ):
    """
    Handles the route to get users.

    Parameters
    ----------
    preloaded : bool
        Flag indicating if preloaded users are requested.
    limit : int, optional
        Maximum number of users to retrieve, by default None
    db : AsyncSession, optional
        Database session, by default Depends(db_manager.async_get_db)

    Returns
    -------
    ResponseCapsule
        Response capsule containing the requested users or error message.
    """
    api_meta_dic = api_meta_to_rc_dic(api_meta=api_meta)
    try:
        if preloaded:
            db_objects = await get_preloaded_users(db, limit=limit)
            api_model = ApiPreloadedUser
        else:
            db_objects = await get_users(db, limit=limit)
            api_model = ApiUser
    except Exception as e:
        msg = f'A problem occured while reading users in database : {e}'
        error_traceback = str(traceback.format_exc())
        return ResponseCapsule(http_code=500, error=error_traceback, **api_meta_dic)
    try:
        ta = TypeAdapter(list[api_model])
        db_objects = ta.validate_python(db_objects)
        typee = get_api_model_name(api_model=api_model)
        return ResponseCapsule(data=db_objects, type=typee)
    except ValidationError as e:
        msg = f'A problem occured during the validation of the orm model : {e}'
        error_traceback = str(traceback.format_exc())
        return ResponseCapsule(http_code=422, error=error_traceback, **api_meta_dic)


async def handler_get_user_by_id(
    user_id: uuid.UUID, 
    preloaded: bool, 
    db: AsyncSession = Depends(db_manager.async_get_db),
    api_meta: ApiMetadata = Depends(get_api_meta),
    ) -> ResponseCapsule:
    """
    Handles the route to get a user by ID.

    Parameters
    ----------
    user_id : uuid.UUID
        ID of the user.
    preloaded : bool
        Flag indicating if a preloaded user is requested.
    db : AsyncSession, optional
        Database session, by default Depends(db_manager.async_get_db)

    Returns
    -------
    ResponseCapsule
        Response capsule containing the requested user or error message.
    """
    api_meta_dic = api_meta_to_rc_dic(api_meta=api_meta)
    try:
        if preloaded:
            db_obj = await get_preloaded_user(db=db, user_id=user_id)
            api_model = ApiPreloadedUser
        else:
            db_obj = await get_user(db=db, user_id=user_id)
            api_model = ApiUser
    except Exception as e:
        msg = f'A problem occured while reading users in database : {e}'
        error_traceback = str(traceback.format_exc())
        return ResponseCapsule(http_code=500, error=error_traceback, **api_meta_dic)
    if not db_obj:
        msg = "Missing user with user_id given"
        return ResponseCapsule(http_code=400, error=msg, **api_meta_dic)
    try:
        db_obj = api_model.model_validate(db_obj)
        typee = get_api_model_name(api_model=api_model)
        return ResponseCapsule(data=db_obj, type=typee)
    except ValidationError as e:
        msg = f'A problem occured during the validation of the orm model : {e}'
        error_traceback = str(traceback.format_exc())
        return ResponseCapsule(http_code=422, error=error_traceback, **api_meta_dic)


async def handler_get_user_by_email(
    email: str, 
    preloaded: bool, 
    db: AsyncSession = Depends(db_manager.async_get_db),
    api_meta: ApiMetadata = Depends(get_api_meta),
    ) -> ResponseCapsule:
    """
    Handles the route to get a user by email.

    Parameters
    ----------
    email : str
        Email of the user.
    preloaded : bool
        Flag indicating if a preloaded user is requested.
    db : AsyncSession, optional
        Database session, by default Depends(db_manager.async_get_db)

    Returns
    -------
    ResponseCapsule
        Response capsule containing the requested user or error message.
    """
    api_meta_dic = api_meta_to_rc_dic(api_meta=api_meta)
    try:
        if preloaded:
            db_obj = await get_preloaded_user_by_email(db=db, email=email)
            api_model = ApiPreloadedUser
        else:
            db_obj = await get_user_by_email(db=db, email=email)
            api_model = ApiUser
    except Exception as e:
        msg = f'A problem occured while reading users in database : {e}'
        error_traceback = str(traceback.format_exc())
        return ResponseCapsule(http_code=500, error=error_traceback, **api_meta_dic)
    if not db_obj:
        msg = "Missing user with email given"
        return ResponseCapsule(http_code=400, error=msg, **api_meta_dic)
    try:
        db_obj = api_model.model_validate(db_obj)
        typee = get_api_model_name(api_model=api_model)
        return ResponseCapsule(data=db_obj, type=typee)
    except ValidationError as e:
        msg = f'A problem occured during the validation of the orm model : {e}'
        error_traceback = str(traceback.format_exc())
        return ResponseCapsule(http_code=422, error=error_traceback, **api_meta_dic)



async def handler_post_user(
    user: ApiUser, 
    db: AsyncSession = Depends(db_manager.async_get_db),
    api_meta: ApiMetadata = Depends(get_api_meta),
    ) -> ResponseCapsule:
    """
    Creates a new user in the database.

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
        db_user = await get_user_by_email(db=db, email=user.email)
        api_model = ApiUser
    except Exception as e:
        msg = f'A problem occured while reading user in database : {e}'
        error_traceback = str(traceback.format_exc())
        return ResponseCapsule(http_code=500, error=error_traceback, **api_meta_dic)
    if db_user:
        msg = "Email is already registered"
        return ResponseCapsule(http_code=400, error=msg, **api_meta_dic)
    db_user = User(**dict(user))
    data = await create_user(db=db, db_user=db_user)
    try:
        data = api_model.model_validate(data)
        typee = get_api_model_name(api_model=api_model)
        
        # Add user in the cache
        # cache_users.USERS_CACHE[data.email] = data
        
        return ResponseCapsule(data=data, type=typee)
    except ValidationError as e:
        msg = f'A problem occured during the validation of the orm model : {e}'
        error_traceback = str(traceback.format_exc())
        return ResponseCapsule(http_code=422, error=error_traceback, **api_meta_dic)
    
    
async def handler_put_user(
    user: ApiUser, 
    user_id: uuid.UUID, 
    db: AsyncSession = Depends(db_manager.async_get_db),
    api_meta: ApiMetadata = Depends(get_api_meta),
    ) -> ResponseCapsule:
    """
    Updates an existing user in the database.

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
        db_user = await get_user(db=db, user_id=user_id)
        api_model = ApiUser
    except Exception as e:
        
        msg = f'A problem occured while reading user in database : {e}'
        error_traceback = str(traceback.format_exc())
        return ResponseCapsule(http_code=500, error=error_traceback, **api_meta_dic)
    if not db_user:
        msg = "Missing user with user_id given"
        return ResponseCapsule(http_code=400, error=msg, **api_meta_dic)
    for key, val in dict(user).items():
        # replace attributes that have values
        if val:
            setattr(db_user, key, val)
    setattr(db_user, 'id', user_id)
    data = await update_user(db=db, db_user=db_user)
    try:
        data = api_model.model_validate(data)
        typee = get_api_model_name(api_model=api_model)
        
        # Update user in the cache
        # cache_users.USERS_CACHE[data.email] = data
        
        return ResponseCapsule(data=data, type=typee)
    except ValidationError as e:
        msg = f'A problem occured during the validation of the orm model : {e}'
        error_traceback = str(traceback.format_exc())
        return ResponseCapsule(http_code=422, error=error_traceback, **api_meta_dic)
                
    
async def handler_soft_delete_user(
    user_id: uuid.UUID, 
    db: AsyncSession = Depends(db_manager.async_get_db),
    api_meta: ApiMetadata = Depends(get_api_meta),
    ) -> ResponseCapsule:
    """
    Soft deletes a user from the database.

    Parameters
    ----------
    user_id : uuid.UUID
        The ID of the user to soft delete.
    db : AsyncSession, optional
        Database session, by default Depends(db_manager.async_get_db)

    Returns
    -------
    ResponseCapsule
        Response capsule containing the soft deleted user or error message.
    """
    api_meta_dic = api_meta_to_rc_dic(api_meta=api_meta)
    try:
        db_user = await get_user(db=db, user_id=user_id)
        api_model = ApiUser
    except Exception as e:
        
        msg = f'A problem occured while reading user in database : {e}'
        error_traceback = str(traceback.format_exc())
        return ResponseCapsule(http_code=500, error=error_traceback, **api_meta_dic)
    if not db_user:
        msg = "Missing user with user_id given"
        return ResponseCapsule(http_code=400, error=msg, **api_meta_dic)

    data = await soft_delete_user(db=db, db_user=db_user)
    try:
        data = api_model.model_validate(data)
        typee = get_api_model_name(api_model=api_model)
        
        # Remove user from the cache
        # cache_users.USERS_CACHE.pop(data.email)
        
        return ResponseCapsule(data=data, type=typee)
    except ValidationError as e:
        msg = f'A problem occured during the validation of the orm model : {e}'
        error_traceback = str(traceback.format_exc())
        return ResponseCapsule(http_code=422, error=error_traceback, **api_meta_dic)
