from fastapi import HTTPException, Path, Query, APIRouter
from boilerplate_shared_module.schemes_pydantic.basic_schemes import ResponseCapsule
from ..metadata.api import get_api_meta
from boilerplate_shared_module.tools.api_infos import ApiMetadata
from ..handlers.cache import *
import uuid
from redis.asyncio import Redis


router = APIRouter()


## USERS ##

@router.get("/users_cache", response_model=ResponseCapsule)
async def route_read_users_cache(
    api_meta: ApiMetadata = Depends(get_api_meta)
    ):
    """
    Route to read users cache.

    Parameters
    ----------

    Returns
    -------
    ResponseCapsule
        Response capsule containing the user information.
    """
    response =  await handler_get_users_cache(api_meta=api_meta)
    if response.http_code != 200:
        raise HTTPException(status_code=response.http_code, detail=response.model_dump())
    return response


@router.get("/user_cache/{user_id}", response_model=ResponseCapsule)
async def route_read_user(user_id: uuid.UUID, 
                           api_meta: ApiMetadata = Depends(get_api_meta)
                          ):
    """
    Route to read a user in the cache by their ID.

    Parameters
    ----------
    user_id : uuid.UUID
        The ID of the user.
    db : AsyncSession, optional
        Asynchronous session to interact with the database, by default Depends(db_manager.async_get_db)

    Returns
    -------
    ResponseCapsule
        Response capsule containing the user information.
    """
    response = await handler_get_user_cache_by_id(user_id=user_id, 
                                                  api_meta=api_meta)
    if response.http_code != 200: 
        raise HTTPException(status_code=response.http_code, detail=response.model_dump())
    return response


@router.get("/user_cache", response_model=ResponseCapsule)
async def route_read_user_by_email(email: str, 
                                     api_meta: ApiMetadata = Depends(get_api_meta)
                                   ):
    """
    Route to read a user in the cache by their email.

    Parameters
    ----------
    email : str
        The email address of the user.

    Returns
    -------
    ResponseCapsule
        Response capsule containing the user information.
    """
    response = await handler_get_user_cache_by_email(email=email, 
                                                     api_meta=api_meta)
    if response.http_code != 200: 
        raise HTTPException(status_code=response.http_code, detail=response.model_dump())
    return response


@router.post("/user_cache", response_model=ResponseCapsule, status_code=201)
async def route_create_new_user(user: ApiUser, 
                                db: AsyncSession = Depends(db_manager.async_get_db),
                                       api_meta: ApiMetadata = Depends(get_api_meta)
                                ):
    """
    Route to create a new user in the cache.

    Parameters
    ----------
    user : ApiUser
        The user data to be created.
    db : AsyncSession, optional
        Asynchronous session to interact with the database, by default Depends(db_manager.async_get_db)

    Returns
    -------
    ResponseCapsule
        Response capsule containing the result of the user creation.
    """
    response =  await handler_post_user_cache(user=user, db=db, 
                                              api_meta=api_meta)
    if response.http_code != 200:
        raise HTTPException(status_code=response.http_code, detail=response.model_dump())
    return response


@router.put("/user_cache/{user_id}", response_model=ResponseCapsule, status_code=201)
async def route_put_user(user: ApiUser, 
                         user_id: uuid.UUID, 
                         db: AsyncSession = Depends(db_manager.async_get_db),
                         api_meta: ApiMetadata = Depends(get_api_meta)
                         ):
    """
    Route to update a user in the cache by their ID.

    Parameters
    ----------
    user : ApiUser
        The updated user data.
    user_id : uuid.UUID
        The ID of the user to be updated.
    db : AsyncSession, optional
        Asynchronous session to interact with the database, by default Depends(db_manager.async_get_db)
    
    Returns
    -------
    ResponseCapsule
        Response capsule containing the result of the user update.
    """
    response =  await handler_put_user_cache(user=user, 
                                             user_id=user_id, 
                                             db=db,
                                             api_meta=api_meta)
    if response.http_code != 200:
        raise HTTPException(status_code=response.http_code, detail=response.model_dump())
    return response


@router.delete("/user_cache/{user_id}", response_model=ResponseCapsule, status_code=201)
async def route_soft_delete_user(user_id: uuid.UUID, 
                                 api_meta: ApiMetadata = Depends(get_api_meta)
                                 ):
    """
    Route to soft delete a user in the cache by their ID.

    Parameters
    ----------
    user_id : uuid.UUID
        The ID of the user to be soft deleted.

    Returns
    -------
    ResponseCapsule
        Response capsule containing the result of the user soft deletion.
    """
    response =  await handler_soft_delete_user_cache(user_id=user_id, 
                                                     api_meta=api_meta)
    if response.http_code != 200:
        # raise HTTPException(status_code=response.http_code, detail=response.model_dump())
        raise HTTPException(status_code=response.http_code, detail=response.model_dump())
    return response
