from fastapi import HTTPException, Path, Query, APIRouter
from boilerplate_shared_module.schemes_pydantic.basic_schemes import ResponseCapsule
from ..metadata.api import get_api_meta
from boilerplate_shared_module.tools.api_infos import ApiMetadata
from ..handlers.users import *
import uuid


not_auth_router = APIRouter()

## USERS ##

@not_auth_router.get("/users", response_model=ResponseCapsule)
async def route_read_users(
    limit: int | None = None, 
    preloaded: bool = True, 
    db: AsyncSession = Depends(db_manager.async_get_db),
    api_meta: ApiMetadata = Depends(get_api_meta),
    ):
    """
    Route to read users.

    Parameters
    ----------
    limit : int, optional
        Maximum number of users to retrieve, by default None
    preloaded : bool, optional
        Whether to load preloaded data, by default True
    db : AsyncSession, optional
        Asynchronous session to interact with the database, by default Depends(db_manager.async_get_db)

    Returns
    -------
    ResponseCapsule
        Response capsule containing the user information.
    """
    response =  await handler_get_users(
        limit=limit, 
        preloaded=preloaded, 
        db=db,
        api_meta=api_meta
        )
    if response.http_code != 200:
        raise HTTPException(status_code=response.http_code, detail=response.model_dump())
    return response


@not_auth_router.get("/user/{user_id}", response_model=ResponseCapsule)
async def route_read_user(
    user_id: uuid.UUID, 
    preloaded: bool = True, 
    db: AsyncSession = Depends(db_manager.async_get_db),
    api_meta: ApiMetadata = Depends(get_api_meta),
    ):
    """
    Route to read a user by their ID.

    Parameters
    ----------
    user_id : uuid.UUID
        The ID of the user.
    preloaded : bool, optional
        Whether to load preloaded data, by default True
    db : AsyncSession, optional
        Asynchronous session to interact with the database, by default Depends(db_manager.async_get_db)

    Returns
    -------
    ResponseCapsule
        Response capsule containing the user information.
    """
    response = await handler_get_user_by_id(
        user_id=user_id, 
        preloaded=preloaded, 
        db=db,
        api_meta=api_meta
        )
    if response.http_code != 200: 
        raise HTTPException(status_code=response.http_code, detail=response.model_dump())
    return response


@not_auth_router.get("/user", response_model=ResponseCapsule)
async def route_read_user_by_email(
    email: str, 
    preloaded: bool = True, 
    db: AsyncSession = Depends(db_manager.async_get_db),
    api_meta: ApiMetadata = Depends(get_api_meta),
    ):
    """
    Route to read a user by their email.

    Parameters
    ----------
    email : str
        The email address of the user.
    preloaded : bool, optional
        Whether to load preloaded data, by default True
    db : AsyncSession, optional
        Asynchronous session to interact with the database, by default Depends(db_manager.async_get_db)

    Returns
    -------
    ResponseCapsule
        Response capsule containing the user information.
    """
    response = await handler_get_user_by_email(
        email=email, 
        preloaded=preloaded, 
        db=db,
        api_meta=api_meta
        )
    if response.http_code != 200: 
        raise HTTPException(status_code=response.http_code, detail=response.model_dump())
    return response


@not_auth_router.post("/user", response_model=ResponseCapsule, status_code=201)
async def route_create_new_user(
    user: ApiUser, 
    db: AsyncSession = Depends(db_manager.async_get_db),
    api_meta: ApiMetadata = Depends(get_api_meta),
    ):
    """
    Route to create a new user.

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
    response =  await handler_post_user(
        user=user, 
        db=db,
        api_meta=api_meta
        )
    if response.http_code != 200:
        raise HTTPException(status_code=response.http_code, detail=response.model_dump())
    return response


@not_auth_router.put("/user/{user_id}", response_model=ResponseCapsule, status_code=201)
async def route_put_user(
    user: ApiUser, 
    user_id: uuid.UUID, 
    db: AsyncSession = Depends(db_manager.async_get_db),
    api_meta: ApiMetadata = Depends(get_api_meta),
    ):
    """
    Route to update a user by their ID.

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
    response =  await handler_put_user(
        user=user, 
        user_id=user_id, 
        db=db,
        api_meta=api_meta
        )
    if response.http_code != 200:
        raise HTTPException(status_code=response.http_code, detail=response.model_dump())
    return response


@not_auth_router.delete("/user/{user_id}", response_model=ResponseCapsule, status_code=201)
async def route_soft_delete_user(
    user_id: uuid.UUID, 
    db: AsyncSession = Depends(db_manager.async_get_db),
    api_meta: ApiMetadata = Depends(get_api_meta),
    ):
    """
    Route to soft delete a user by their ID.

    Parameters
    ----------
    user_id : uuid.UUID
        The ID of the user to be soft deleted.
    db : AsyncSession, optional
        Asynchronous session to interact with the database, by default Depends(db_manager.async_get_db)

    Returns
    -------
    ResponseCapsule
        Response capsule containing the result of the user soft deletion.
    """
    response =  await handler_soft_delete_user(
        user_id=user_id, 
        db=db,
        api_meta=api_meta
        )
    if response.http_code != 200:
        # raise HTTPException(status_code=response.http_code, detail=response.model_dump())
        raise HTTPException(status_code=response.http_code, detail=response.model_dump())
    return response
