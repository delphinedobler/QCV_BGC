from fastapi import HTTPException, Path, Query, APIRouter
from ..metadata.api import get_api_meta
from boilerplate_shared_module.tools.api_infos import ApiMetadata
from ..handlers.objects import *
import uuid
from .users import get_authenticated_user


not_auth_router = APIRouter()

## OBJECTS ##

@not_auth_router.get("/objects", response_model=ResponseCapsule)
async def route_read_objects(
    request: Request,
    limit: int | None = None, 
    preloaded: bool = True, 
    db: AsyncSession = Depends(db_manager.async_get_db),
    api_meta: ApiMetadata = Depends(get_api_meta),
    ):
    """
    Retrieves objects from the database via HTTP GET request.

    Parameters
    ----------
    limit : int, optional
        The maximum number of objects to retrieve, by default None.
    preloaded : bool, optional
        Whether to retrieve preloaded data, by default True.
    db : AsyncSession
        The asynchronous database session.

    Returns
    -------
    ResponseCapsule
        The response capsule containing the retrieved objects (list[ApiObject]).

    Raises
    ------
    HTTPException
        If an error occurs during the retrieval process.

    """
    response =  await handler_get_objects(
        limit=limit, 
        preloaded=preloaded, 
        db=db,
        api_meta=api_meta,
        )
    if response.http_code != 200:
        raise HTTPException(status_code=response.http_code, detail=response.model_dump())
    return response


@not_auth_router.get("/object/{object_id}", response_model=ResponseCapsule)
async def route_read_object(
    request: Request, 
    object_id: uuid.UUID, 
    preloaded: bool = True, 
    db: AsyncSession = Depends(db_manager.async_get_db),
    api_meta: ApiMetadata = Depends(get_api_meta),
    ):
    """
    Retrieves a specific object from the database by ID via HTTP GET request.

    Parameters
    ----------
    object_id : UUID
        The ID of the object to retrieve.
    preloaded : bool, optional
        Whether to retrieve preloaded data, by default True.
    db : AsyncSession
        The asynchronous database session.

    Returns
    -------
    ResponseCapsule
        The response capsule containing the retrieved object (ApiObject).

    Raises
    ------
    HTTPException
        If the object with the specified ID is not found or an error occurs during the retrieval process.

    See Also
    --------
    route_read_object, route_read_object_by_name, route_create_new_object, route_put_object, route_soft_delete_object, boilerplate_fs_api.shared_models.orm.ApiObject
    """
    response = await handler_get_object_by_id(
        object_id=object_id, 
        preloaded=preloaded, 
        db=db,
        api_meta=api_meta,
        )
    if response.http_code != 200: 
        raise HTTPException(status_code=response.http_code, detail=response.model_dump())
    return response


@not_auth_router.get("/object", response_model=ResponseCapsule)
async def route_read_object_by_name(
    name: str, 
    preloaded: bool = True, 
    db: AsyncSession = Depends(db_manager.async_get_db),
    api_meta: ApiMetadata = Depends(get_api_meta),
    ):
    """
    Retrieves a specific object from the database by name via HTTP GET request.

    Parameters
    ----------
    name : str
        The name of the object to retrieve.
    preloaded : bool, optional
        Whether to retrieve preloaded data, by default True.
    db : AsyncSession
        The asynchronous database session.

    Returns
    -------
    ResponseCapsule
        The response capsule containing the retrieved object (ApiObject).

    Raises
    ------
    HTTPException
        If the object with the specified name is not found or an error occurs during the retrieval process.

    See Also
    --------
    route_read_object, route_read_object_by_name, route_create_new_object, route_put_object, route_soft_delete_object, boilerplate_fs_api.shared_models.orm.ApiObject

    """
    response = await handler_get_object_by_name(
        name=name, 
        preloaded=preloaded, 
        db=db,
        api_meta=api_meta,
        )
    if response.http_code != 200: 
        raise HTTPException(status_code=response.http_code, detail=response.model_dump())
    return response


@not_auth_router.post("/object", response_model=ResponseCapsule, status_code=201)
async def route_create_new_object(
    object: ApiObject, 
    db: AsyncSession = Depends(db_manager.async_get_db),
    db_user: User = Depends(get_authenticated_user),
    api_meta: ApiMetadata = Depends(get_api_meta),
    ):
    """
    Creates a new object in the database via HTTP POST request.

    Parameters
    ----------
    object : ApiObject
        The object to be created.
    db : AsyncSession
        The asynchronous database session.
    db_user : User
        The authenticated user.

    Returns
    -------
    ResponseCapsule
        The response capsule containing the created object (ApiObject).

    Raises
    ------
    HTTPException
        If the creation process fails.

    See Also
    --------
    route_read_object, route_read_object_by_name, route_create_new_object, route_put_object, route_soft_delete_object, boilerplate_fs_api.shared_models.orm.ApiObject

    """
    response =  await handler_post_object(
        object=object, 
        db=db, 
        db_user=db_user,
        api_meta=api_meta,
        )
    if response.http_code != 200:
        raise HTTPException(status_code=response.http_code, detail=response.model_dump())
    return response


@not_auth_router.put("/object/{object_id}", response_model=ResponseCapsule, status_code=201)
async def route_put_object(
    object: ApiObject, 
    object_id: uuid.UUID, 
    db: AsyncSession = Depends(db_manager.async_get_db),
    db_user: User = Depends(get_authenticated_user),
    api_meta: ApiMetadata = Depends(get_api_meta),
    ):
    """
    Updates an existing object in the database via HTTP PUT request.

    Parameters
    ----------
    object : ApiObject
        The object with updated data.
    object_id : UUID
        The ID of the object to update.
    db : AsyncSession
        The asynchronous database session.
    db_user : User
        The authenticated user.

    Returns
    -------
    ResponseCapsule
        The response capsule containing the updated object (ApiObject).

    Raises
    ------
    HTTPException
        If the update process fails.

    See Also
    --------
    route_read_object, route_read_object_by_name, route_create_new_object, route_put_object, route_soft_delete_object, boilerplate_fs_api.shared_models.orm.ApiObject

    """
    response =  await handler_put_object(
        object=object, 
        object_id=object_id, 
        db=db, 
        db_user=db_user,
        api_meta=api_meta,
        )
    if response.http_code != 200:
        raise HTTPException(status_code=response.http_code, detail=response.model_dump())
    return response


@not_auth_router.delete("/object/{object_id}", response_model=ResponseCapsule, status_code=201)
async def route_soft_delete_object(
    object_id: uuid.UUID, 
    db: AsyncSession = Depends(db_manager.async_get_db),
    db_user: User = Depends(get_authenticated_user),
    api_meta: ApiMetadata = Depends(get_api_meta),
    ):
    """
    Soft deletes an object in the database via HTTP DELETE request.

    Parameters
    ----------
    object_id : UUID
        The ID of the object to delete.
    db : AsyncSession
        The asynchronous database session.
    db_user : User
        The authenticated user.

    Returns
    -------
    ResponseCapsule
        The response capsule containing the soft deleted object (ApiObject).

    Raises
    ------
    HTTPException
        If the deletion process fails.
        

    See Also
    --------
    route_read_object, route_read_object_by_name, route_create_new_object, route_put_object, route_soft_delete_object, boilerplate_fs_api.shared_models.orm.ApiObject

    """
    response =  await handler_soft_delete_object(
        object_id=object_id, 
        db=db, 
        db_user=db_user,
        api_meta=api_meta,
        )
    if response.http_code != 200:
        # raise HTTPException(status_code=response.http_code, detail=response.model_dump())
        raise HTTPException(status_code=response.http_code, detail=response.model_dump())
    return response
