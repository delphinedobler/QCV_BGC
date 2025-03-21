from fastapi import HTTPException, Path, Query, APIRouter
from boilerplate_shared_module.schemes_pydantic.basic_schemes import ResponseCapsule
from ..metadata.api import get_api_meta
from boilerplate_shared_module.tools.api_infos import ApiMetadata
from ..handlers.sessions import *
import uuid
from .users import get_authenticated_user


router = APIRouter()

## SESSIONS ##

@router.get("/session_logs", response_model=ResponseCapsule)
async def route_read_session_logs(
    limit: int | None = None, 
    preloaded: bool = True, 
    db: AsyncSession = Depends(db_manager.async_get_db),
    api_meta: ApiMetadata = Depends(get_api_meta),
    ):
    """
    Retrieves session_logs from the database via HTTP GET request.

    Parameters
    ----------
    limit : int, optional
        The maximum number of session_logs to retrieve, by default None.
    preloaded : bool, optional
        Whether to retrieve preloaded data, by default True.
    db : AsyncSession
        The asynchronous database session.

    Returns
    -------
    ResponseCapsule
        The response capsule containing the retrieved session_logs (list[ApiSessionLog]).

    Raises
    ------
    HTTPException
        If an error occurs during the retrieval process.

    """
    response =  await handler_get_session_logs(
        limit=limit, 
        preloaded=preloaded, 
        db=db,
        api_meta=api_meta
        )
    if response.http_code != 200:
        raise HTTPException(status_code=response.http_code, detail=response.model_dump())
    return response


@router.get("/session_log/{session_log_id}", response_model=ResponseCapsule)
async def route_read_session_log(
    session_log_id: uuid.UUID, 
    preloaded: bool = True, 
    db: AsyncSession = Depends(db_manager.async_get_db),
    api_meta: ApiMetadata = Depends(get_api_meta),
    ):
    """
    Retrieves a specific session_log from the database by ID via HTTP GET request.

    Parameters
    ----------
    session_log_id : UUID
        The ID of the session_log to retrieve.
    preloaded : bool, optional
        Whether to retrieve preloaded data, by default True.
    db : AsyncSession
        The asynchronous database session.

    Returns
    -------
    ResponseCapsule
        The response capsule containing the retrieved session_log (ApiSessionLog).

    Raises
    ------
    HTTPException
        If the session_log with the specified ID is not found or an error occurs during the retrieval process.

    See Also
    --------
    route_read_session_log, route_read_session_log_by_token, route_create_new_session_log, route_put_session_log, route_soft_delete_session_log, boilerplate_fs_api.shared_models.orm.ApiSessionLog
    """
    response = await handler_get_session_log_by_id(
        session_log_id=session_log_id, 
        preloaded=preloaded, 
        db=db,
        api_meta=api_meta
        )
    if response.http_code != 200: 
        raise HTTPException(status_code=response.http_code, detail=response.model_dump())
    return response


@router.get("/session_log", response_model=ResponseCapsule)
async def route_read_session_log_by_token(
    token: str, 
    preloaded: bool = True, 
    db: AsyncSession = Depends(db_manager.async_get_db),
    api_meta: ApiMetadata = Depends(get_api_meta),
    ):
    """
    Retrieves a specific session_log from the database by token via HTTP GET request.

    Parameters
    ----------
    token : str
        The token of the session_log to retrieve.
    preloaded : bool, optional
        Whether to retrieve preloaded data, by default True.
    db : AsyncSession
        The asynchronous database session.

    Returns
    -------
    ResponseCapsule
        The response capsule containing the retrieved session_log (ApiSessionLog).

    Raises
    ------
    HTTPException
        If the session_log with the specified token is not found or an error occurs during the retrieval process.

    See Also
    --------
    route_read_session_log, route_read_session_log_by_token, route_create_new_session_log, route_put_session_log, route_soft_delete_session_log, boilerplate_fs_api.shared_models.orm.ApiSessionLog

    """
    response = await handler_get_session_log_by_token(
        token=token, 
        preloaded=preloaded, 
        db=db,
        api_meta=api_meta
        )
    if response.http_code != 200: 
        raise HTTPException(status_code=response.http_code, detail=response.model_dump())
    return response


@router.post("/session_log", response_model=ResponseCapsule, status_code=201)
async def route_create_new_session_log(
    session_log: ApiSessionLog, 
    db: AsyncSession = Depends(db_manager.async_get_db),
    db_user: User = Depends(get_authenticated_user),
    api_meta: ApiMetadata = Depends(get_api_meta),
    ):
    """
    Creates a new session_log in the database via HTTP POST request.

    Parameters
    ----------
    session_log : ApiSessionLog
        The session_log to be created.
    db : AsyncSession
        The asynchronous database session.
    db_user : User
        The authenticated user.

    Returns
    -------
    ResponseCapsule
        The response capsule containing the created session_log (ApiSessionLog).

    Raises
    ------
    HTTPException
        If the creation process fails.

    See Also
    --------
    route_read_session_log, route_read_session_log_by_token, route_create_new_session_log, route_put_session_log, route_soft_delete_session_log, boilerplate_fs_api.shared_models.orm.ApiSessionLog

    """
    response =  await handler_post_session_log(
        session_log=session_log, 
        db=db, 
        db_user=db_user,
        api_meta=api_meta
        )
    if response.http_code != 200:
        raise HTTPException(status_code=response.http_code, detail=response.model_dump())
    return response


@router.put("/session_log/{session_log_id}", response_model=ResponseCapsule, status_code=201)
async def route_put_session_log(
    session_log: ApiSessionLog, 
    session_log_id: uuid.UUID, 
    db: AsyncSession = Depends(db_manager.async_get_db),
    db_user: User = Depends(get_authenticated_user),
    api_meta: ApiMetadata = Depends(get_api_meta),
    ):
    """
    Updates an existing session_log in the database via HTTP PUT request.

    Parameters
    ----------
    session_log : ApiSessionLog
        The session_log with updated data.
    session_log_id : UUID
        The ID of the session_log to update.
    db : AsyncSession
        The asynchronous database session.
    db_user : User
        The authenticated user.

    Returns
    -------
    ResponseCapsule
        The response capsule containing the updated session_log (ApiSessionLog).

    Raises
    ------
    HTTPException
        If the update process fails.

    See Also
    --------
    route_read_session_log, route_read_session_log_by_token, route_create_new_session_log, route_put_session_log, route_soft_delete_session_log, boilerplate_fs_api.shared_models.orm.ApiSessionLog

    """
    response =  await handler_put_session_log(
        session_log=session_log, 
        session_log_id=session_log_id, 
        db=db, 
        db_user=db_user,
        api_meta=api_meta
        )
    if response.http_code != 200:
        raise HTTPException(status_code=response.http_code, detail=response.model_dump())
    return response


@router.delete("/session_log/{session_log_id}", response_model=ResponseCapsule, status_code=201)
async def route_soft_delete_session_log(
    session_log_id: uuid.UUID, 
    db: AsyncSession = Depends(db_manager.async_get_db),
    db_user: User = Depends(get_authenticated_user),
    api_meta: ApiMetadata = Depends(get_api_meta),
    ):
    """
    Soft deletes an session_log in the database via HTTP DELETE request.

    Parameters
    ----------
    session_log_id : UUID
        The ID of the session_log to delete.
    db : AsyncSession
        The asynchronous database session.
    db_user : User
        The authenticated user.

    Returns
    -------
    ResponseCapsule
        The response capsule containing the soft deleted session_log (ApiSessionLog).

    Raises
    ------
    HTTPException
        If the deletion process fails.
        

    See Also
    --------
    route_read_session_log, route_read_session_log_by_token, route_create_new_session_log, route_put_session_log, route_soft_delete_session_log, boilerplate_fs_api.shared_models.orm.ApiSessionLog

    """
    response =  await handler_soft_delete_session_log(
        session_log_id=session_log_id, 
        db=db, 
        db_user=db_user,
        api_meta=api_meta
        )
    if response.http_code != 200:
        # raise HTTPException(status_code=response.http_code, detail=response.model_dump())
        raise HTTPException(status_code=response.http_code, detail=response.model_dump())
    return response
