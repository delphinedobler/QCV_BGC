from fastapi import HTTPException, Path, Query, APIRouter
from boilerplate_shared_module.schemes_pydantic.basic_schemes import ResponseCapsule
from ..metadata.api import get_api_meta
from boilerplate_shared_module.tools.api_infos import ApiMetadata
from ..handlers.workspaces import *
import uuid
from .users import get_authenticated_user


router = APIRouter()

## SESSIONS ##

@router.get("/workspaces", response_model=ResponseCapsule)
async def route_read_workspaces(
    limit: int | None = None, 
    api_meta: ApiMetadata = Depends(get_api_meta),
    ):
    """
    Retrieves workspaces from the database via HTTP GET request.

    Parameters
    ----------
    limit : int, optional
        The maximum number of workspaces to retrieve, by default None.
    api_meta : ApiMetadata
        Metadata for the API request.

    Returns
    -------
    ResponseCapsule
        The response capsule containing the retrieved workspaces (list[WorkSpaceModel]).

    Raises
    ------
    HTTPException
        If an error occurs during the retrieval process.

    """
    response =  await handler_get_workspaces(
        limit=limit, 
        api_meta=api_meta,
        )
    if response.http_code != 200:
        raise HTTPException(status_code=response.http_code, detail=response.model_dump())
    return response


@router.get("/workspace/{workspace_id}", response_model=ResponseCapsule)
async def route_read_workspace(
    workspace_id: uuid.UUID, 
    api_meta: ApiMetadata = Depends(get_api_meta),
    ):
    """
    Retrieves a specific workspace from the database by ID via HTTP GET request.

    Parameters
    ----------
    workspace_id : UUID
        The ID of the workspace to retrieve.
    api_meta : ApiMetadata
        Metadata for the API request.

    Returns
    -------
    ResponseCapsule
        The response capsule containing the retrieved workspace (WorkSpaceModel).

    Raises
    ------
    HTTPException
        If the workspace with the specified ID is not found or an error occurs during the retrieval process.

    See Also
    --------
    route_read_workspace_by_username, route_create_new_workspace, route_put_workspace, route_soft_delete_workspace, boilerplate_fs_api.shared_models.orm.WorkSpaceModel
    """
    response = await handler_get_workspace_by_id(
        workspace_id=workspace_id, 
        api_meta=api_meta
        )
    if response.http_code != 200: 
        raise HTTPException(status_code=response.http_code, detail=response.model_dump())
    return response


@router.get("/workspace", response_model=ResponseCapsule)
async def route_read_workspace_by_username(
    username: str, 
    api_meta: ApiMetadata = Depends(get_api_meta),
    ):
    """
    Retrieves a specific workspace from the database by username via HTTP GET request.

    Parameters
    ----------
    username : str
        The username associated with the workspace to retrieve.
    api_meta : ApiMetadata
        Metadata for the API request.

    Returns
    -------
    ResponseCapsule
        The response capsule containing the retrieved workspace (WorkSpaceModel).

    Raises
    ------
    HTTPException
        If the workspace with the specified username is not found or an error occurs during the retrieval process.

    See Also
    --------
    route_read_workspace, route_read_workspace_by_id, route_create_new_workspace, route_put_workspace, route_soft_delete_workspace, boilerplate_fs_api.shared_models.orm.WorkSpaceModel
    """
    response = await handler_get_workspace_by_username(
        username=username, 
        api_meta=api_meta
        )
    if response.http_code != 200: 
        raise HTTPException(status_code=response.http_code, detail=response.model_dump())
    return response


@router.post("/workspace", response_model=ResponseCapsule, status_code=201)
async def route_create_new_workspace(
    workspace: WorkSpaceModel, 
    db_user: User = Depends(get_authenticated_user),
    api_meta: ApiMetadata = Depends(get_api_meta),
    ):
    """
    Creates a new workspace in the database via HTTP POST request.

    Parameters
    ----------
    workspace : WorkSpaceModel
        The workspace to be created.
    db_user : User
        The authenticated user.
    api_meta : ApiMetadata
        Metadata for the API request.

    Returns
    -------
    ResponseCapsule
        The response capsule containing the created workspace (WorkSpaceModel).

    Raises
    ------
    HTTPException
        If the creation process fails.

    See Also
    --------
    route_read_workspace, route_read_workspace_by_username, route_put_workspace, route_soft_delete_workspace, boilerplate_fs_api.shared_models.orm.WorkSpaceModel
    """
    response =  await handler_post_workspace(
        workspace=workspace, 
        db_user=db_user,
        api_meta=api_meta
        )
    if response.http_code != 200:
        raise HTTPException(status_code=response.http_code, detail=response.model_dump())
    return response


@router.put("/workspace/{workspace_id}", response_model=ResponseCapsule, status_code=201)
async def route_put_workspace(
    workspace: WorkSpaceModel, 
    workspace_id: uuid.UUID, 
    db_user: User = Depends(get_authenticated_user),
    api_meta: ApiMetadata = Depends(get_api_meta),
    ):
    """
    Updates an existing workspace in the database via HTTP PUT request.

    Parameters
    ----------
    workspace : WorkSpaceModel
        The workspace with updated data.
    workspace_id : UUID
        The ID of the workspace to update.
    db_user : User
        The authenticated user.

    Returns
    -------
    ResponseCapsule
        The response capsule containing the updated workspace (WorkSpaceModel).

    Raises
    ------
    HTTPException
        If the update process fails.

    See Also
    --------
    route_read_workspace, route_read_workspace_by_token, route_create_new_workspace, route_put_workspace, route_soft_delete_workspace, boilerplate_fs_api.shared_models.orm.WorkSpaceModel

    """
    response =  await handler_put_workspace(
        workspace=workspace, 
        workspace_id=workspace_id, 
        db_user=db_user,
        api_meta=api_meta
        )
    if response.http_code != 200:
        raise HTTPException(status_code=response.http_code, detail=response.model_dump())
    return response


@router.delete("/workspace/{workspace_id}", response_model=ResponseCapsule, status_code=201)
async def route_soft_delete_workspace(
    workspace_id: uuid.UUID, 
    db_user: User = Depends(get_authenticated_user),
    api_meta: ApiMetadata = Depends(get_api_meta),
    ):
    """
    Soft deletes an workspace in the database via HTTP DELETE request.

    Parameters
    ----------
    workspace_id : UUID
        The ID of the workspace to delete.
    db_user : User
        The authenticated user.

    Returns
    -------
    ResponseCapsule
        The response capsule containing the soft deleted workspace (WorkSpaceModel).

    Raises
    ------
    HTTPException
        If the deletion process fails.
        

    See Also
    --------
    route_read_workspace, route_read_workspace_by_token, route_create_new_workspace, route_put_workspace, route_soft_delete_workspace, boilerplate_fs_api.shared_models.orm.WorkSpaceModel

    """
    response =  await handler_soft_delete_workspace(
        workspace_id=workspace_id, 
        db_user=db_user,
        api_meta=api_meta
        )
    if response.http_code != 200:
        # raise HTTPException(status_code=response.http_code, detail=response.model_dump())
        raise HTTPException(status_code=response.http_code, detail=response.model_dump())
    return response
