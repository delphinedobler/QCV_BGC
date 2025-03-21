from fastapi import HTTPException, Path, Query, Depends, Request
from ..logic.orm import *
from boilerplate_shared_module.schemes_pydantic.orm import *
from boilerplate_shared_module.schemes_pydantic.basic_schemes import *
from boilerplate_shared_module.tools.models import *
from ..metadata.api import get_api_meta
from boilerplate_shared_module.tools.api_infos import ApiMetadata, api_meta_to_rc_dic
from pydantic import ValidationError
from pydantic import TypeAdapter
import traceback
from ..logic.workspaces import *
from boilerplate_shared_module.schemes_pydantic.frontend import WorkSpaceModel



async def handler_get_workspaces(
    limit: int | None = None,
    mongo_client: AsyncMongoClient = mongo_manager.async_client,
    api_meta: ApiMetadata = Depends(get_api_meta),
) -> ResponseCapsule:
    api_meta_dic = api_meta_to_rc_dic(api_meta=api_meta)
    try:
        json_workspaces = await get_workspaces(
            mongo_client=mongo_client, 
            limit=limit,
            )
    except Exception as e:
        msg = f'A problem occurred while reading workspaces in database: {e}'
        error_traceback = str(traceback.format_exc())
        return ResponseCapsule(
            message=msg,
            http_code=500, 
            error=error_traceback, 
            **api_meta_dic,
            )
    try:
        ta = TypeAdapter(list[WorkSpaceModel])
        json_workspaces = ta.validate_python(json_workspaces)
        typee = get_api_model_name(api_model=WorkSpaceModel)
        return ResponseCapsule(data=json_workspaces, type=typee)
    except ValidationError as e:
        msg = f'A problem occurred during the validation of the logic model: {e}'
        error_traceback = str(traceback.format_exc())
        return ResponseCapsule(
            message=msg,
            http_code=422, 
            error=error_traceback, 
            **api_meta_dic,
            )


async def handler_get_workspace_by_id(
    workspace_id: uuid.UUID,
    mongo_client: AsyncMongoClient = mongo_manager.async_client,
    api_meta: ApiMetadata = Depends(get_api_meta),
) -> ResponseCapsule:
    api_meta_dic = api_meta_to_rc_dic(api_meta=api_meta)
    try:
        json_obj = await get_workspace(
            mongo_client=mongo_client, 
            workspace_id=workspace_id,
            )
    except Exception as e:
        msg = f'A problem occurred while reading workspace in database: {e}'
        error_traceback = str(traceback.format_exc())
        return ResponseCapsule(
            message=msg,
            http_code=500, 
            error=error_traceback, 
            **api_meta_dic,
            )
    if not json_obj:
        msg = "Missing workspace with workspace_id given"
        return ResponseCapsule(
            message=msg,
            data=None,
            http_code=200, 
            error=msg, 
            **api_meta_dic,
            )
    try:
        json_obj = WorkSpaceModel.model_validate(json_obj)
        typee = get_api_model_name(api_model=WorkSpaceModel)
        return ResponseCapsule(data=json_obj, type=typee)
    except ValidationError as e:
        msg = f'A problem occurred during the validation of the logic model: {e}'
        error_traceback = str(traceback.format_exc())
        return ResponseCapsule(
            message=msg,
            http_code=422, 
            error=error_traceback, 
            **api_meta_dic,
            )
        

async def handler_get_workspace_by_username(
    username: str,
    mongo_client: AsyncMongoClient = mongo_manager.async_client,
    api_meta: ApiMetadata = Depends(get_api_meta),
) -> ResponseCapsule:
    api_meta_dic = api_meta_to_rc_dic(api_meta=api_meta)
    try:
        json_obj = await get_workspace_by_username(
            mongo_client=mongo_client, 
            username=username,
            )
    except Exception as e:
        msg = f'A problem occurred while reading workspace in database: {e}'
        error_traceback = str(traceback.format_exc())
        return ResponseCapsule(
            message=msg,
            http_code=500, 
            error=error_traceback, 
            **api_meta_dic,
            )
    if not json_obj:
        msg = "Missing workspace with username given"
        return ResponseCapsule(
            message=msg,
            data=None,
            http_code=200, 
            error=msg, 
            **api_meta_dic,
            )
    try:
        json_obj = WorkSpaceModel.model_validate(json_obj)
        typee = get_api_model_name(api_model=WorkSpaceModel)
        return ResponseCapsule(data=json_obj, type=typee)
    except ValidationError as e:
        msg = f'A problem occurred during the validation of the logic model: {e}'
        error_traceback = str(traceback.format_exc())
        return ResponseCapsule(
            message=msg,
            http_code=422, 
            error=error_traceback, 
            **api_meta_dic,
            )


async def handler_post_workspace(
    workspace: WorkSpaceModel,
    db_user: User,
    mongo_client: AsyncMongoClient = mongo_manager.async_client,
    api_meta: ApiMetadata = Depends(get_api_meta),
) -> ResponseCapsule:
    api_meta_dic = api_meta_to_rc_dic(api_meta=api_meta)
    try:
        existing_json_workspace = await get_workspace(
            mongo_client=mongo_client, 
            workspace_id=workspace.id,
            )
        
    except Exception as e:
        msg = f'A problem occurred while reading workspace in database: {e}'
        error_traceback = str(traceback.format_exc())
        return ResponseCapsule(
            message=msg,
            http_code=500, 
            error=error_traceback, 
            **api_meta_dic,
            )
    if existing_json_workspace:
        msg = "Workspace id is already registered"
        return ResponseCapsule(
            message=msg,
            http_code=400, 
            error=msg, 
            **api_meta_dic,
            )
    json_workspace = await create_workspace(
        mongo_client=mongo_client,
        workspace_dic=workspace.model_dump(mode='json'), 
        user_id=db_user.id,
        )
    try:
        ws = WorkSpaceModel.model_validate(json_workspace)
        typee = get_api_model_name(api_model=WorkSpaceModel)
        return ResponseCapsule(data=ws, type=typee)
    except ValidationError as e:
        msg = f'A problem occurred during the validation of the logic model: {e}'
        error_traceback = str(traceback.format_exc())
        return ResponseCapsule(
            message=msg,
            http_code=422, 
            error=error_traceback, 
            **api_meta_dic,
            )


async def handler_put_workspace(
    workspace: WorkSpaceModel,
    workspace_id: uuid.UUID,
    db_user: User,
    mongo_client: AsyncMongoClient = mongo_manager.async_client,
    api_meta: ApiMetadata = Depends(get_api_meta),
) -> ResponseCapsule:
    api_meta_dic = api_meta_to_rc_dic(api_meta=api_meta)
    try:
        existing_json_workspace = await get_workspace(
            mongo_client=mongo_client,
            workspace_id=workspace_id,
            )
    except Exception as e:
        msg = f'A problem occurred while reading workspace in database: {e}'
        error_traceback = str(traceback.format_exc())
        return ResponseCapsule(
            message=msg,
            http_code=500, 
            error=error_traceback, 
            **api_meta_dic,
            )
    if not existing_json_workspace:
        msg = "Missing workspace with workspace_id given"
        return ResponseCapsule(
            message=msg,
            http_code=400, 
            error=msg, 
            **api_meta_dic,
            )
    # for key, val in dict(workspace).items():
    #     if val:
    #         setattr(json_workspace, key, val)
    # setattr(json_workspace, 'id', workspace_id)
    try:
        json_workspace = workspace.model_dump(mode='json')
        updated_json_ws = await update_workspace(
            mongo_client=mongo_client, 
            workspace_id=workspace_id,
            workspace_dic=json_workspace, 
            user_id=db_user.id,
            )
    except Exception as e:
        msg = f'A problem occurred while upating workspace in database: {e}'
        error_traceback = str(traceback.format_exc())
        return ResponseCapsule(
            message=msg,
            http_code=500, 
            error=error_traceback, 
            **api_meta_dic,
            )
    
    try:
        ws = WorkSpaceModel.model_validate(updated_json_ws)
        typee = get_api_model_name(api_model=WorkSpaceModel)
        return ResponseCapsule(data=ws, type=typee)
    except ValidationError as e:
        msg = f'A problem occurred during the validation of the logic model: {e}'
        error_traceback = str(traceback.format_exc())
        return ResponseCapsule(
            message=msg,
            http_code=422, 
            error=error_traceback, 
            **api_meta_dic,
            )


async def handler_soft_delete_workspace(
    workspace_id: uuid.UUID,
    db_user: User,
    mongo_client: AsyncMongoClient = mongo_manager.async_client,
    api_meta: ApiMetadata = Depends(get_api_meta),
    ) -> ResponseCapsule:
    api_meta_dic = api_meta_to_rc_dic(api_meta=api_meta)
    try:
        json_workspace = await get_workspace(
            mongo_client=mongo_client, 
            workspace_id=workspace_id,
            )
    except Exception as e:
        msg = f'A problem occurred while reading workspace in database: {e}'
        error_traceback = str(traceback.format_exc())
        return ResponseCapsule(
            message=msg,
            http_code=500, 
            error=error_traceback, 
            **api_meta_dic,
            )
    if not json_workspace:
        msg = "Missing workspace with workspace_id given"
        return ResponseCapsule(
            message=msg,
            http_code=400, 
            error=msg, 
            **api_meta_dic,
            )
    try:
        deleted_json_ws = await soft_delete_workspace(
            mongo_client=mongo_client, 
            workspace_id=workspace_id, 
            user_id=db_user.id,
            )
    except Exception as e:
        msg = f'A problem occurred while soft deleting workspace in database: {e}'
        error_traceback = str(traceback.format_exc())
        return ResponseCapsule(
            message=msg,
            http_code=500, 
            error=error_traceback, 
            **api_meta_dic,
            )
    try:
        ws = WorkSpaceModel.model_validate(deleted_json_ws)
        typee = get_api_model_name(api_model=WorkSpaceModel)
        return ResponseCapsule(data=ws, type=typee)
    except ValidationError as e:
        msg = f'A problem occurred during the validation of the logic model: {e}'
        error_traceback = str(traceback.format_exc())
        return ResponseCapsule(
            message=msg,
            http_code=422, 
            error=error_traceback, 
            **api_meta_dic,
            )
