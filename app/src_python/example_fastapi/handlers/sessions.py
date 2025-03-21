from sqlalchemy.ext.asyncio import AsyncSession
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
from ..logic.storage.cache import manager as redis_management
from ..logic.storage.database import manager as db_management


redis_manager: redis_management.RedisManager = redis_management.redis_manager
db_manager: db_management.DatabaseManager = db_management.db_manager


## SESSIONS ##

async def handler_get_session_logs(
    preloaded: bool, 
    limit: int | None = None, 
    db: AsyncSession = Depends(db_manager.async_get_db), 
    api_meta: ApiMetadata = Depends(get_api_meta),
    ):
    """
    Handles the GET request to retrieve session_logs.

    Parameters
    ----------
    preloaded : bool
        Whether to retrieve preloaded data.
    limit : int, optional
        The maximum number of session_logs to retrieve, by default None.
    db : AsyncSession
        The asynchronous database session.

    Returns
    -------
    ResponseCapsule
        The response capsule containing the retrieved session_logs (list[ApiSessionLog]).

    Raises
    ------
    ResponseCapsule
        If any error occurs during the process.

    See Also
    --------
    get_preloaded_session_logs, get_session_logs, boilerplate_fs_api.shared_models.orm.ApiSessionLog
    """
    api_meta_dic = api_meta_to_rc_dic(api_meta=api_meta)
    try:
        if preloaded:
            db_session_logs = await get_preloaded_session_logs(db, limit=limit)
            api_model = ApiPreloadedSessionLog
        else:
            db_session_logs = await get_session_logs(db, limit=limit)
            api_model = ApiSessionLog
    except Exception as e:
        msg = f'A problem occured while reading session_logs in database : {e}'
        error_traceback = str(traceback.format_exc())
        return ResponseCapsule(http_code=500, error=error_traceback, **api_meta_dic)
    try:
        ta = TypeAdapter(list[api_model])
        db_session_logs = ta.validate_python(db_session_logs)
        typee = get_api_model_name(api_model=api_model)
        return ResponseCapsule(data=db_session_logs, type=typee)
    except ValidationError as e:
        msg = f'A problem occured during the validation of the orm model : {e}'
        error_traceback = str(traceback.format_exc())
        return ResponseCapsule(http_code=422, error=error_traceback, **api_meta_dic)


async def handler_get_session_log_by_id(
    session_log_id: uuid.UUID, 
    preloaded: bool, 
    db: AsyncSession = Depends(db_manager.async_get_db),
    api_meta: ApiMetadata = Depends(get_api_meta),
    ) -> ResponseCapsule:
    """
    Handles the GET request to retrieve a specific session_log by ID.

    Parameters
    ----------
    session_log_id : UUID
        The ID of the session_log to retrieve.
    preloaded : bool
        Whether to retrieve preloaded data.
    db : AsyncSession
        The asynchronous database session.

    Returns
    -------
    ResponseCapsule
        The response capsule containing the retrieved session_log (ApiSessionLog).

    Raises
    ------
    ResponseCapsule
        If any error occurs during the process.

    See Also
    --------
    get_preloaded_session_log, get_session_log, boilerplate_fs_api.shared_models.orm.ApiSessionLog
    """
    api_meta_dic = api_meta_to_rc_dic(api_meta=api_meta)
    try:
        if preloaded:
            db_obj = await get_preloaded_session_log(db=db, session_log_id=session_log_id)
            api_model = ApiPreloadedSessionLog
        else:
            db_obj = await get_session_log(db=db, session_log_id=session_log_id)
            api_model = ApiSessionLog
    except Exception as e:
        msg = f'A problem occured while reading session_logs in database : {e}'
        error_traceback = str(traceback.format_exc())
        return ResponseCapsule(http_code=500, error=error_traceback, **api_meta_dic)
    if not db_obj:
        msg = "Missing session_log with session_log_id given"
        return ResponseCapsule(http_code=400, error=msg, **api_meta_dic)
    try:
        db_obj = api_model.model_validate(db_obj)
        typee = get_api_model_name(api_model=api_model)
        return ResponseCapsule(data=db_obj, type=typee)
    except ValidationError as e:
        msg = f'A problem occured during the validation of the orm model : {e}'
        error_traceback = str(traceback.format_exc())
        return ResponseCapsule(http_code=422, error=error_traceback, **api_meta_dic)


async def handler_get_session_log_by_token(
    token: str, 
    preloaded: bool, 
    db: AsyncSession = Depends(db_manager.async_get_db),
    api_meta: ApiMetadata = Depends(get_api_meta),
    ) -> ResponseCapsule:
    """
    Handles the GET request to retrieve a specific session_log by token.

    Parameters
    ----------
    token : str
        The token of the session_log to retrieve.
    preloaded : bool
        Whether to retrieve preloaded data.
    db : AsyncSession
        The asynchronous database session.

    Returns
    -------
    ResponseCapsule
        The response capsule containing the retrieved session_log (ApiSessionLog).

    Raises
    ------
    ResponseCapsule
        If any error occurs during the process.

    See Also
    --------
    get_preloaded_session_log_by_token, get_session_log_by_token, boilerplate_fs_api.shared_models.orm.ApiSessionLog
    """
    api_meta_dic = api_meta_to_rc_dic(api_meta=api_meta)
    try:
        if preloaded:
            db_obj = await get_preloaded_session_log_by_token(db=db, token=token)
            api_model = ApiPreloadedSessionLog
        else:
            db_obj = await get_session_log_by_token(db=db, token=token)
            api_model = ApiSessionLog
    except Exception as e:
        msg = f'A problem occured while reading session_logs in database : {e}'
        error_traceback = str(traceback.format_exc())
        return ResponseCapsule(http_code=500, error=error_traceback, **api_meta_dic)
    if not db_obj:
        msg = "Missing session_log with token given"
        return ResponseCapsule(http_code=400, error=msg, **api_meta_dic)
    try:
        db_obj = api_model.model_validate(db_obj)
        typee = get_api_model_name(api_model=api_model)
        return ResponseCapsule(data=db_obj, type=typee)
    except ValidationError as e:
        msg = f'A problem occured during the validation of the orm model : {e}'
        error_traceback = str(traceback.format_exc())
        return ResponseCapsule(http_code=422, error=error_traceback, **api_meta_dic)



async def handler_post_session_log(
    session_log: ApiSessionLog, 
    db_user: User, 
    db: AsyncSession = Depends(db_manager.async_get_db),
    api_meta: ApiMetadata = Depends(get_api_meta),
    ) -> ResponseCapsule:
    """
    Handles the POST request to create a new session_log.

    Parameters
    ----------
    session_log : ApiSessionLog
        The session_log to be created.
    db_user : User
        The user performing the operation.
    db : AsyncSession
        The asynchronous database session.

    Returns
    -------
    ResponseCapsule
        The response capsule containing the created session_log (ApiSessionLog).

    Raises
    ------
    ResponseCapsule
        If any error occurs during the process.

    See Also
    --------
    get_session_log_by_token, create_session_log, boilerplate_fs_api.shared_models.orm.ApiSessionLog
    """
    api_meta_dic = api_meta_to_rc_dic(api_meta=api_meta)
    try:
        db_session_log = await get_session_log_by_token(db=db, token=session_log.token)
        api_model = ApiSessionLog
    except Exception as e:
        msg = f'A problem occured while reading session_log in database : {e}'
        error_traceback = str(traceback.format_exc())
        return ResponseCapsule(http_code=500, error=error_traceback, **api_meta_dic)
    if db_session_log:
        msg = "Shiny id is already registered"
        return ResponseCapsule(http_code=400, error=msg, **api_meta_dic)
    db_session_log = SessionLog(**dict(session_log))
    data = await create_session_log(db=db, db_session_log=db_session_log, db_user=db_user)
    try:
        data = api_model.model_validate(data)
        typee = get_api_model_name(api_model=api_model)
        
        return ResponseCapsule(data=data, type=typee)
    except ValidationError as e:
        msg = f'A problem occured during the validation of the orm model : {e}'
        error_traceback = str(traceback.format_exc())
        return ResponseCapsule(http_code=422, error=error_traceback, **api_meta_dic)
    
    
async def handler_put_session_log(
    session_log: ApiSessionLog, 
    session_log_id: uuid.UUID, 
    db_user: User, 
    db: AsyncSession = Depends(db_manager.async_get_db),
    api_meta: ApiMetadata = Depends(get_api_meta),
    ) -> ResponseCapsule:
    """
    Handles the PUT request to update an existing session_log.

    Parameters
    ----------
    session_log : ApiSessionLog
        The session_log with updated data.
    session_log_id : UUID
        The ID of the session_log to be updated.
    db_user : User
        The user performing the operation.
    db : AsyncSession
        The asynchronous database session.

    Returns
    -------
    ResponseCapsule
        The response capsule containing the updated session_log (ApiSession).

    Raises
    ------
    ResponseCapsule
        If any error occurs during the process.

    See Also
    --------
    get_session_log, update_session_log, boilerplate_fs_api.shared_models.orm.ApiSession
    """
    api_meta_dic = api_meta_to_rc_dic(api_meta=api_meta)
    try:
        db_session_log = await get_session_log(db=db, session_log_id=session_log_id)
        api_model = ApiSessionLog
    except Exception as e:
        
        msg = f'A problem occured while reading session_log in database : {e}'
        error_traceback = str(traceback.format_exc())
        return ResponseCapsule(http_code=500, error=error_traceback, **api_meta_dic)
    if not db_session_log:
        msg = "Missing session_log with session_log_id given"
        return ResponseCapsule(http_code=400, error=msg, **api_meta_dic)
    for key, val in dict(session_log).items():
        # replace attributes that have values
        if val:
            setattr(db_session_log, key, val)
    setattr(db_session_log, 'id', session_log_id)
    data = await update_session_log(db=db, db_session_log=db_session_log, db_user=db_user)
    try:
        data = api_model.model_validate(data)
        typee = get_api_model_name(api_model=api_model)
        
        return ResponseCapsule(data=data, type=typee)
    except ValidationError as e:
        msg = f'A problem occured during the validation of the orm model : {e}'
        error_traceback = str(traceback.format_exc())
        return ResponseCapsule(http_code=422, error=error_traceback, **api_meta_dic)
                
    
async def handler_soft_delete_session_log(
    session_log_id: uuid.UUID, 
    db_user: User, 
    db: AsyncSession = Depends(db_manager.async_get_db),
    api_meta: ApiMetadata = Depends(get_api_meta),
    ) -> ResponseCapsule:
    """
    Handles the soft delete request for an session_log.

    Parameters
    ----------
    session_log_id : UUID
        The ID of the session_log to be soft deleted.
    db_user : User
        The user performing the operation.
    db : AsyncSession
        The asynchronous database session.

    Returns
    -------
    ResponseCapsule
        The response capsule containing the deleted session_log (ApiSessionLog).

    Raises
    ------
    ResponseCapsule
        If any error occurs during the process.

    See Also
    --------
    get_session_log, soft_delete_session_log, boilerplate_fs_api.shared_models.orm.ApiSessionLog
    """
    api_meta_dic = api_meta_to_rc_dic(api_meta=api_meta)
    try:
        db_session_log = await get_session_log(db=db, session_log_id=session_log_id)
        api_model = ApiSessionLog
    except Exception as e:
        
        msg = f'A problem occured while reading session_log in database : {e}'
        error_traceback = str(traceback.format_exc())
        return ResponseCapsule(http_code=500, error=error_traceback, **api_meta_dic)
    if not db_session_log:
        msg = "Missing session_log with session_log_id given"
        return ResponseCapsule(http_code=400, error=msg, **api_meta_dic)

    data = await soft_delete_session_log(db=db, db_session_log=db_session_log, db_user=db_user)
    try:
        data = api_model.model_validate(data)
        typee = get_api_model_name(api_model=api_model)
        
        return ResponseCapsule(data=data, type=typee)
    except ValidationError as e:
        msg = f'A problem occured during the validation of the orm model : {e}'
        error_traceback = str(traceback.format_exc())
        return ResponseCapsule(http_code=422, error=error_traceback, **api_meta_dic)
