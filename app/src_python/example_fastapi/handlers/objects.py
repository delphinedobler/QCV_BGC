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
# from ..logic.storage.cache import manager as redis_management
from ..logic.storage.database import manager as db_management


# redis_manager: redis_management.RedisManager = redis_management.redis_manager
db_manager: db_management.DatabaseManager = db_management.db_manager


## OBJECTS ##

async def handler_get_objects(
    preloaded: bool, 
    limit: int | None = None, 
    db: AsyncSession = Depends(db_manager.async_get_db), 
    api_meta: ApiMetadata = Depends(get_api_meta),
    ):
    """
    Handles the GET request to retrieve objects.

    Parameters
    ----------
    preloaded : bool
        Whether to retrieve preloaded data.
    limit : int, optional
        The maximum number of objects to retrieve, by default None.
    db : AsyncSession
        The asynchronous database session.

    Returns
    -------
    ResponseCapsule
        The response capsule containing the retrieved objects (list[ApiObject]).

    Raises
    ------
    ResponseCapsule
        If any error occurs during the process.

    See Also
    --------
    get_preloaded_objects, get_objects, boilerplate_fs_api.shared_models.orm.ApiObject
    """
    api_meta_dic = api_meta_to_rc_dic(api_meta=api_meta)
    try:
        if preloaded:
            db_objects = await get_preloaded_objects(db, limit=limit)
            api_model = ApiPreloadedObject
        else:
            db_objects = await get_objects(db, limit=limit)
            api_model = ApiObject
    except Exception as e:
        msg = f'A problem occured while reading objects in database : {e}'
        error_traceback = str(traceback.format_exc())
        return ResponseCapsule(http_code=500, error=error_traceback, message=msg, **api_meta_dic)
    try:
        ta = TypeAdapter(list[api_model])
        db_objects = ta.validate_python(db_objects)
        typee = get_api_model_name(api_model=api_model)
        return ResponseCapsule(data=db_objects, type=typee)
    except ValidationError as e:
        msg = f'A problem occured during the validation of the orm model : {e}'
        error_traceback = str(traceback.format_exc())
        return ResponseCapsule(http_code=422, error=error_traceback, message=msg, **api_meta_dic)


async def handler_get_object_by_id(
    object_id: uuid.UUID, 
    preloaded: bool, 
    db: AsyncSession = Depends(db_manager.async_get_db),
    api_meta: ApiMetadata = Depends(get_api_meta),
    ) -> ResponseCapsule:
    """
    Handles the GET request to retrieve a specific object by ID.

    Parameters
    ----------
    object_id : UUID
        The ID of the object to retrieve.
    preloaded : bool
        Whether to retrieve preloaded data.
    db : AsyncSession
        The asynchronous database session.

    Returns
    -------
    ResponseCapsule
        The response capsule containing the retrieved object (ApiObject).

    Raises
    ------
    ResponseCapsule
        If any error occurs during the process.

    See Also
    --------
    get_preloaded_object, get_object, boilerplate_fs_api.shared_models.orm.ApiObject
    """
    api_meta_dic = api_meta_to_rc_dic(api_meta=api_meta)
    try:
        if preloaded:
            db_obj = await get_preloaded_object(db=db, object_id=object_id)
            api_model = ApiPreloadedObject
        else:
            db_obj = await get_object(db=db, object_id=object_id)
            api_model = ApiObject
    except Exception as e:
        msg = f'A problem occured while reading objects in database : {e}'
        error_traceback = str(traceback.format_exc())
        return ResponseCapsule(http_code=500, error=error_traceback, message=msg, **api_meta_dic)
    if not db_obj:
        msg = "Missing object with object_id given"
        return ResponseCapsule(http_code=400, error=msg, message=msg, **api_meta_dic)
    try:
        db_obj = api_model.model_validate(db_obj)
        typee = get_api_model_name(api_model=api_model)
        return ResponseCapsule(data=db_obj, type=typee)
    except ValidationError as e:
        msg = f'A problem occured during the validation of the orm model : {e}'
        error_traceback = str(traceback.format_exc())
        return ResponseCapsule(http_code=422, error=error_traceback, message=msg, **api_meta_dic)


async def handler_get_object_by_name(
    name: str, 
    preloaded: bool, 
    db: AsyncSession = Depends(db_manager.async_get_db),
    api_meta: ApiMetadata = Depends(get_api_meta),
    ) -> ResponseCapsule:
    """
    Handles the GET request to retrieve a specific object by name.

    Parameters
    ----------
    name : str
        The name of the object to retrieve.
    preloaded : bool
        Whether to retrieve preloaded data.
    db : AsyncSession
        The asynchronous database session.

    Returns
    -------
    ResponseCapsule
        The response capsule containing the retrieved object (ApiObject).

    Raises
    ------
    ResponseCapsule
        If any error occurs during the process.

    See Also
    --------
    get_preloaded_object_by_name, get_object_by_name, boilerplate_fs_api.shared_models.orm.ApiObject
    """
    api_meta_dic = api_meta_to_rc_dic(api_meta=api_meta)
    try:
        if preloaded:
            db_obj = await get_preloaded_object_by_name(db=db, name=name)
            api_model = ApiPreloadedObject
        else:
            db_obj = await get_object_by_name(db=db, name=name)
            api_model = ApiObject
    except Exception as e:
        msg = f'A problem occured while reading objects in database : {e}'
        error_traceback = str(traceback.format_exc())
        return ResponseCapsule(http_code=500, error=error_traceback, message=msg, **api_meta_dic)
    if not db_obj:
        msg = "Missing object with name given"
        return ResponseCapsule(http_code=400, error=msg, message=msg, **api_meta_dic)
    try:
        db_obj = api_model.model_validate(db_obj)
        typee = get_api_model_name(api_model=api_model)
        return ResponseCapsule(data=db_obj, type=typee)
    except ValidationError as e:
        msg = f'A problem occured during the validation of the orm model : {e}'
        error_traceback = str(traceback.format_exc())
        return ResponseCapsule(http_code=422, error=error_traceback, message=msg, **api_meta_dic)



async def handler_post_object(
    object: ApiObject, 
    db_user: User, 
    db: AsyncSession = Depends(db_manager.async_get_db),
    api_meta: ApiMetadata = Depends(get_api_meta),
    ) -> ResponseCapsule:
    """
    Handles the POST request to create a new object.

    Parameters
    ----------
    object : ApiObject
        The object to be created.
    db_user : User
        The user performing the operation.
    db : AsyncSession
        The asynchronous database session.

    Returns
    -------
    ResponseCapsule
        The response capsule containing the created object (ApiObject).

    Raises
    ------
    ResponseCapsule
        If any error occurs during the process.

    See Also
    --------
    get_object_by_name, create_object, boilerplate_fs_api.shared_models.orm.ApiObject
    """
    api_meta_dic = api_meta_to_rc_dic(api_meta=api_meta)
    try:
        db_object = await get_object_by_name(db=db, name=object.name)
        api_model = ApiObject
    except Exception as e:
        msg = f'A problem occured while reading object in database : {e}'
        error_traceback = str(traceback.format_exc())
        return ResponseCapsule(http_code=500, error=error_traceback, message=msg, **api_meta_dic)
    if db_object:
        msg = "Name is already registered"
        return ResponseCapsule(http_code=400, error=msg, message=msg, **api_meta_dic)
    db_object = Object(**dict(object))
    data = await create_object(db=db, db_object=db_object, db_user=db_user)
    try:
        data = api_model.model_validate(data)
        typee = get_api_model_name(api_model=api_model)
        
        return ResponseCapsule(data=data, type=typee)
    except ValidationError as e:
        msg = f'A problem occured during the validation of the orm model : {e}'
        error_traceback = str(traceback.format_exc())
        return ResponseCapsule(http_code=422, error=error_traceback, message=msg, **api_meta_dic)
    
    
async def handler_put_object(
    object: ApiObject, 
    object_id: uuid.UUID, 
    db_user: User, 
    db: AsyncSession = Depends(db_manager.async_get_db),
    api_meta: ApiMetadata = Depends(get_api_meta),
    ) -> ResponseCapsule:
    """
    Handles the PUT request to update an existing object.

    Parameters
    ----------
    object : ApiObject
        The object with updated data.
    object_id : UUID
        The ID of the object to be updated.
    db_user : User
        The user performing the operation.
    db : AsyncSession
        The asynchronous database session.

    Returns
    -------
    ResponseCapsule
        The response capsule containing the updated object (ApiObject).

    Raises
    ------
    ResponseCapsule
        If any error occurs during the process.

    See Also
    --------
    get_object, update_object, boilerplate_fs_api.shared_models.orm.ApiObject
    """
    api_meta_dic = api_meta_to_rc_dic(api_meta=api_meta)
    try:
        db_object = await get_object(db=db, object_id=object_id)
        api_model = ApiObject
    except Exception as e:
        
        msg = f'A problem occured while reading object in database : {e}'
        error_traceback = str(traceback.format_exc())
        return ResponseCapsule(http_code=500, error=error_traceback, message=msg, **api_meta_dic)
    if not db_object:
        msg = "Missing object with object_id given"
        return ResponseCapsule(http_code=400, error=msg, message=msg, **api_meta_dic)
    for key, val in dict(object).items():
        # replace attributes that have values
        if val:
            setattr(db_object, key, val)
    setattr(db_object, 'id', object_id)
    data = await update_object(db=db, db_object=db_object, db_user=db_user)
    try:
        data = api_model.model_validate(data)
        typee = get_api_model_name(api_model=api_model)
        
        return ResponseCapsule(data=data, type=typee)
    except ValidationError as e:
        msg = f'A problem occured during the validation of the orm model : {e}'
        error_traceback = str(traceback.format_exc())
        return ResponseCapsule(http_code=422, error=error_traceback, message=msg, **api_meta_dic)
                
    
async def handler_soft_delete_object(
    object_id: uuid.UUID, 
    db_user: User, 
    db: AsyncSession = Depends(db_manager.async_get_db),
    api_meta: ApiMetadata = Depends(get_api_meta),
    ) -> ResponseCapsule:
    """
    Handles the soft delete request for an object.

    Parameters
    ----------
    object_id : UUID
        The ID of the object to be soft deleted.
    db_user : User
        The user performing the operation.
    db : AsyncSession
        The asynchronous database session.

    Returns
    -------
    ResponseCapsule
        The response capsule containing the deleted object (ApiObject).

    Raises
    ------
    ResponseCapsule
        If any error occurs during the process.

    See Also
    --------
    get_object, soft_delete_object, boilerplate_fs_api.shared_models.orm.ApiObject
    """
    api_meta_dic = api_meta_to_rc_dic(api_meta=api_meta)
    try:
        db_object = await get_object(db=db, object_id=object_id)
        api_model = ApiObject
    except Exception as e:
        
        msg = f'A problem occured while reading object in database : {e}'
        error_traceback = str(traceback.format_exc())
        return ResponseCapsule(http_code=500, error=error_traceback, message=msg, **api_meta_dic)
    if not db_object:
        msg = "Missing object with object_id given"
        return ResponseCapsule(http_code=400, error=msg, message=msg, **api_meta_dic)

    data = await soft_delete_object(db=db, db_object=db_object, db_user=db_user)
    try:
        data = api_model.model_validate(data)
        typee = get_api_model_name(api_model=api_model)
        
        return ResponseCapsule(data=data, type=typee)
    except ValidationError as e:
        msg = f'A problem occured during the validation of the orm model : {e}'
        error_traceback = str(traceback.format_exc())
        return ResponseCapsule(http_code=422, error=error_traceback, message=msg, **api_meta_dic)
