from pymongo import MongoClient, AsyncMongoClient, ReturnDocument
from bson.objectid import ObjectId
from datetime import datetime, timezone
from .storage.mongodb.manager import mongo_manager
from uuid import uuid4, UUID


"""
Here are stored functions used by the handlers to apply actions on the database. These functions don't manage errors, these are
managed at the handler level.
"""

db_name = mongo_manager.db_name
collection_name = 'workspaces'


async def get_workspaces(
    mongo_client: AsyncMongoClient, 
    limit: int | None = None,
    ) -> dict | None:
    db = mongo_manager.get_db(mongo_client=mongo_client)
    collection = db[collection_name]
    workspaces = []
    cursor = collection.find(
        {"deleted_at": None}
        )
    if limit is not None:
        cursor = cursor.limit(limit)
    async for workspace in cursor:
        workspaces.append(workspace)
    return workspaces


async def get_workspace(
    mongo_client: AsyncMongoClient, 
    workspace_id: str | UUID,
    ) -> dict | None:
    db = mongo_manager.get_db(mongo_client=mongo_client)
    collection = db[collection_name]
    if isinstance(workspace_id, str):
        workspace_id = UUID(workspace_id)
    workspace = await collection.find_one(
        {
            "id": workspace_id,
            "deleted_at": None,
            }
        )
    return workspace


async def get_workspace_by_username(
    mongo_client: AsyncMongoClient, 
    username: str,
    ) -> dict | None:
    db = mongo_manager.get_db(mongo_client=mongo_client)
    collection = db[collection_name]
    workspace = await collection.find_one(
        {
            "username": username,
            "deleted_at": None,
            }
        )
    return workspace


async def create_workspace(
    mongo_client: AsyncMongoClient, 
    workspace_dic: dict,
    user_id: str | UUID | None = None,
    ) -> dict:
    db = mongo_manager.get_db(mongo_client=mongo_client)
    collection = db[collection_name]
    date_now = datetime.now(tz=timezone.utc)
    
    ws_id = workspace_dic.get('id')
    if ws_id is None:
        ws_id = uuid4()
    if isinstance(ws_id, str):
        ws_id = UUID(ws_id) 
    if isinstance(user_id, str):
        user_id = UUID(user_id)    
    
    workspace_dic |= {
        "id": ws_id,
        "_id": ws_id,
        "created_by": user_id,
        "created_at": date_now,
        "updated_by": user_id,
        "updated_at": date_now,
        "deleted_at": None,
        "active": True,
    }
    # username = workspace_json_dict.get('username')
    # try:
    #     existing_workspace = await collection.find_one({"username": username})
    #     if existing_workspace:
    #         raise ValueError(f"Workspace with {username=} already exists.")
    #     await collection.insert_one(workspace_json_dict)
    # except Exception as e:
    #     raise e
    await collection.insert_one(workspace_dic)
    return workspace_dic


async def update_workspace(
    mongo_client: AsyncMongoClient, 
    workspace_id: str | UUID, 
    workspace_dic: dict,
    user_id: str | UUID | None = None,
    ) -> dict | None:
    db = mongo_manager.get_db(mongo_client=mongo_client)
    collection = db[collection_name]
    
    if isinstance(workspace_id, str):
        workspace_id = UUID(workspace_id)
    if isinstance(user_id, str):
        user_id = UUID(user_id)
    
    date_now = datetime.now(tz=timezone.utc)
    workspace_dic['id'] = workspace_id
    workspace_dic |= {
        "updated_by": user_id,
        "updated_at": date_now,
        "active": True,
    }
    
    result = await collection.find_one_and_update(
        {
            'id': workspace_id,
            "deleted_at": None,
            }, 
        {"$set": workspace_dic},
        return_document=ReturnDocument.AFTER
        )
    # if result.modified_count == 1:
    #     updated_workspace = await collection.find_one({"id": id})
    #     return updated_workspace
    return result


async def soft_delete_workspace(
    mongo_client: AsyncMongoClient, 
    workspace_id: str,
    user_id: str | UUID | None = None,
    ) -> dict | None:
    db = mongo_manager.get_db(mongo_client=mongo_client)
    collection = db[collection_name]
    
    if isinstance(workspace_id, str):
        workspace_id = UUID(workspace_id)
    if isinstance(user_id, str):
        user_id = UUID(user_id)
    
    date_now = datetime.now(tz=timezone.utc)
    result = await collection.find_one_and_update(
        {
            'id': workspace_id,
            "deleted_at": None,
            }, 
        {
            "$set": {
                "updated_at": date_now,
                "updated_by": user_id,
                "deleted_at": date_now,
                "deleted_by": user_id,
                "active": False,
            }
        },
        return_document=ReturnDocument.AFTER
        )
    # if result.modified_count == 1:
    #     updated_workspace = await collection.find_one({"id": id})
    #     return updated_workspace
    return result



async def delete_workspace(
    mongo_client: AsyncMongoClient, 
    workspace_id: str | UUID,
    ) -> bool:
    db = mongo_manager.get_db(mongo_client=mongo_client)
    collection = db[collection_name]
    
    if isinstance(workspace_id, str):
        workspace_id = UUID(workspace_id)
    
    result = await collection.delete_one({"id": workspace_id})
    return result
