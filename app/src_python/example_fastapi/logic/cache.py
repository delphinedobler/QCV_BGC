import traceback
from pydantic import TypeAdapter
import asyncio
import logging
import uuid

from boilerplate_shared_module.schemes_pydantic.orm import *
from boilerplate_shared_module.tools.models import update_obj_fields

from .orm import *
from .storage.database import manager as db_management
from .storage.cache import manager as redis_management

logger = logging.getLogger('cache-logic')

redis_manager: redis_management.RedisManager = redis_management.redis_manager
db_manager: db_management.DatabaseManager = db_management.db_manager


async def get_users_cache_from_db():
    async with db_manager.AsyncSessionLocal() as db:
        users = await get_users(db)
        ta = TypeAdapter(List[ApiUser])
        users = ta.validate_python(users)
    return {
        user.email: user for user in users
    }
    
    
async def add_or_refresh_user_in_cache(api_user: ApiUser, 
                                       db: AsyncSession | None = None,
                                       ) -> ApiUser:
    # Set user last seen field
    users_cache = await redis_manager.get_users_cache()
    email = api_user.email
    existing_user = users_cache.get(email, None)
    
    date_now = datetime.now(tz=timezone.utc)
    
    # Modify the user in the cache if it already exists
    if existing_user:
        new_user = update_obj_fields(old_obj=existing_user, new_obj=api_user)
        users_cache[email] = new_user
    
    # If user doesn't exists in cache, so it must be created (using databse if one exists, else manually)
    else:
        if db is not None:
            # If database exists => use database to create user
            db_user = User(**api_user.model_dump())
            # Create the user in the database
            orm_user = await create_user(db=db, db_user=db_user)
            # Validate the orm user into an api user
            api_user = ApiUser.model_validate(orm_user)
                
        else:
            # Instanciate some fields of the user manually (directly done by database when it exists)
            api_user.created_at = date_now
            api_user.updated_at = date_now
            api_user.id = uuid.uuid4()
            
        # Add the user in the cache dic    
        users_cache[email] = api_user
            
    await redis_manager.set_users_cache(value=users_cache)
    return users_cache[email]
        
        
async def sync_cache_user_with_db(db: AsyncSession):
    users_cache = await redis_manager.get_users_cache()
    new_users_cache = {}
    for email, api_user in users_cache.items():
        # Get the user from the database and update the object
        existing_db_user = await get_user(db=db, 
                                          user_id=api_user.id)
        db_user = update_obj_fields(old_obj=existing_db_user,
                                     new_obj=api_user)
        db_user = await update_user(db=db, db_user=db_user)
        
        api_user = ApiUser.model_validate(db_user)
        
        # Keep the session only if it hasn't been deleted
        if api_user.deleted_at is None:
            new_users_cache[email] = api_user
    # Update the cache users
    await redis_manager.set_users_cache(value=new_users_cache)


async def create_session_in_db(db: AsyncSession, 
                               api_session: ApiSessionLog,
                               api_user: ApiUser) -> ApiSessionLog:
    # Translate the api user to an orm user
    db_session = SessionLog(**dict(api_session))
    # Create the user in the database
    orm_session = await create_session_log(db=db, db_session_log=db_session, db_user=api_user)
    # Validate the orm user into an api user
    new_session = ApiSessionLog.model_validate(orm_session)
    return new_session
    
async def add_or_refresh_session_in_cache(api_session: ApiSessionLog, 
                                          api_user: ApiUser, 
                                          db: AsyncSession | None = None, 
                                          ) -> ApiSessionLog:
    # Set session last seen field
    last_seen = datetime.now(tz=timezone.utc)
    api_session.last_seen = last_seen
    
    sessions_cache = await redis_manager.get_sessions_cache()
    existing_sessions = sessions_cache.get(api_user.email, [])
    session_id = api_session.id    
    # If a session id exists, it means that the session already exists in the cache
    if session_id:
        searched_session = await redis_manager.get_session_cache(
            user_email=api_user.email,
            session_id=session_id,
        )
        # searched_session = api_session.model_copy()
        # searched_session.id = session_id
        index = existing_sessions.index(searched_session.session)
        existing_session = existing_sessions[index]
        new_session = update_obj_fields(
            old_obj=existing_session,
            new_obj=api_session
            )
        existing_sessions[index] = new_session
        
    # Session must be created (cache and, if possible database)
    else:
        if db is not None:
            # If database exists => use database to create user
            new_session = await create_session_in_db(
                db=db,
                api_session=api_session,
                api_user=api_user
                )
        else:
            # Instanciate some fields of the user manually (directly done by database when it exists)
            new_session = api_session
            new_session.id = uuid.uuid4()
            # Use fields of the user (created before) to instanciate some fields of the session
            new_session.created_at = api_user.created_at
            new_session.updated_at = api_user.updated_at
            new_session.user_id = api_user.id
    
        existing_sessions.append(new_session)
            
    await redis_manager.set_sessions_cache_for_user(
        user_email=api_user.email, 
        value=existing_sessions
        )
    return new_session

        
async def sync_cache_session_with_db(db: AsyncSession):
    # global SESSIONS_CACHE
    sessions_cache = await redis_manager.get_sessions_cache()
    for email, api_sessions in sessions_cache.items():
        new_api_sessions = []
        for api_session in api_sessions:
            # Get the session from the database and update the object
            existing_db_session = await get_session_log(db=db, 
                                               session_log_id=api_session.id)
            existing_db_user = await get_user_by_email(db=db, email=email)
            db_session = update_obj_fields(
                old_obj=existing_db_session,
                new_obj=api_session
                )
            db_session = await update_session_log(db=db, db_session_log=db_session, db_user=existing_db_user)
            api_session = ApiSessionLog.model_validate(db_session)
            # Keep the session only if it hasn't been deleted
            if api_session.deleted_at is None:
                new_api_sessions.append(api_session)
        # Update the cache sessions
        await redis_manager.set_sessions_cache_for_user(
            user_email=email,
            value=new_api_sessions
            )
        # sessions_cache[email] = new_api_sessions
    # await set_sessions_cache(sessions_cache)


async def initialize_users_cache(
    from_db: bool = False
    ):
        if not from_db:
            await redis_manager.set_users_cache(value={})
        else:
            async with db_manager.AsyncSessionLocal() as db:
                orm_users = await get_users(db=db, limit=None)
            ta = TypeAdapter(list[ApiUser])
            users: list[ApiUser] = ta.validate_python(orm_users)
            users_dic = {user.email: user for user in users}
            await redis_manager.set_users_cache(value=users_dic)
            
            
async def initialize_sessions_cache(
    from_db: bool = False
    ):
        if not from_db:
            await redis_manager.set_sessions_cache(value={})
        else:
            users_dic = await redis_manager.get_users_cache()
            async with db_manager.AsyncSessionLocal() as db:
                orm_sessions = await get_session_logs(db=db, limit=None)
            ta = TypeAdapter(list[ApiSessionLog])
            sessions: list[ApiSessionLog] = ta.validate_python(orm_sessions)
            sessions_dic = {}
            for user in users_dic.values():
                sessions_dic[user.email] = []
                for session in sessions:
                    if session.user_id == user.id:
                        sessions_dic[user.email] += [session]
            await redis_manager.set_sessions_cache(value=sessions_dic)
            

async def initialize_cache(from_db: bool = False):
    try:
    #     if not os.path.exists(redis_manager.init_lock_path):
    #         with open(redis_manager.init_lock_path, 'w') as lock_file:
        await initialize_users_cache(from_db=from_db)
        await initialize_sessions_cache(from_db=from_db)
    except Exception as e:
        print(e)
        print(traceback.format_exc())
    #             lock_file.write("Cache initialized")
    #         print('Cache initialized')
    # except Exception as e:
    #     print(e)


async def reset_cache():
    await redis_manager.reset_users_cache()
    await redis_manager.reset_sessions_cache()


async def add_or_refresh_user_infos_in_cache_and_db(user_infos: ApiUserInfos, 
                                             db: AsyncSession | None = None, 
                                             ) -> ApiUserInfos:
    api_user = user_infos.user
    session = user_infos.session
    new_api_user = await add_or_refresh_user_in_cache(
        api_user=api_user, 
        db=db
        )
    new_api_session = await add_or_refresh_session_in_cache(
        api_session=session, 
        api_user=new_api_user,
        db=db,
        )
    return ApiUserInfos(user=new_api_user, session=new_api_session)


async def sync_cache_with_db():
    logger.info('Syncing database with cache...')
    try:
        async with db_manager.AsyncSessionLocal() as db:
            await sync_cache_user_with_db(db=db)
            await sync_cache_session_with_db(db=db)
        logger.info('Database has been successfully synced with cache...')
    except Exception as e:
        logger.error(e)
        logger.error(str(traceback.format_exc()))


async def sync_cache_with_db_periodically(
    time_s=30, 
    # stop_event: asyncio.Event = None
    ):
    try:
        # while not stop_event.is_set():
            await sync_cache_with_db()
            await asyncio.sleep(time_s)
    except Exception as e:
        logger.error(e)
        logger.error(str(traceback.format_exc()))