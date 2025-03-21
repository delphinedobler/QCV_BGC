from sqlalchemy.orm import joinedload, selectinload
from sqlalchemy.ext.asyncio import AsyncSession
from sqlalchemy.future import select
from sqlalchemy import or_, and_, func, desc
from datetime import datetime, timezone

from .storage.database.orm_classes import *


"""
Here are stored functions used by the handlers to apply actions on the database. These functions don't manage errors, these are
managed at the handler level.
"""


## USERS ##

async def get_users(db: AsyncSession, limit: int | None = None):
    """
    Get a list of users

    Parameters
    ----------
    db : AsyncSession
        Database async session
    limit : int | None, optional
        Maximum number of users to return, by default None

    Returns
    -------
    List[User]
        List of User objects
    """
    statement = select(User) \
                .where(User.deleted_at == None)
    if isinstance(limit, int):
        statement = statement.limit(limit)
    results = await db.execute(statement)
    results = results.unique().scalars().all()
    return results


async def get_user(db: AsyncSession, user_id: uuid.UUID):
    """
    Get a user by ID

    Parameters
    ----------
    db : AsyncSession
        Database async session
    user_id : uuid.UUID
        ID of the user

    Returns
    -------
    User | None
        User object if found, None otherwise
    """
    result = await db.execute(select(User)
                              .where(and_(User.id == user_id, User.deleted_at == None))
                              )
    result = result.unique().scalar_one_or_none()
    return result


async def get_preloaded_users(db: AsyncSession, limit: int | None = None):
    """
    Get a list of preloaded users with associated floats, f_rules, and p_rules

    Parameters
    ----------
    db : AsyncSession
        Database async session
    limit : int | None, optional
        Maximum number of users to return, by default None

    Returns
    -------
    List[User]
        List of preloaded User objects
    """
    # joinedLoad must be chained because `lazy='joined'` isn't described in the orm classes relationship
    # -> but this missing statement in the ORM description permit to preload OR NOT
    statement = select(User) \
                .where(User.deleted_at == None) \
                .options(joinedload(User.objects), 
                         joinedload(User.session_logs))
    if isinstance(limit, int):
        statement = statement.limit(limit)
    results = await db.execute(statement)
    results = results.unique().scalars().all()
    return results


async def get_preloaded_user(db: AsyncSession, user_id: uuid.UUID):
    """
    Get a preloaded user by ID with associated floats, f_rules, and p_rules

    Parameters
    ----------
    db : AsyncSession
        Database async session
    user_id : uuid.UUID
        ID of the user

    Returns
    -------
    User | None
        Preloaded User object if found, None otherwise
    """
    # joinedLoad must be chained because `lazy='joined'` isn't described in the orm classes relationship
    # -> but this missing statement in the ORM description permit to preload OR NOT
    statement = select(User) \
                .where(and_(User.id == user_id, User.deleted_at == None)) \
                .options(joinedload(User.objects), 
                         joinedload(User.session_logs))
    result = await db.execute(statement)
    result = result.unique().scalar_one_or_none()
    return result


async def get_preloaded_user_by_email(db: AsyncSession, email: str):
    """
    Get a preloaded user by email with associated floats, f_rules, and p_rules

    Parameters
    ----------
    db : AsyncSession
        Database async session
    email : str
        Email address of the user

    Returns
    -------
    User | None
        Preloaded User object if found, None otherwise
    """
    # joinedLoad must be chained because `lazy='joined'` isn't described in the orm classes relationship
    # -> but this missing statement in the ORM description permit to preload OR NOT
    statement = select(User) \
                .where(and_(User.email == email, User.deleted_at == None)) \
                .options(joinedload(User.objects), 
                         joinedload(User.session_logs))
    result = await db.execute(statement)
    result = result.unique().scalar_one_or_none()
    return result


async def get_user_by_email(db: AsyncSession, email: str):
    """
    Get a user by email

    Parameters
    ----------
    db : AsyncSession
        Database async session
    email : str
        Email address of the user

    Returns
    -------
    User | None
        User object if found, None otherwise
    """
    statement = select(User) \
                .where(and_(User.email == email, User.deleted_at == None))
    result = await db.execute(statement)
    result = result.unique().scalar_one_or_none()
    return result


async def create_user(db: AsyncSession, db_user: User):
    """
    Create a new user

    Parameters
    ----------
    db : AsyncSession
        Database async session
    db_user : User
        User object to create

    Returns
    -------
    User
        Created User object
    """
    db.add(db_user)
    await db.commit()
    await db.refresh(db_user)
    return db_user


async def update_user(db: AsyncSession, db_user: User):
    """
    Update an existing user

    Parameters
    ----------
    db : AsyncSession
        Database async session
    db_user : User
        User object to update

    Returns
    -------
    User
        Updated User object
    """
    await db.commit()
    await db.refresh(db_user)
    return db_user


async def soft_delete_user(db: AsyncSession, db_user: User):
    """
    Soft delete a user by setting deleted_at

    Parameters
    ----------
    db : AsyncSession
        Database async session
    db_user : User
        User object to soft delete

    Returns
    -------
    User
        Soft deleted User object
    """
    setattr(db_user, 'deleted_at', datetime.now(tz=timezone.utc))
    await db.commit()
    return db_user


## OBJECTS ##

async def get_objects(db: AsyncSession, limit: int | None = None):
    """
    Retrieves objects from the database.

    Parameters
    ----------
    db : AsyncSession
        The asynchronous database session.
    limit : int, optional
        The maximum number of objects to retrieve, by default None.

    Returns
    -------
    list[Object]
        A list of retrieved objects.

    """
    statement = select(Object) \
                .where(Object.deleted_at == None)
    if isinstance(limit, int):
        statement = statement.limit(limit)
    results = await db.execute(statement)
    results = results.unique().scalars().all()
    return results


async def get_object(db: AsyncSession, object_id: uuid.UUID):
    """
    Retrieves a specific object from the database by ID.

    Parameters
    ----------
    db : AsyncSession
        The asynchronous database session.
    object_id : UUID
        The ID of the object to retrieve.

    Returns
    -------
    Object
        The retrieved object, if found, otherwise None.

    """
    statement = select(Object) \
                .where(and_(Object.id == object_id, Object.deleted_at == None))
    result = await db.execute(statement)
    result = result.unique().scalar_one_or_none()
    return result


async def get_preloaded_objects(db: AsyncSession, limit: int | None = None):
    """
    Retrieves preloaded objects from the database.

    Parameters
    ----------
    db : AsyncSession
        The asynchronous database session.
    limit : int, optional
        The maximum number of objects to retrieve, by default None.

    Returns
    -------
    list[Object]
        A list of retrieved preloaded objects.

    """
    # joinedLoad must be chained because `lazy='joined'` isn't described in the orm classes relationship
    # -> but this missing statement in the ORM description permit to preload OR NOT
    statement = select(Object) \
                .where(Object.deleted_at == None) \
                # .options(
                #     joinedload(Object.f_rules),
                #     joinedload(Object.p_rules)
                # )
    if isinstance(limit, int):
        statement = statement.limit(limit)
    results = await db.execute(statement)
    results = results.unique().scalars().all()
    return results


async def get_preloaded_object(db: AsyncSession, object_id: uuid.UUID):
    """
    Retrieves a specific preloaded object from the database by ID.

    Parameters
    ----------
    db : AsyncSession
        The asynchronous database session.
    object_id : UUID
        The ID of the object to retrieve.

    Returns
    -------
    Object
        The retrieved preloaded object, if found, otherwise None.

    """
    # joinedLoad must be chained because `lazy='joined'` isn't described in the orm classes relationship
    # -> but this missing statement in the ORM description permit to preload OR NOT
    statement = select(Object) \
                .where(and_(Object.id == object_id, Object.deleted_at == None)) \
                # .options(
                #     joinedload(Object.f_rules),
                #     joinedload(Object.p_rules)
                # )
    result = await db.execute(statement)
    result = result.unique().scalar_one_or_none()
    return result


async def get_preloaded_object_by_name(db: AsyncSession, name: str):
    """
    Retrieves a specific preloaded object from the database by name.

    Parameters
    ----------
    db : AsyncSession
        The asynchronous database session.
    name : str
        The name of the object to retrieve.

    Returns
    -------
    Object
        The retrieved preloaded object, if found, otherwise None.

    """
    # joinedLoad must be chained because `lazy='joined'` isn't described in the orm classes relationship
    # -> but this missing statement in the ORM description permit to preload OR NOT
    statement = select(Object) \
                .where(and_(Object.name == name, Object.deleted_at == None)) \
                # .options(
                #     joinedload(Object.f_rules),
                #     joinedload(Object.p_rules)
                # )
    result = await db.execute(statement)
    result = result.unique().scalar_one_or_none()
    return result


async def get_object_by_name(db: AsyncSession, name: str):
    """
    Retrieves a specific object from the database by name.

    Parameters
    ----------
    db : AsyncSession
        The asynchronous database session.
    name : str
        The name of the object to retrieve.

    Returns
    -------
    Object
        The retrieved object, if found, otherwise None.

    """
    statement = select(Object) \
                .where(and_(Object.name == name, Object.deleted_at == None))
    result = await db.execute(statement)
    result = result.unique().scalar_one_or_none()
    return result


async def create_object(db: AsyncSession, db_object: Object, db_user: User):
    """
    Creates a new object in the database.

    Parameters
    ----------
    db : AsyncSession
        The asynchronous database session.
    db_object : Object
        The object to be created.
    db_user : User
        The user performing the operation.

    Returns
    -------
    Object
        The created object.

    """
    db_object.created_by = db_user.id
    db_object.updated_by = db_user.id
    db_object.user_id = db_user.id

    db.add(db_object)
    await db.commit()
    await db.refresh(db_object)
    return db_object


async def update_object(db: AsyncSession, db_object: Object, db_user: User):
    """
    Updates an existing object in the database.

    Parameters
    ----------
    db : AsyncSession
        The asynchronous database session.
    db_object : Object
        The object with updated data.
    db_user : User
        The user performing the operation.

    Returns
    -------
    Object
        The updated object.

    """
    db_object.updated_by = db_user.id

    await db.commit()
    await db.refresh(db_object)
    return db_object


async def soft_delete_object(db: AsyncSession, db_object: Object, db_user: User):
    """
    Soft deletes an object in the database.

    Parameters
    ----------
    db : AsyncSession
        The asynchronous database session.
    db_object : Object
        The object to be soft deleted.
    db_user : User
        The user performing the operation.

    Returns
    -------
    Object
        The soft deleted object.

    """
    setattr(db_object, 'deleted_at', datetime.now(tz=timezone.utc))
    setattr(db_object, 'updated_by', db_user.id)
    setattr(db_object, 'deleted_by', db_user.id)
    await db.commit()
    return db_object


## SESSIONS ##

async def get_session_logs(db: AsyncSession, limit: int | None = None):
    """
    Retrieves session_logs from the database.

    Parameters
    ----------
    db : AsyncSession
        The asynchronous database session.
    limit : int, optional
        The maximum number of session_logs to retrieve, by default None.

    Returns
    -------
    list[SessionLog]
        A list of retrieved session_logs.

    """
    statement = select(SessionLog) \
                .where(SessionLog.deleted_at == None)
    if isinstance(limit, int):
        statement = statement.limit(limit)
    results = await db.execute(statement)
    results = results.unique().scalars().all()
    return results


async def get_session_log(db: AsyncSession, session_log_id: uuid.UUID):
    """
    Retrieves a specific session_log from the database by ID.

    Parameters
    ----------
    db : AsyncSession
        The asynchronous database session.
    session_log_id : UUID
        The ID of the session_log to retrieve.

    Returns
    -------
    SessionLog
        The retrieved session_log, if found, otherwise None.

    """
    statement = select(SessionLog) \
                .where(and_(SessionLog.id == session_log_id, SessionLog.deleted_at == None))
    result = await db.execute(statement)
    result = result.unique().scalar_one_or_none()
    return result


async def get_preloaded_session_logs(db: AsyncSession, limit: int | None = None):
    """
    Retrieves preloaded session_logs from the database.

    Parameters
    ----------
    db : AsyncSession
        The asynchronous database session.
    limit : int, optional
        The maximum number of session_logs to retrieve, by default None.

    Returns
    -------
    list[SessionLog]
        A list of retrieved preloaded session_logs.

    """
    # joinedLoad must be chained because `lazy='joined'` isn't described in the orm classes relationship
    # -> but this missing statement in the ORM description permit to preload OR NOT
    statement = select(SessionLog) \
                .where(SessionLog.deleted_at == None) \
                # .options(
                #     joinedload(SessionLog.f_rules),
                #     joinedload(SessionLog.p_rules)
                # )
    if isinstance(limit, int):
        statement = statement.limit(limit)
    results = await db.execute(statement)
    results = results.unique().scalars().all()
    return results


async def get_preloaded_session_log(db: AsyncSession, session_log_id: uuid.UUID):
    """
    Retrieves a specific preloaded session_log from the database by ID.

    Parameters
    ----------
    db : AsyncSession
        The asynchronous database session.
    session_log_id : UUID
        The ID of the session_log to retrieve.

    Returns
    -------
    SessionLog
        The retrieved preloaded session_log, if found, otherwise None.

    """
    # joinedLoad must be chained because `lazy='joined'` isn't described in the orm classes relationship
    # -> but this missing statement in the ORM description permit to preload OR NOT
    statement = select(SessionLog) \
                .where(and_(SessionLog.id == session_log_id, SessionLog.deleted_at == None)) \
                # .options(
                #     joinedload(SessionLog.f_rules),
                #     joinedload(SessionLog.p_rules)
                # )
    result = await db.execute(statement)
    result = result.unique().scalar_one_or_none()
    return result


async def get_preloaded_session_log_by_token(db: AsyncSession, token: str):
    """
    Retrieves a specific preloaded session_log from the database by token.

    Parameters
    ----------
    db : AsyncSession
        The asynchronous database session.
    token : str
        The token of the session_log to retrieve.

    Returns
    -------
    SessionLog
        The retrieved preloaded session_log, if found, otherwise None.

    """
    # joinedLoad must be chained because `lazy='joined'` isn't described in the orm classes relationship
    # -> but this missing statement in the ORM description permit to preload OR NOT
    statement = select(SessionLog) \
                .where(and_(SessionLog.token == token, SessionLog.deleted_at == None)) \
                # .options(
                #     joinedload(SessionLog.f_rules),
                #     joinedload(SessionLog.p_rules)
                # )
    result = await db.execute(statement)
    result = result.unique().scalar_one_or_none()
    return result


async def get_session_log_by_token(db: AsyncSession, token: str):
    """
    Retrieves a specific session_log from the database by token.

    Parameters
    ----------
    db : AsyncSession
        The asynchronous database session.
    token : str
        The token of the session_log to retrieve.

    Returns
    -------
    SessionLog
        The retrieved session_log, if found, otherwise None.

    """
    statement = select(SessionLog) \
                .where(and_(SessionLog.token == token, SessionLog.deleted_at == None))
    result = await db.execute(statement)
    result = result.unique().scalar_one_or_none()
    return result


async def create_session_log(db: AsyncSession, db_session_log: SessionLog, db_user: User):
    """
    Creates a new session_log in the database.

    Parameters
    ----------
    db : AsyncSession
        The asynchronous database session.
    db_session_log : SessionLog
        The session_log to be created.
    db_user : User
        The user performing the operation.

    Returns
    -------
    SessionLog
        The created session_log.

    """
    db_user =await get_user(db=db, user_id=db_user.id)
    db_session_log.created_by = db_user.id
    db_session_log.updated_by = db_user.id
    db_session_log.user_id = db_user.id

    db.add(db_session_log)
    await db.commit()
    await db.refresh(db_session_log)
    if db_user:
        db_user.last_seen = db_session_log.last_seen
        await db.commit()
        await db.refresh(db_user)
    
    return db_session_log


async def update_session_log(db: AsyncSession, db_session_log: SessionLog, db_user: User = None):
    """
    Updates an existing session_log in the database.

    Parameters
    ----------
    db : AsyncSession
        The asynchronous database session.
    db_session_log : SessionLog
        The session_log with updated data.
    db_user : User | None
        The user performing the operation.

    Returns
    -------
    SessionLog
        The updated session_log.

    """
    if db_user:
        user_id = db_user.id
        db_user =await get_user(db=db, user_id=user_id)
        db_session_log.updated_by = user_id
    
    await db.commit()
    await db.refresh(db_session_log)
    if db_user:
        db_user.last_seen = db_session_log.last_seen
        await db.commit()
        await db.refresh(db_user)
    return db_session_log


async def soft_delete_session_log(db: AsyncSession, db_session_log: SessionLog, db_user: User):
    """
    Soft deletes an session_log in the database.

    Parameters
    ----------
    db : AsyncSession
        The asynchronous database session.
    db_session_log : SessionLog
        The session_log to be soft deleted.
    db_user : User
        The user performing the operation.

    Returns
    -------
    SessionLog
        The soft deleted session_log.

    """
    setattr(db_session_log, 'deleted_at', datetime.now(tz=timezone.utc))
    setattr(db_session_log, 'updated_by', db_user.id)
    setattr(db_session_log, 'deleted_by', db_user.id)
    await db.commit()
    return db_session_log
