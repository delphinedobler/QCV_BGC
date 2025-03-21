import inflect
from sqlalchemy import create_engine, MetaData, Engine
import psycopg2
from sqlalchemy.engine.url import make_url, URL


def tablename2plural(tablename):
    """
    Convert a table name from singular to plural form.

    Parameters
    ----------
    tablename : str
        The singular form of the table name.

    Returns
    -------
    str
        The plural form of the table name.
    """
    #define engine for singular-plural conversion
    p = inflect.engine()
    
    splitted_string = tablename.split('_')
    splitted_string[-1] = p.plural(splitted_string[-1])
    
    return '_'.join(splitted_string)


def get_url(username: str, 
            password: str, 
            host: str, 
            port: int, 
            db_name: str, 
            asynchrone=False
            ):
    """
    Generate a URL for connecting to the database.

    Parameters
    ----------
    asynchrone : bool, optional
        Whether to generate an asynchronous URL, by default False.

    Returns
    -------
    sqlalchemy.engine.url.URL
        The URL object for connecting to the database.
    """
    # config_params = load_config(config_file)
    prefix = 'postgresql'
    if asynchrone:
        prefix += '+asyncpg'
    
    url_object = URL.create(
        prefix,
        username=username,
        password=password,  # plain (unescaped) text
        host=host,
        database=db_name,
        port=port
    )
    
    return make_url(
        # f"{prefix}://{username}:{pwd}@{db_host}:{db_port}/{db_name}"
        url_object
        )


def check_database_existence(db_name: str, postgres_conn, close_connection: bool = True):
    """
    Check if the database exists.

    Parameters
    ----------
    db_name : str
        Name of the database which existence must be checked
    postgres_conn : psycopg2.extensions.connection
        The connection to the PostgreSQL database.
    close_connection : bool, optional
        Whether to close the connection after execution, by default True.

    Returns
    -------
    bool
        True if the database exists, False otherwise.
    """        
    cur = postgres_conn.cursor()
    try:
        # Execute SQL query to check if the database exists
        cur.execute(f"SELECT datname FROM pg_database WHERE datname='{db_name}';")
        # Fetch the result
        result = cur.fetchone()
        # Check if the database exists
        if result:
            return True
        else:
            return False
    except Exception as e:
        print("Error:", e)
        return False
    finally:
        # Close the cursor
        cur.close()
        
        if close_connection:
            # Close the connection
            postgres_conn.close()


def drop_everything(engine: Engine):
    """
    Drop all tables and constraints from a database.

    Parameters
    ----------
    engine : Engine
        The SQLAlchemy engine.
    """
    from sqlalchemy.inspection import inspect
    from sqlalchemy.schema import (
        DropConstraint,
        DropTable,
        MetaData,
        Table,
        ForeignKeyConstraint,
    )

    con = engine.connect()
    trans = con.begin()
    inspector = inspect(engine)

    # We need to re-create a minimal metadata with only the required things to
    # successfully emit drop constraints and tables commands for postgres (based
    # on the actual schema of the running instance)
    meta = MetaData()
    tables = []
    all_fkeys = []

    for table_name in inspector.get_table_names():
        fkeys = []

        for fkey in inspector.get_foreign_keys(table_name):
            if not fkey["name"]:
                continue

            fkeys.append(ForeignKeyConstraint((), (), name=fkey["name"]))

        tables.append(Table(table_name, meta, *fkeys))
        all_fkeys.extend(fkeys)

    for fkey in all_fkeys:
        con.execute(DropConstraint(fkey))

    for table in tables:
        con.execute(DropTable(table))

    trans.commit()