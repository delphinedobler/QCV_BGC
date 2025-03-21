import os
from sqlalchemy.orm import sessionmaker
from sqlalchemy.ext.asyncio import create_async_engine, AsyncSession
import psycopg2
import yaml
from .tools import *
from . import orm_classes


class DatabaseManager:
    def __init__(
        self, 
        base_postgres_db_name = 'postgres'
        ):
        self.base_postgres_db_name = base_postgres_db_name
        self.host, self.port, self.db_name, self.username, self.password = \
            self.parse_db_metadata()
        self.url = get_url(username=self.username,
                            host=self.host,
                            password=self.password,
                            port=self.port,
                            db_name=self.db_name,
                            asynchrone=False
                            )
        self.async_url = self.url = get_url(username=self.username,
                            password=self.password,
                            host=self.host,
                            port=self.port,
                            db_name=self.db_name,
                            asynchrone=True
                            )
        self.engine = create_engine(self.url, connect_args={}, future=True)
        self.async_engine = create_async_engine(self.async_url, 
                                        #echo=True, 
                                        future=True)
        self.SessionLocal = sessionmaker(
            autocommit=False, autoflush=False, bind=self.engine, future=True
            )
        self.AsyncSessionLocal = sessionmaker(
            self.async_engine, class_=AsyncSession, expire_on_commit=False
        )
        
    def parse_db_metadata(self):
        host = os.environ.get('POSTGRES_DB_HOST', None)
        port = os.environ.get('POSTGRES_DB_PORT', None)
        usn = os.environ.get('POSTGRES_DB_USER', None)
        pwd = os.environ.get('POSTGRES_DB_PASSWORD', None)
        db_name = os.environ.get('POSTGRES_DB_NAME', None)
        return host, port, db_name, usn, pwd
    
    def get_connection_postgres_db(self):
        """
        Get a connection to the PostgreSQL database.

        Returns
        -------
        psycopg2.extensions.connection
            The connection to the PostgreSQL database.
        """
        conn = psycopg2.connect(
                dbname=self.base_postgres_db_name,
                user=self.username,
                password=self.password,
                host=self.host,
                port=self.port 
            )
        return conn
    
    def get_connection_db(self):
        """
        Get a connection to the database.

        Returns
        -------
        psycopg2.extensions.connection
            The connection to the other database than postgres one.
        """
        conn = psycopg2.connect(
                dbname=self.db_name,
                user=self.username,
                password=self.password,
                host=self.host,
                port=self.port 
            )
        return conn
    
    def create_database_if_not_exists(self):
        """
        Create a database if it does not exist.
        """
        
        # Creation of a cursor object to execute PostgreSQL command in a database session
        base_postgres_conn = self.get_connection_postgres_db()
        cur = base_postgres_conn.cursor()
        

        if not check_database_existence(
            postgres_conn=base_postgres_conn, 
            close_connection=False,
            db_name=self.db_name
            ):
            # # Désactivation de l'autocommit pour éviter les transactions implicites
            # base_postgres_conn.autocommit = True
            # Fermez le curseur et la connexion avant de créer la base de données
            cur.close()
            base_postgres_conn.close()

            # Réouvrez la connexion avec autocommit activé pour créer la base de données
            base_postgres_conn = self.get_connection_postgres_db()
            base_postgres_conn.autocommit = True
            cur = base_postgres_conn.cursor()
            
            # Creation of the database with SQL query
            cur.execute(f"CREATE DATABASE {self.db_name};")

            # Rétablissement de l'autocommit à sa valeur par défaut
            base_postgres_conn.autocommit = False
            # Print a success message
            print("Database creation successful !")
            
        else:
            print("Database already exists !")
        
        # Close the cursor
        cur.close()
        # Close the connection
        base_postgres_conn.close()
    
    def get_db(self):
        """
        Get a database session.

        Yields
        ------
        SessionLocal
            A SQLAlchemy session.
        """
        db = self.SessionLocal()
        try:
            yield db
        finally:
            db.close()

            
    async def async_get_db(self):
        """
        Get an asynchronous database session.

        Yields
        ------
        session
            The asynchronous database session.
        """
        async with self.AsyncSessionLocal() as db:
            # try:
            yield db
            await db.commit()
            # finally:
            #     await db.close()
    
    
    async def initialize_database(self, erase_tables_before=False):
        """
        Initialize the database.

        Parameters
        ----------
        erase_tables_before : bool, optional
            Whether to erase tables before initialization, by default False.
        """
        async with self.AsyncSessionLocal() as db:
            self.create_database_if_not_exists()
            if erase_tables_before:
                drop_everything(engine=self.engine)
                await db.commit()
        async with self.async_engine.begin() as conn:
            # Create tables
            await conn.run_sync(orm_classes.Base.metadata.create_all) #(bind=self.engine, checkfirst=True)
            
            # Create user if not exists (temporary)
            # user_email = 'yann.reynaud@pokapok.org'
            # user = await get_user_by_email(db=db, email=user_email)
            # if not user:
            #     user = User(first_name='Yann', last_name='REYNAUD', email=user_email)
            #     await create_user(db=db, db_user=user)


db_manager: DatabaseManager = DatabaseManager()