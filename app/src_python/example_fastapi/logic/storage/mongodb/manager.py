import os
from urllib.parse import quote_plus
from boilerplate_shared_module.schemes_pydantic.frontend import *
from pymongo import MongoClient, AsyncMongoClient
from bson.binary import Binary, UuidRepresentation
import bson



class MongoDbManager:
    def __init__(
        self, 
        ):
        self.host, self.port, self.db_name, self.username, self.password = \
            self.parse_db_metadata()
        self.client: MongoClient = MongoClient(self.url, uuidRepresentation="standard")
        self.async_client: AsyncMongoClient = AsyncMongoClient(self.url, uuidRepresentation="standard")
        
    def parse_db_metadata(self):
        host = os.environ.get('MONGODB_HOST', None)
        port = os.environ.get('MONGODB_PORT', None)
        usn = os.environ.get('MONGODB_USER', None)
        pwd = os.environ.get('MONGODB_PASSWORD', None)
        db_name = os.environ.get('MONGODB_NAME', None)
        return host, port, db_name, usn, pwd
    
    @property
    def url(self):
        uri = f"mongodb://{quote_plus(self.username)}:{quote_plus(self.password)}@{self.host}:{self.port}"
        return uri
    
    def get_db(self, mongo_client: MongoClient | AsyncMongoClient):
        db = mongo_client.get_database(                              
            mongo_manager.db_name,                                                          
            bson.codec_options.CodecOptions(
                uuid_representation=UuidRepresentation.STANDARD,
                ),
        ) 
        return db
    


mongo_manager: MongoDbManager = MongoDbManager()