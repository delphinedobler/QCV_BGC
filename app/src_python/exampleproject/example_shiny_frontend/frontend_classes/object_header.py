from uuid import UUID, uuid4
from datetime import datetime, timezone
from pydantic import BaseModel, UUID4
from typing import Optional, Type, Any
from pydantic.alias_generators import to_camel
from ..tools.class_tools import classproperty, get_object_model
import inspect
from shiny import Session
from boilerplate_shared_module.schemes_pydantic.basic_schemes import ApiHeaderModel
from boilerplate_shared_module.tools.cookies import wrap_request_session_into_cookies
from boilerplate_shared_module.consumers.main_api.workspaces import *


# class CamelModel(BaseModel):
#     class Config:
#         alias_generator = to_camel
#         populate_by_name = True
#         from_attributes = True
#         revalidate_instances = 'always'


class ObjectHeader:
    """
    A class representing the header information of an object.

    Attributes
    ----------
    created_at : datetime
        The timestamp when the object was created.
    updated_at : datetime
        The timestamp when the object was last updated.
    deleted_at : datetime or None
        The timestamp when the object was deleted, or None if it has not been deleted.
    id : UUID
        The unique identifier of the object.
    created_by : UUID or None
        The UUID of the user who created the object, or None if unknown.
    updated_by : UUID or None
        The UUID of the user who last updated the object, or None if unknown.
    deleted_by : UUID or None
        The UUID of the user who deleted the object, or None if not deleted.
    """
    def __init__(self):
        """
        Initializes an instance of the ObjectHeader class.
        """
        self._created_at = datetime.now(timezone.utc)
        self._updated_at = self._created_at
        self.id: UUID = uuid4()
        self._archived: bool = False
        self._deleted_at: datetime | None = None
        self._archived_at: datetime | None = None
        self._active: bool = True
        self.locked: bool = False
        self.created_by: UUID | None = None
        self.updated_by: UUID | None = None
        self.deleted_by: UUID | None = None
        self.column_names_mapper = {
            'id': {
                'user_friendly_name': 'Id',
                'attr': 'id'
                },
            'created_at': {
                'user_friendly_name': 'Created At',
                'attr': 'created_at'
                },
            'updated_at': {
                'user_friendly_name': 'Updated At',
                'attr': 'updated_at'
                },
            'archived_at': {
                'user_friendly_name': 'Archived Date',
                'attr': 'archived_at'
                },
            'archived': {
                'user_friendly_name': 'Archived',
                'attr': 'archived'
                },
        }
    
    @classmethod       
    def get_model_from_json(cls, path: str | None = None):
        if path is None:
            return None
        with open(path, 'r') as f:
            json_str = f.read()
        model = get_object_model(cls.models_mapper, cls)
        if model is not None:
            m = model.model_validate_json(json_data=json_str)
            return m
    
    def set_attrs_from_model(self, model):
        if model is not None:
            for attr_name in model.model_fields.keys():
                model_val = getattr(model, attr_name)
                # Look if a matching class object is found for the model attr and convert when one is found
                # matching_type = self.models_mapper.get(type(model_val))
                # val = matching_type.from_model(model=model_val) if matching_type else model_val
                # setattr(self, attr_name, val)
                converted_val = self._convert_model_value(model_val)
                setattr(self, attr_name, converted_val)    
    
    def _convert_model_value(self, value):
        """
        Méthode générique pour convertir les objets modélisés, y compris ceux
        qui sont dans des wrappers (listes, tuples, dictionnaires).
        """
        if isinstance(value, list):
            return [self._convert_model_value(v) for v in value]
        elif isinstance(value, tuple):
            return tuple(self._convert_model_value(v) for v in value)
        elif isinstance(value, dict):
            return {k: self._convert_model_value(v) for k, v in value.items()}
        else:
            # Si c'est un modèle, on cherche la classe correspondante dans models_mapper
            matching_type = self.models_mapper.get(type(value))
            if matching_type:
                return matching_type.from_model(model=value)
            return value
    
    def save_as_model(
        self, 
        session: Session | None = None,
        path: str| None = None,
        ):
        model_cls = get_object_model(self.models_mapper, self.__class__)
        model = model_cls.model_validate(self)
        if path is not None:
            if path.endswith('.json'):
                json_str = model.model_dump_json()
                with open(path, 'w') as f:
                    f.write(json_str)
            else:
                print(NotImplementedError(f"{path} format isn't supported for exporting the controller as model"))
        else:
            cookies = wrap_request_session_into_cookies(
                http_session=session.http_conn.session,
            )
            response = request_put_workspace(
                workspace=model,
                workspace_id=str(self.id),
                cookies=cookies,
                )
            if response.error is not None:
                raise ValueError(response.error)
    
    @classmethod
    def get_init_kwargs_from_model(cls, model):
        init_signature = inspect.signature(cls.__init__)
        # init_args_with_types = {param.name: param.annotation for param in init_signature.parameters.values() if param.name != 'self'}
        init_args = [param for param in init_signature.parameters.keys() if param != 'self']
        if 'kwargs' in init_args:
            kwargs = model.dict()
        else:
            kwargs = model.dict(include=init_args)
        return kwargs
        
    
    @classmethod
    def from_model(cls, **kwargs):
        """
        Initialize an instance using a model
        
        Keyword Arguments
        -----------------
        - For the first option (initializing the todo list item using an existing one saved in a json file):
            model_path : str | None
                The path of the json file describing an existing model of the item
        - For the second option (initialisation using an instance of the model)
            model : Type[BaseModel] | None
                An instance of the pydantic model
        Other keyword arguments can be given when needed to initialize the class, but not stored in the model
        """
        model = kwargs.pop('model', None)
        model_path = kwargs.pop('model_path', None)
        if model is None:
            model = cls.get_model_from_json(path=model_path)
        if model is not None:
            try:
                # Get kwargs from the model attributes
                init_kwargs = cls.get_init_kwargs_from_model(model=model)
                
                # Add kwargs given as entry
                init_kwargs |= kwargs
                
                # Create an instance
                instance = cls(**init_kwargs)
                
                # Override instance attrs using model
                instance.set_attrs_from_model(model=model)
                return instance
            except Exception as e:
                raise e
        else:
            raise ValueError(f"{cls} couldn' be loaded using a model. Please, ensure that arguments contain a model or a model_path")
    
    ## Important for parsing: Mapper between pydantic models and class objects ##
    @classproperty
    def models_mapper(cls) -> dict[Type[BaseModel], Any]:
        return {ApiHeaderModel: cls}
    
    def add_column_prefix(self, prefix: str, mapper_name: str='column_names_mapper'):
        user_friendly_prefix = prefix.replace('_', ' ').capitalize() 
        mapper = getattr(self, mapper_name)
        updated_mapper = {}
        for key, value in mapper.items():
            updated_key = f"{prefix}_{key}"
            updated_value = {k: f"{user_friendly_prefix} {v}" if k == 'user_friendly_name' else v for k, v in value.items()}
            updated_mapper[updated_key] = updated_value
        setattr(self, mapper_name, updated_mapper)
        
    @property
    def archived(self):
        return self._archived
    
    @archived.setter
    def archived(self, value: bool):
        self._archived = value
        
    @property
    def date_format(self):
        """
        Returns the date format used for timestamps.

        Returns
        -------
        str
            The date format string.
        """
        return '%Y-%m-%dT%H:%M:%S.%f'
    
    @property
    def created_at(self):
        """
        Get the timestamp when the object was last updated.

        Returns
        -------
        datetime
            The timestamp of the last update.
        """
        return self._created_at

    @created_at.setter
    def created_at(self, value: datetime):
        """
        Set the timestamp when the object was last updated.

        Parameters
        ----------
        value : datetime
            The new timestamp of the last update.
        """
        self._created_at = value
    
    @property
    def updated_at(self):
        """
        Get the timestamp when the object was last updated.

        Returns
        -------
        datetime
            The timestamp of the last update.
        """
        return self._updated_at

    @updated_at.setter
    def updated_at(self, value: datetime):
        """
        Set the timestamp when the object was last updated.

        Parameters
        ----------
        value : datetime
            The new timestamp of the last update.
        """
        self._updated_at = value
    
    @property
    def deleted_at(self):
        """
        Get the timestamp when the object was deleted, if applicable.

        Returns
        -------
        datetime or None
            The timestamp of deletion or None if not deleted.
        """
        return self._deleted_at

    @deleted_at.setter
    def deleted_at(self, value: datetime):
        """
        Set the timestamp when the object was deleted.

        Parameters
        ----------
        value : datetime
            The timestamp of deletion.
        """
        self._deleted_at = value
    
    @property
    def archived_at(self):
        return self._archived_at

    @archived_at.setter
    def archived_at(self, value: datetime):
        self._archived_at = value
    
    @property
    def active(self):
        return self._active

    @active.setter
    def active(self, value: bool):
        self._active = value

    @property
    def locked(self):
        return self._locked

    @locked.setter
    def locked(self, value: bool):
        self._locked = value
