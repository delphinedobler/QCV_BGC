from pydantic import BaseModel
from typing import Type, Any
from pydantic import BaseModel


class classproperty:
    def __init__(self, func):
        self.fget = func
    def __get__(self, instance, owner):
        return self.fget(owner)

def get_object_model(dic: dict[Type[BaseModel], Any], object: Any) -> Type[BaseModel]:
    for key, value in dic.items():
        if value == object:
            return key
    return None