from boilerplate_shared_module.schemes_pydantic.basic_schemes import ApiHeaderModel
from pydantic import BaseModel, Field, model_validator, field_serializer
from typing import Optional, List, Annotated, Tuple
from shiny import Session
from datetime import datetime
import os
import pandas as pd
import xarray as xr
from pathlib import Path
from boilerplate_shared_module.schemes_pydantic.orm import ApiSessionLog



class KeycloakSessionModel(BaseModel):
    token: str
    refresh_token: str
    token_expiration_date: datetime


class SessionModel(BaseModel):
    shiny_sess: Session
    keycloak_sess: Optional[ApiSessionLog] = None
    
    class Config:
        arbitrary_types_allowed=True