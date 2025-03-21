from fastapi import HTTPException, Path, Query, APIRouter
from boilerplate_shared_module.schemes_pydantic.basic_schemes import ResponseCapsule
from ..metadata.api import get_api_meta
from boilerplate_shared_module.tools.api_infos import ApiMetadata, api_meta_to_rc_dic
from boilerplate_shared_module.consumers.tools import get_app_route, ApiTypeEnum
from ..handlers.objects import *
import uuid
from .users import get_authenticated_user


not_auth_router = APIRouter()


@not_auth_router.get("/ping", response_model=ResponseCapsule)
async def route_read_ping(
    request: Request,
    api_meta: ApiMetadata = Depends(get_api_meta),
    ):
    """
    Make a ping on to verify if a non authenticated route can be reached

    Returns
    -------
    ResponseCapsule
        The response capsule containing the reponse.

    """
    api_meta_dic = api_meta_to_rc_dic(api_meta=api_meta)
    login_page_url = get_app_route(route='login', api_type=ApiTypeEnum.main)
    return ResponseCapsule(data=f'PING OK, login url is located at {login_page_url}', **api_meta_dic)