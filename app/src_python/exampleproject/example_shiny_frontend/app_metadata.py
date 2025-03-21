import os
from .tools.proxy import *
from boilerplate_shared_module.schemes_pydantic.api_metadata import *
from boilerplate_shared_module.tools.api_infos import *


def get_app_meta(software_used: SoftwareUsageEnum, app_type: ApiTypeEnum | None = None):
    if software_used == SoftwareUsageEnum.simple_shiny:

        app_meta = ApiMetadata(
            type=None,
            name=None,
            host=HOST,
            port=os.environ.get('NETWORK_FRONT_INTERNAL_PORT'),
            http_secured=False,
            reverse_proxy=REVERSE_PROXY
            )
    elif software_used == SoftwareUsageEnum.advanced_shiny:
        app_meta = get_api_metadata(api_type=app_type)
    else:
        raise NotImplementedError(f"`get_app_meta` is not yet implemented for {software_used=}")
    return app_meta


software_used = SoftwareUsageEnum.advanced_shiny
app_meta = get_app_meta(
    software_used=software_used,
    app_type=ApiTypeEnum.main_front # This line can be commented if software_used=SoftwareUsageEnum.simple_shiny
    )