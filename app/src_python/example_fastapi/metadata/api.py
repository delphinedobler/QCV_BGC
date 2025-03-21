from functools import partial
from boilerplate_shared_module.tools.api_infos import get_api_metadata, ApiTypeEnum

app_type = ApiTypeEnum.main
get_api_meta = partial(get_api_metadata, api_type=app_type)