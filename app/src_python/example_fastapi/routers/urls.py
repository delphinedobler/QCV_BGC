from ..metadata.api import get_api_meta
from boilerplate_shared_module.schemes_pydantic.api_metadata import ApiMetadata


api_meta: ApiMetadata = get_api_meta()
common_prefix = f"{api_meta.prefix}/{api_meta.name}/{api_meta.version}"
docs_url = f"{common_prefix}/docs"
open_api_url = f"{common_prefix}/openapi.json"