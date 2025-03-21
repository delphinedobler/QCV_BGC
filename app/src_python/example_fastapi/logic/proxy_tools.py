import os
from urllib.parse import urljoin
from ..metadata.api import get_api_meta
from boilerplate_shared_module.schemes_pydantic.api_metadata import ApiMetadata


api_meta: ApiMetadata = get_api_meta()


REVERSE_PROXY = api_meta.reverse_proxy
HOST = api_meta.host
# FRONTEND_PORT_INTERNAL = os.environ.get('INTERNAL_PORT_FRONTEND')


# def get_frontend_url_prefix(http_secured=False, reverse_proxy=REVERSE_PROXY):
#     if http_secured:
#         prefix = 'https://'
#     else:
#         prefix = 'http://'
#     prefix = urljoin(prefix, HOST)
#     if not reverse_proxy:
#         prefix += f':{FRONTEND_PORT_INTERNAL}'
#     final_prefix = prefix + '/'
#     return final_prefix
    