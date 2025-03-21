import yaml
from pathlib import Path
import os

from boilerplate_shared_module.schemes_pydantic.api_metadata import ApiMetadata


api_config_file = os.path.join(str(Path(__file__).resolve().parent), 'api_config.yaml')


def str_to_bool(s):
    return s.lower() in ['true', '1', 't', 'y', 'yes']


REVERSE_PROXY: bool = str_to_bool(os.environ.get('NETWORK_REVERSE_PROXY'))
HTTP_SECURED: bool = False


def get_api_metadata(api_config_file: str = api_config_file):
    """
    Get the metadata of an API.

    Args:
        api_type (ApiTypeEnum): Type of the Api (example: ingester, main_api, ...).

    Returns:
        ApiMetadata: The metadata of the API.
    """
    with open(api_config_file, 'r') as file:
        data = yaml.safe_load(file)
        data = {
            'host': os.environ.get('NETWORK_HOSTNAME'),
            'port': None if REVERSE_PROXY else os.environ.get('NETWORK_BACK_EXPOSED_PORT'),
            'http_secured': HTTP_SECURED,
            'reverse_proxy': REVERSE_PROXY,
            'prefix_base': '/api',
        } | data
        meta_model = ApiMetadata(**data)
    return meta_model
