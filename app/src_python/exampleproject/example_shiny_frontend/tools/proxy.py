from .others import *
from urllib.parse import urljoin


def is_app_reverse_proxy():
    REVERSE_PROXY = os.environ.get('NETWORK_REVERSE_PROXY', None)
    if REVERSE_PROXY is None:
        return False
    else:
        return str_to_bool(REVERSE_PROXY)
    
REVERSE_PROXY = is_app_reverse_proxy()


def get_host():
    if not REVERSE_PROXY:
        #  TODO : put back localhost
        # return 'localhost'
        return os.environ.get('NETWORK_HOSTNAME')
    else:
        return os.environ.get('NETWORK_HOSTNAME')

HOST = get_host()
