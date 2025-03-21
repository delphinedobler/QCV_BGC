import os


UNAUTHENTICATED_USERNAME = 'anonymous'


def str_to_bool(s):
    return s.lower() in ['true', '1', 't', 'y', 'yes']


def get_user_data_dir(username: str):
    dir = os.path.join('/', 'runtime', 'data', username)
    #Ensure existence of the dir
    os.makedirs(dir, exist_ok=True)
    return dir


def get_workspace_file_path(username: str):
    extension = 'json'
    filename = f"workspace.{extension}"
    user_dir = get_user_data_dir(username=username)
    return os.path.join(user_dir, filename)