from .workspace import Workspace
from shiny import Session
from .object_header import ObjectHeader
from ..tools.others import *
from .models import SessionModel
from boilerplate_shared_module.managers.authentication import *
from boilerplate_shared_module.managers.redis import RedisManager
from boilerplate_shared_module.schemes_pydantic.api_metadata import *
from boilerplate_shared_module.consumers.main_api.authentication import *
from boilerplate_shared_module.tools.cookies import wrap_request_session_into_cookies
from boilerplate_shared_module.consumers.main_api.workspaces import *


class ObjectStore(ObjectHeader):
    """
    A class representing a store for objects.

    This class manages workspaces and sessions for users.

    Attributes
    ----------
    workspaces : dict[str, Workspace]
        A dictionary mapping usernames to workspace instances.
    sessions : dict[str, Session]
        A dictionary mapping usernames to session instances.
    keycloak_config_filepath: str
        Path of JSON configuration file used for keycloak. 
        Must be given if `software_used` is `SoftwareUsageEnum.simple_shiny`
    software_used: SoftwareUsageEnum
        Type of usage of the software (whether it is a simple usage or an advanced one for example). 
        Advanced usage means that there is a backend connected with the frontend
    backend_type: ApiTypeEnum | None
        Type of the backend used for authentication. 
        Must be given if `software_used` is `SoftwareUsageEnum.advanced_shiny`
    """
    def __init__(self, 
                 front_ROUTE_AUTH_HOMEPAGE: str,
                 front_app_meta: ApiMetadata,
                 redis_manager: RedisManager | None = None,
                 keycloak_config_filepath: str = '/home/serviceuser/service/dev-secrets/ServiceClient.json',
                 software_used: SoftwareUsageEnum = SoftwareUsageEnum.simple_shiny,
                 backend_type: ApiTypeEnum | None = None,
                 ):
        """
        Initializes an instance of the ObjectStore class.
        """
        super().__init__()
        self.app_meta: ApiMetadata = front_app_meta
        self.front_ROUTE_AUTH_HOMEPAGE = front_ROUTE_AUTH_HOMEPAGE
        self.unauthenticated_username = UNAUTHENTICATED_USERNAME
        self.workspaces: dict[str, Workspace] = {}
        self.keycloak_config_filepath: str = keycloak_config_filepath
        self.software_used: SoftwareUsageEnum = software_used
        self.backend_type: ApiTypeEnum | None = backend_type
        self.keycloak_manager: KeycloakManager = self.init_keycloak_manager()
        self.redis_manager: RedisManager | None = redis_manager
        # self.sessions: dict[str, Session] = {}
        # self.refresh_token_task_dic: dict[str, Task] = {}
    
    def init_keycloak_manager(self):
        # It is mandatory to have the keycloak manager to decode the access token 
        # A simple decode function (jwt.decode) doesn't work because it needs a public ssh key to decode it
            if not os.path.exists(self.keycloak_config_filepath):
                raise FileNotFoundError(
                    f'Keycloak config file was not found at {self.keycloak_config_filepath}. Please provide a config file and / or ensure that its path is correct !'
                    )
            return KeycloakManager(
                keycloak_config_path=self.keycloak_config_filepath,
            )
    
    def add_workspace_if_not_exists(self, ws: Workspace):
        """
        Add a workspace to the store.

        Parameters
        ----------
        ws : Workspace
            The workspace to add.
        """
        username = ws.username
        # Map workspace to username (only if user is authenticated => not anonymous)
        if (username not in self.workspaces.keys()) and (self.unauthenticated_username not in username):
            self.workspaces[username] = ws
    
    
    def remove_workspace(self, username: str):
        """
        Remove a workspace from the store.

        Parameters
        ----------
        ws : Workspace
            The workspace to remove.
        """
        if username in self.workspaces.keys():
            workspace = self.workspaces.get(username)
            workspace.on_shutdown()
            self.workspaces.pop(username)
    
    def get_keycloak_session(self, session: Session) -> ApiSessionInfo | None:
        if self.software_used == SoftwareUsageEnum.simple_shiny:
            authorization = session.http_conn.session.get('authorization')
            if not authorization:
                return None
            else:
                # Retrieve keycloak session model
                access_token = authorization.get('access_token')
                refresh_token = authorization.get('refresh_token')
                keycloak_session = self.keycloak_manager.get_keycloak_session_model(
                    access_token=access_token,
                    refresh_token=refresh_token
                    )
                return keycloak_session
        elif self.software_used == SoftwareUsageEnum.advanced_shiny:
            access_token = session.http_conn.session.get('access_token')
            if not access_token:
                return None
            else:
                keycloak_session = self.redis_manager.sync_get_session_cache(
                        access_token=access_token,
                    )
                return keycloak_session
        else:
            raise NotImplementedError(
                f"`ObjectStore.get_keycloak_session_model` isn't implemented yet for {self.software_used=}"
            )
    
    def get_user_identity_dic(
        self, 
        session: Session | None = None,
        session_infos: ApiSessionInfo | None = None,
        )-> dict[str, str]:
        unauthenticated_identity_dic = {
            'username': self.unauthenticated_username,
            'first_name': self.unauthenticated_username,
            'last_name': '',
        }
        # Retrieve the access token
        if self.software_used == SoftwareUsageEnum.simple_shiny:
            if session_infos is None:
                authorization = session.http_conn.session.get('authorization')
                if not authorization:
                    access_token = None
                
                else:
                    access_token = authorization.get('access_token')
            else:
                access_token = session_infos.session.token
        elif self.software_used == SoftwareUsageEnum.advanced_shiny:
            # Retrieve the access token with the shiny session or in the session_infos
            if session_infos is None:
                access_token = session.http_conn.session.get('access_token')
            else:
                access_token = session_infos.session.token
        else:
            raise NotImplementedError(
                f"`ObjectStore.get_user_identity_dic` isn't implemented yet for {self.software_used=}"
            )
        
        # Retrieve user identity information using access token
        if access_token is not None:
            # NOTE : User information will always be available at this level because token validity has been filtered previously by middleware
            user = self.keycloak_manager.get_user_infos(
                access_token=access_token,
                verify_expiration=False,
                )
            if user is None:
                identity_dic = unauthenticated_identity_dic
            else:
                username = user.email
                first_name = user.first_name
                last_name = user.last_name
                identity_dic = {
                    'username': username,
                    'first_name': first_name,
                    'last_name': last_name,
                }
        else:
            identity_dic = unauthenticated_identity_dic
        return identity_dic
    
    def get_or_create_workspace(self, session: Session):
        """
        Get the workspace associated with the specified username.

        If the workspace does not exist, it will be created.

        Parameters
        ----------
        username : str
            The username associated with the workspace.

        Returns
        -------
        Workspace or None
            The workspace instance, or None if not found.
        """
        # authorization = session.http_conn.session.get('authorization')
        # # TODO : Think about where keycloak session is initialized (from cookies? in shiny server?)
        # # session = self.get_session(username=username)
        # if authorization is not None:
        #     # Retrieve main user information
        #     access_token = authorization.get('access_token')
            
        #     # NOTE : User information will always be available at this level because token validity has been filtered previously by middleware
        #     user = self.keycloak_manager.get_user_infos(access_token=access_token)
        #     username = user.email
        #     first_name = user.first_name
        #     last_name = user.last_name
            
        #     # Retrieve keycloak session model
        #     access_token = authorization.get('access_token')
        #     refresh_token = authorization.get('refresh_token')
        #     keycloak_session = self.keycloak_manager.get_keycloak_session_model(
        #         access_token=access_token,
        #         refresh_token=refresh_token
        #         )
        # else:
        #     username = self.unauthenticated_username
        #     first_name = self.unauthenticated_username
        #     last_name = ''
        #     keycloak_session = None
        keycloak_session_infos = self.get_keycloak_session(
            session=session,
            )
        user_identity_dic = self.get_user_identity_dic(
            session=session,
            session_infos=keycloak_session_infos,
            )
        username = user_identity_dic.get('username')
        ws = self.workspaces.get(username, None)
        keycloak_sess = None if keycloak_session_infos is None else keycloak_session_infos.session
        session_model = SessionModel(shiny_sess=session, keycloak_sess=keycloak_sess)
        if ws is None:
            # ws = Workspace(username=username, session=session)
            ws = self.init_workspace(
                session_model=session_model,
                **user_identity_dic,
                )
            self.add_workspace_if_not_exists(ws=ws)
        ws.add_session_if_not_exists(session_model=session_model)
        return ws
    
    def init_workspace(self, username: str, first_name: str, last_name: str, session_model: SessionModel):
        shiny_sess = session_model.shiny_sess
        if self.software_used == SoftwareUsageEnum.advanced_shiny:
            
            cookies = wrap_request_session_into_cookies(http_session=shiny_sess.http_conn.session)
            
            def create_ws_and_get_model():
                response = request_post_workspace(
                    workspace=WorkSpaceModel(
                        username=username,
                        first_name=first_name,
                        last_name=last_name,
                        ),
                    cookies=cookies,
                )
                return response.data
            try:
                response = request_get_workspace_by_username(
                    username=username,
                    cookies=cookies,
                    )
                ws_model = response.data
            except Exception as e:
                print(e)
                ws_model = create_ws_and_get_model()
                
            if ws_model is None:
                ws_model = create_ws_and_get_model()
                
            # if os.path.exists(file_path):
            ws = Workspace.from_model(
                model=ws_model, 
                session_model=session_model,
                )
        else:
            file_path = get_workspace_file_path(username=username)
            if os.path.exists(file_path):
                ws = Workspace.from_model(
                    model_path=file_path, 
                    session_model=session_model,
                    )
            else:
                ws = Workspace(
                    session_model=session_model, 
                    username=username,
                    first_name=first_name,
                    last_name=last_name,
                    )
        return ws
        
    async def remove_workspace_session(self, session: Session, username: str | None = None):
        if username is None:
            username = self.unauthenticated_username
        workspace = self.get_or_create_workspace(session=session)
        await workspace.remove_session(session=session)
        if not workspace.sessions_models:
            self.remove_workspace(username=username)
    
    def get_keycloak_auth_url(
        self, 
        request_state: str = '',
        ):
        if self.software_used == SoftwareUsageEnum.simple_shiny:
            auth_callback_url = get_app_route(
                app_meta=self.app_meta,
                route='auth_callback',
                )
            auth_url = self.keycloak_manager.auth_url(
                authentication_callback_url=auth_callback_url,
                request_state=request_state,
                )
            return auth_url
            
        elif self.software_used == SoftwareUsageEnum.advanced_shiny:
            auth_url = get_app_route(
                api_type=self.backend_type,
                route='auth',
                )
            return auth_url
        else:
            raise NotImplementedError(
                f"`ObjectStore.get_auth_url` isn't implemented yet for {self.software_used=}"
            )
            
    def refresh_token_in_session(self, session: Session, user_workspace: Workspace):
        # Retrieve workspace and authorization
        # user_workspace = self.get_or_create_workspace(session=session)        
        # Get session model
        session_model = user_workspace.get_session_model_from_shiny_sess(session=session)
        old_refresh_token = session_model.keycloak_sess.refresh_token
        if self.software_used == SoftwareUsageEnum.simple_shiny:
            ## In this case, we don't need to refresh cache because it doesn't exist.
            ## => Only refresh authorization in headers
            # Refresh keycloak token
            refreshed_auth = self.keycloak_manager.refresh_authorization_token(
                refresh_token=old_refresh_token,
                )
        
            # Update cookies with the new token
            session.http_conn.session['authorization'] = refreshed_auth
            
            # Update session model with new keycloak session
            keycloak_sess = self.keycloak_manager.get_keycloak_session_model(
                access_token=refreshed_auth.get('access_token'),
                refresh_token=refreshed_auth.get('refresh_token'),
                )
        elif self.software_used == SoftwareUsageEnum.advanced_shiny:
            ## In this case, cache must be refreshed with the new authorization
            ## We use the backend refresh route to refresh it
            cookies = wrap_request_session_into_cookies(
                http_session=session.http_conn.session,
                )
            try:
                response = request_get_refreshed_token(
                    cookies=cookies,
                    )
                if response.error is None:
                    # Retrieve token from the response
                    authorization = response.data
                    access_token = authorization.get('access_token')
                    
                    # Update cookies with the new token
                    session.http_conn.session['access_token'] = access_token
                    
                    # Retrieve the updated session from the cache
                    keycloak_sess: ApiSessionInfo = self.redis_manager.sync_get_session_cache(
                        access_token=access_token,
                        )
                
                else:
                    raise Exception(response.error)
            except Exception as e:
                raise e
            
        else:
            raise NotImplementedError(
                f"`ObjectStore.refresh_token_in_session` isn't implemented yet for {self.software_used=}"
            )
        setattr(session_model, 'keycloak_sess', keycloak_sess.session)
        # Update user workspace session model
        user_workspace.add_or_edit_keycloak_session(session_model=session_model)
    
    def start_refresh_token_task(self, session: Session, 
                                 refresh_time_s=30,
                                 ):
        user_workspace = self.get_or_create_workspace(session=session)
        # authorization = session.http_conn.session.get('authorization')
        # session.http_conn.session['test'] = 'test'
        # func = self.keycloak_manager.refresh_authorization_token
        # # args = [object_store, username, session]
        # kwargs = {
        #           'authorization': authorization,
        #          }
        func = self.refresh_token_in_session
        kwargs = {
            'user_workspace': user_workspace,
            'session': session,
        }
        user_workspace.background_scheduler.add_job(func=func, 
                                        kwargs=kwargs, 
                                        #  args=args,
                                        id=f'refresh_token_task_{id(session)}', 
                                        trigger='interval', 
                                        seconds=refresh_time_s)
    
    def cancel_refresh_token_task(self, session: Session):
        user_workspace = self.get_or_create_workspace(session=session)
        user_workspace.background_scheduler.remove_job(f'refresh_token_task_{id(session)}')
        
