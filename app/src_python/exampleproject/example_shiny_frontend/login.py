from shiny import ui, Inputs, Outputs, Session, reactive, App, render
from .frontend_classes.objects_store import ObjectStore
from boilerplate_shared_module.consumers.tools import get_app_route
from boilerplate_shared_module.tools.api_infos import *

# import os
# from boilerplate_fs_api.consumers.tools import get_backend_route
# from boilerplate_fs_api.consumers.routers import *


# # internal_port_backend = os.environ.get('INTERNAL_PORT_BACKEND')
# AUTH_URL = get_backend_route(api_type=ApiTypeEnum.main, route='auth')

class LoginPage:
    def __init__(self, 
                 object_store: ObjectStore,
                 ):
        self.object_store = object_store
        self.auth_callback_url = get_app_route(
            app_meta=self.object_store.app_meta, 
            route='/auth',
            )
        self.app = App(self.app_ui, self.server, 
                       debug=True
                       )
        
    @property
    def app_ui(self):
        login_page = ui.page_fluid(
            ui.head_content(
                ui.tags.style("""
                    body, html {
                        height: 100%;
                        margin: 0;
                        display: flex;
                        justify-content: center;
                        align-items: center;
                        background-color: #f5f5f5;  /* Couleur de fond */
                    }
                    .login-container {
                        text-align: center;
                        padding: 40px;
                        background: #ffffff;  /* Fond blanc pour le contenu centré */
                        border-radius: 10px;
                        box-shadow: 0 4px 8px rgba(0, 0, 0, 0.1);  /* Ombre autour de la boîte */
                    }
                    .login-title {
                        font-size: 24px;
                        margin-bottom: 20px;
                    }
                    .login-button {
                        margin-top: 20px;
                    }
                    .logo {
                        position: absolute;
                        top: 20px;
                        left: 20px;
                        width: 100px;
                    }
                """)
                ),
            ui.a(
                ui.output_image('pokapok_image', 
                                {"class": "logo"},
                                height='100px'),
                                href="https://www.pokapok.org/index.php/en/pokapok-en/",
                ),
            ui.div(
                {"class": "login-container"},
                ui.div("Python Boilerplate Login Page", {"class": "login-title"},),
                # ui.a(ui.input_action_button(label="SSO with Keycloak", id="sso_button", style="margin-top: 10px;"), 
                #      href=AUTH_URL, target="_blank"),
                ui.input_action_button(label="SSO with Keycloak", id="sso_button", class_="login-button"),
            ),
            ui.HTML("""
                <script type="text/javascript">
                    Shiny.addCustomMessageHandler('redirect', function(message) {
                        window.location.href = message.url;
                    });
                </script>
                """),
            title="Login"
        )
        return login_page

    
    def server(self, input: Inputs, output: Outputs, session: Session):
        """
        Configures the Shiny server for the login page.

        Parameters
        ----------
        input : Inputs
            The input elements of the application.
        output : Outputs
            The output elements of the application.
        session : Session
            The user's session.
        """
        @reactive.effect
        @reactive.event(input.sso_button)
        async def sign_in():
            backend_auth_route = get_app_route(
                app_meta=self.object_store.app_meta,
                route='auth',
                )
            # Catch redirect url if one exists in the session
            redirect_url = session.http_conn.session.get('redirect_url')
            if redirect_url is None:
                redirect_url = self.object_store.front_ROUTE_AUTH_HOMEPAGE
            backend_auth_route += f"?redirect_url={redirect_url}"
            await session.send_custom_message(
                'redirect', 
                {
                    'url': backend_auth_route,
                },
            )
        
        @render.image
        def pokapok_image():
            image_path = "/embedded-resources/images/Pokapok.png"
            img = {"src": image_path, "height": "80px"}
            return img
