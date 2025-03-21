from shiny import ui, Inputs, Outputs, Session, reactive, App, render, module
from shiny.types import SilentException
from ..frontend_classes.objects_store import ObjectStore
from ..frontend_classes.workspace import *
from shiny.ui import *
import os
from pathlib import Path
import re
from ..modules.functional import *
from ..tools.others import *
from ..tools.shiny import *


class BaseShinyTmpl:
    """
    Base shiny template.

    Attributes
    ----------
    object_store : ObjectStore
        The object store to store user sessions and workspaces.
    app : App
        The Shiny application.
    panel_sidebar: list
        List of shiny ui components for the sidebar
    panel_main: list
        List of shiny ui components for the main layout
    
    Note
    ----
    This class initializes self attributes for a shiny class. There are 2 ways of defining the app.
    First, you can use the template, but you must define properties for `panel_main` and `panel_sidebar` that are lists of shiny ui components :
    ```
    class FjMensualReportPage(BaseShinyTmpl):
        def __init__(self, object_store: ObjectStore):
            title = 'Forfait Jour Mensual Report'
            super().__init__(object_store=object_store, 
                             title=title, 
                             panel_sidebar=self.panel_sidebar,
                             panel_main=self.panel_main
                            )
    ```
    Secondly, if you don't want the template, you can create your custom app_ui by overriding the `app_ui` property, and so you don't need to define 
    properties for `panel_main` and `panel_sidebar`.
    
    In Both cases, don't forget to define the `server` function
    """
    def __init__(self, object_store: ObjectStore, title: str, panel_sidebar: list = [], panel_main: list = []):
        self.title = title
        self.object_store: ObjectStore = object_store
        self._panel_sidebar = panel_sidebar
        self._panel_main = panel_main
        self.app = App(self.app_ui, self.server, 
                       debug=True
                       )
        # self.app.on_shutdown(self.pool.shutdown)
    
    @property
    def snake_classname(self):
        name = self.__class__.__name__
        snake_name = re.sub(r'(?<!^)(?=[A-Z])', '_', name).lower()
        return snake_name   
            
    @property
    def app_ui(self):
        """
        Defines the user interface of the application.

        Returns
        -------
        ui.Page
            The user interface.
        """
        app_ui = ui.page_fluid(  
            ui.head_content(
                ui.tags.style("""
                                .title {
                                    font-size: 24px;
                                    margin-bottom: 20px;
                                }}
                                .logo {
                                    position: absolute;
                                    top: 20px;
                                    left: 20px;
                                    height: 100px
                                }
                                .center {
                                    margin: auto;
                                    padding: 10px;
                                }
                            """)
                        ),
            ui.a(
                ui.output_image('pokapok_image', 
                                {"class": "logo"},
                                height='100px'),
                                href="https://www.pokapok.org/index.php/en/pokapok-en/",
                ),
            ui.panel_title(ui.h2(self.title),
                           {"class": "title"}
                           ),
            # ui.div(
            #     ui.input_action_button("logout_btn", "Logout", class_="btn-danger"),
            #     style="float: right;"
            # ),
            dev_module_ui(self.snake_classname),
            sync_controller_as_model_ui(self.snake_classname),
            ui.layout_sidebar(
                ui.sidebar(
                    ui.output_ui(id='username_output_text'),
                    *self.panel_sidebar,
                    ui.output_ui('dynamic_sidebar_ui'),
                    open='open',
                    position='left',
                    bg='lightblue',
                ),
                *self.panel_main,
                ui.output_ui('dynamic_main_ui'),
            ),
            ui.HTML("""
                    <script type="text/javascript">
                        Shiny.addCustomMessageHandler('redirect', function(message) {
                            window.location.href = message.url;
                        });
                    </script>
                    """),
            )
        return app_ui
    
    @property
    def panel_sidebar(self):
        return self._panel_sidebar
    
    @property
    def panel_main(self):
        return self._panel_main
    
    def server(self, input: Inputs, output: Outputs, session: Session):
        pass


class BaseServerTmpl:
    """
    This class initializes server functions for an authenticated shiny page. Call it in the shiny server function. Example:
    ```
    username = session.http_conn.session.get('user_email')
        
    user_session: Session = self.object_store.get_session(username=username)
    BaseServerTmpl(input=input, output=output, session=session, 
                   object_store=self.object_store, username=username)
    ```
    
    Note
    ----
    Make sure that an action button named `logout_button` and an output text `username_output_text` are defined in the app_ui !!
    """
    def __init__(self, input: Inputs, output: Outputs, username: str, object_store: ObjectStore, 
                 session: Session, page_name: str
                 ):
        self.input = input
        self.output = output
        self.username = username
        self.page_name = page_name
        self.object_store = object_store
        self.setup_base_callbacks(session=session)

    
    def setup_base_callbacks(self, session: Session):
        
        # username = session.http_conn.session.get('user_email')
        user_workspace: Workspace = self.object_store.get_or_create_workspace(session=session)
        # def fill_dynamic_sidebar_div():
        #     try:
        #         sidebar_info = getattr(user_workspace, f"shiny_sidebar_{self.page_name}")
        #         for widget_id, widget_dic in sidebar_info.items():
        #             func_name = widget_dic.pop('func')
        #             func = globals()[func_name]
        #             _ = widget_dic.pop('value_attr')
        #             widget = func(id=widget_id, **widget_dic)
        #             ui.insert_ui(widget, 
        #                         selector='#dynamic_sidebar_div', 
        #                         session=session)
        #     except:
        #         pass
        
        # fill_dynamic_sidebar_div()
        
        # start_refresh_token_task(username=username,
        #                          object_store=self.object_store,
        #                          session=session,
        #                          refresh_time_s=30)
        
        
        @reactive.effect
        @reactive.event(user_workspace.jobs_queue.reactive, 
                        user_workspace.running_jobs_nb.reactive)
        async def queue_jobs_consumer():
            running_jobs_nb = user_workspace.running_jobs_nb.reactive()
            jobs_queue_attr = user_workspace.jobs_queue
            if (running_jobs_nb == 0) and (not jobs_queue_attr.is_empty()):
                # Get job kwargs for the scheduler
                default_kwargs = {
                    'max_instances': 1,
                    'executor': 'threadpool'
                }
                job_kwargs = await jobs_queue_attr.get()
                kwargs = default_kwargs | job_kwargs
                
                # Schedule the job
                user_workspace.scheduler.add_job(**kwargs)
                
                # Label the job as running
                await user_workspace.increase_running_jobs(nb=1)
                
        @reactive.effect
        @reactive.event(user_workspace.responses_queue.reactive)
        async def queue_response_consumer():
            # while len(user_workspace.responses_queue.value) > 0:
            while not user_workspace.responses_queue.is_empty():
                response = await user_workspace.responses_queue.get()
                
                # Mettre à jour les éléments réactifs de Shiny i
                # await user_workspace.sync_reactives()
            
                print(f"Response from background task: {response}")
                post_process_val = post_process_func(response=response)
            # reactive.invalidate_later(2, session=session)
        
        @render.image
        def pokapok_image():
            from pathlib import Path

            image_path = "/embedded-resources/images/Pokapok.png"
            img = {"src": image_path, "height": "80px"}
            return img
        
        
        @render.ui
        def username_output_text():
            return ui.strong(f"Welcome !")
        