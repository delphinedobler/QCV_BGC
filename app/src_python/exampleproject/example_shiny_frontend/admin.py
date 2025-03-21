import copy
from shiny import (
    ui,
    Inputs,
    Outputs,
    Session,
    reactive,
    render,
)
from .frontend_classes.objects_store import ObjectStore
from .frontend_classes.workspace import Workspace
from .page_templates.authenticated_page_tmpl import *
from .modules.functional import *
from pytabulator import (
    # TableOptions,
    # render_tabulator,
    # output_tabulator,
    # Tabulator,
    theme,
)
import faicons as fa
from pathlib import Path
# from boilerplate_fs_api.shared_models.data import *
import numpy as np
import matplotlib.pyplot as plt
from shinywidgets import (
    output_widget,
    render_widget,
    # render_plotly,
    # render_bokeh,
    bokeh_dependency,
)
# from boilerplate_fs_api.shared_models.data import *
from bokeh.layouts import row
from .modules.admin import (
    main_admin_page_server, 
    admin_page_ui, 
    MODULE_NAME as ADMIN_MODULE_NAME,
)

theme.tabulator_simple()


class AdminPortal(BaseShinyTmpl):
    """
    Main class for the Authenticated Page  application.

    Attributes
    ----------
    title : str
        The title of the application.
    username : str
        The currently logged-in user's username.
    object_store : ObjectStore
        The object store to store user sessions and workspaces.
    authorization : Any
        The authorization of the currently logged-in user.
    loop : asyncio.AbstractEventLoop
        The asynchronous event loop of the application.
    executor : concurrent.futures.ThreadPoolExecutor
        The executor for parallel task execution.
    app : App
        The Shiny application.
    """

    def __init__(self, object_store: ObjectStore):
        """
        Initializes an instance of the AuthenticatedPage class.

        Parameters
        ----------
        object_store : ObjectStore
            The object store to store user sessions and workspaces.
        """
        title = "QCV Workflow Manager"
        self.icons = {
                # "user": fa.icon_svg("user", "regular"),
                # "gear": fa.icon_svg("gear"),
                'home': fa.icon_svg('house'),
                'trash': fa.icon_svg('trash'),
                'plus': fa.icon_svg('plus'),
                'help': fa.icon_svg('circle-info'),
                }
        super().__init__(
            object_store=object_store,
            title=title,
            panel_sidebar=self.panel_sidebar,
            panel_main=self.panel_main,
        )

    @property
    def panel_sidebar(self):
        panel_sidebar = [
            # *self.widgets_vars_value_input(),
        ]
        return panel_sidebar

    @property
    def panel_main(self):
        panel_main = [
            ui.output_ui(id='main_panel_header'),
            ui.output_ui(id='main_panel'),
        ]
        return panel_main

    def server(self, input: Inputs, output: Outputs, session: Session):
        """
        Configures the Shiny server for the authenticated page.

        Parameters
        ----------
        input : Inputs
            The input elements of the application.
        output : Outputs
            The output elements of the application.
        session : Session
            The user's session.
        """
        # user_session: Session = self.object_store.get_session(username=username)
        user_workspace: Workspace = self.object_store.get_or_create_workspace(
            session=session
        )
        username = user_workspace.username
        BaseServerTmpl(
            input=input,
            output=output,
            username=username,
            object_store=self.object_store,
            session=session,
            page_name=self.snake_classname
        )
        
        # ======= Dynamic UI modules import ======= #
        
        
        # ======= Template basics modules import ======= #
        sync_controller_as_model_server(self.snake_classname,
                                        global_session=session, 
                                        object_store=self.object_store,
                                        workspace=user_workspace,
                                        )
        dev_module_server(self.snake_classname, global_session=session, workspace=user_workspace)
        
        
        # ======= Add HERE your reactive vars ======= #
        # NOTE: Must be prefixed by `reactive_`
        # IMPORTANT : Init reactive values using those of workspace saved in the attr `view_elements`
        
        # ======= DON'T TOUCH ======= #
        # NOTE : manage on close session callback
        session_management(self.snake_classname,
                           workspace=user_workspace,
                           object_store=self.object_store,
                           global_session=session,
                           authenticated_session=True)
        
        # ======= Add here other modules import ======= #
        main_ui = reactive.value(admin_page_ui(ADMIN_MODULE_NAME))
        
        main_admin_page_server(
            ADMIN_MODULE_NAME, 
            workspace=user_workspace,
            main_locals=locals(),
            object_store=self.object_store,
            )
        
        
        # ======= Add here your server components ======= #
        
        # @render.ui
        # def dynamic_main_ui():
        #     return ui.div()

        
        @render.ui
        def main_panel_header():
            
            return ui.div(
               style="background-color: lightblue", 
            )

            
        @render.ui
        def main_panel():
            return main_ui()
