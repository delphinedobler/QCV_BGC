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
from .page_templates.not_authenticated_page_tmpl import *
from .modules.functional import *
from .modules.example import (ui_sidebar as example_ui_sidebar, 
                              ui_main as example_ui_main,
                              server as example_server, 
                              MODULE_NAME as EXAMPLE_MODULE_NAME)
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


theme.tabulator_simple()


class NotAuthenticatedPage(BaseShinyTmpl):
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
        title = "Example of Not Authenticated shiny page"
        self.icons = {
            "trash": fa.icon_svg("trash"),
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
            ui.output_ui('sidebar_panel')
            # *self.widgets_vars_value_input(),
        ]
        return panel_sidebar

    @property
    def panel_main(self):
        panel_main = [
            ui.output_ui(id='main_panel'),
        ]
        return panel_main

    # def widgets_vars_value_input(self):
    #     widget_list = []
    #     input_world_stat_var_name = "input_world_stat_var"
    #     for k in VarsEnum.__members__.keys():
    #         widget_input_name = f"{k}_text_input"
    #         widget = ui.input_text(
    #             id=widget_input_name, label=f"Enter a value for {getattr(VarsEnum, k)}"
    #         )
    #         cond = f"input.{input_world_stat_var_name}.includes('{k}')"
    #         panel_cond = ui.panel_conditional(cond, widget)
    #         widget_list.append(panel_cond)
    #     return widget_list

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
        
        # Load Example server module
        example_server(EXAMPLE_MODULE_NAME, 
                       workspace=user_workspace)
        
        
        # ======= Add HERE your reactive vars ======= #
        # NOTE: Must be prefixed by `reactive_`
        # IMPORTANT : Init reactive values using those of workspace saved in the attr `view_elements`
                
        
        # ======= DON'T TOUCH ======= #
        # NOTE : manage on close session callback
        session_management(self.snake_classname,
                           workspace=user_workspace,
                           object_store=self.object_store,
                           global_session=session,
                           authenticated_session=False,
                           )
        
        # ======= Add here other modules import ======= #
        
        
        
        # ======= Add here your server components ======= #
        
        @render.ui
        def sidebar_panel():
            # Get sidebar ui of the example module
            return example_ui_sidebar(EXAMPLE_MODULE_NAME, 
                                      workspace=user_workspace)
        
        @render.ui
        def main_panel():
            # Get main ui of the example module
            return example_ui_main(EXAMPLE_MODULE_NAME)
                