from shiny import ui, Inputs, Outputs, Session, reactive, App, render, module
from ..frontend_classes.objects_store import ObjectStore
from ..frontend_classes.workspace import *
from shiny.ui import *
from ..tools.others import *
from ..tools.authentication import *
import os
import faicons as fa


# @module.server
# def dynamic_ui(input: Inputs, output: Outputs, session: Session, 
#                workspace: Workspace, page_name: str):
#     def fill_dynamic_sidebar_div():
#         try:
#             sidebar_info = getattr(workspace, f"shiny_sidebar_{page_name}")
#             for widget_id, widget_dic in sidebar_info.items():
#                 func_name = widget_dic.pop('func')
#                 func = globals()[func_name]
#                 _ = widget_dic.pop('value_attr')
#                 widget = func(id=widget_id, **widget_dic)
#                 ui.insert_ui(widget, 
#                             selector='#dynamic_sidebar_div', 
#                             session=session)
#         except:
#             pass
        
#     def fill_dynamic_main_div():
#         try:
#             sidebar_info = getattr(workspace, f"shiny_main_{page_name}")
#             for widget_id, widget_dic in sidebar_info.items():
#                 func_name = widget_dic.pop('func')
#                 func = globals()[func_name]
#                 _ = widget_dic.pop('value_attr')
#                 widget = func(id=widget_id, **widget_dic)
#                 ui.insert_ui(widget, 
#                             selector='#dynamic_main_div', 
#                             session=session)
#         except:
#             pass
        
#     fill_dynamic_sidebar_div()
#     fill_dynamic_main_div()

@module.server
def session_management(input: Inputs, output: Outputs, session: Session, 
                        object_store: ObjectStore, 
                        workspace: Workspace, 
                        global_session: Session,
                        authenticated_session: bool = False):
    
    async def on_close_session_func():
        if authenticated_session:
            object_store.cancel_refresh_token_task(session=global_session)

        await object_store.remove_workspace_session(session=global_session)
        
        if len(workspace.sessions_models) == 0:
            workspace.active = False
        
        # Save model only when user is authenticated
        if object_store.unauthenticated_username not in workspace.username:
            # Save the model in file system if simple shiny, else save into mongo
            if object_store.software_used == SoftwareUsageEnum.simple_shiny:
                save_kwargs = {'path': get_workspace_file_path(username=workspace.username)}
            else:
                save_kwargs = {'session': global_session}
            workspace.save_as_model(**save_kwargs)
        
        
    workspace.active = True
    # Begin routine of token refresh if authenticated session
    if authenticated_session == True:
        object_store.start_refresh_token_task(session=global_session,
                                              refresh_time_s=30)
        
    # Enregistrer la fonction de fin de session
    global_session.on_ended(on_close_session_func)

@module.ui
def sync_controller_as_model_ui():
    return ui.div(
                  ui.input_action_button("sync_into_model_btn", fa.icon_svg('floppy-disk'),),
                  style="float: right;"
                )


@module.server
def sync_controller_as_model_server(input: Inputs, output: Outputs, session: Session, 
                                    global_session: Session,
                                    object_store: ObjectStore,
                                    workspace: Workspace, 
                                    ):
    
    @reactive.effect
    @reactive.event(input.sync_into_model_btn)
    async def sync_react_context_action():
        # Save model only when user is authenticated
        if object_store.unauthenticated_username not in workspace.username:
            # Save the model in file system if simple shiny, else save into mongo
            if object_store.software_used == SoftwareUsageEnum.simple_shiny:
                save_kwargs = {'path': get_workspace_file_path(username=workspace.username)}
            else:
                save_kwargs = {'session': global_session}
            workspace.save_as_model(**save_kwargs)


@module.ui
def dev_module_ui():
    app_exec_mode = os.environ.get('APP_EXEC_MODE')
    if app_exec_mode.upper() == 'DEV':
        return ui.div(
            ui.h5('Developper Panel'),
            ui.output_text("username_text"),
            ui.output_text("session_id_text"),
            ui.output_text("workspace_shiny_sessions_text"),
            ui.output_text("workspace_authenticated_sessions_text"),
        )


@module.server
def dev_module_server(input: Inputs, output: Outputs, session: Session, 
                      global_session: Session,
                      workspace: Workspace):
    @output
    @render.text
    def username_text():
        return f'Your username is {workspace.username}'
    
    @output
    @render.text
    def session_id_text():
        return f'You actual session id is "{id(global_session)}"'
    
    def get_ws_shiny_sessions_id():
        sessions_ids = [id(sess.shiny_sess) for sess in workspace.sessions_models]
        return sessions_ids
    
    @reactive.poll(get_ws_shiny_sessions_id, 1)
    def ws_sessions():
        return workspace.sessions_models
    
    @output
    @render.text
    def workspace_shiny_sessions_text():
        return f'Workspace actual shiny sessions are "{[id(session_model.shiny_sess) for session_model in ws_sessions()]}"'
    
    @output
    @render.text
    def workspace_authenticated_sessions_text():
        return f'Workspace actual authenticated sessions are "{[id(session_model.shiny_sess) for session_model in ws_sessions() if session_model.keycloak_sess is not None]}"'

    
# @module.server
# def session_management_div(input: Inputs, output: Outputs, session: Session):
#         @reactive.effect
#         @reactive.event(input.logout_btn)
#         async def logout():
#             try:
#                 await logout_session(session=session)
#             except Exception as e:
#                 display_error(title="Logout Failed", error=str(e), session=session)
