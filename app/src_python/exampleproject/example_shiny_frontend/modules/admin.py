from shiny import ui, Inputs, Outputs, Session, reactive, App, render, module
from ..frontend_classes.objects_store import ObjectStore
from ..frontend_classes.workspace import *
from pytabulator import (
    TableOptions,
    render_tabulator,
    output_tabulator,
    Tabulator,
    theme,
)
from boilerplate_shared_module.consumers.main_api.tests import request_get_ping as request_ingester_get_ping
import asyncio

MODULE_NAME = 'admin_portal'

@module.ui
def admin_page_ui():
    
    theme.tabulator_simple()
    
    return ui.div(
            ui.output_ui('admin_page_main_ui'),
            id='admin_main',
            style='padding: 6px'
        )


@module.server
def main_admin_page_server(
    input: Inputs, output: Outputs, session: Session, 
    workspace: Workspace, 
    object_store: ObjectStore,
    main_locals: dict[str, Any],
                  ):
    main_session = main_locals.get('session')
    
    
    def get_backends_status() -> dict[str, bool]:
        
        backends_ping_func_mapper = {
            'ingester': request_ingester_get_ping
        }
        
        def get_backend_status(backend_name) -> bool:
            ping_func = backends_ping_func_mapper.get(backend_name)
            try:
                response = ping_func()
                if response.error is not None:
                    return False
                else:
                    return True
            except Exception as e:
                return False
        
        return {
            backend_name: get_backend_status(backend_name=backend_name) \
                for backend_name in backends_ping_func_mapper.keys()
                }
    
    @reactive.poll(get_backends_status, 10)
    def backends_status():
        return get_backends_status()
    
    def get_users_shiny_sessions_id():
        sessions_ids = []
        for ws in object_store.workspaces.values():
            sessions_ids += [id(sess.shiny_sess) for sess in ws.sessions_models]
        return sessions_ids
    
    @reactive.poll(get_users_shiny_sessions_id, 1)
    def users_workspaces():
        return object_store.workspaces
    
    @render.ui
    def admin_page_main_ui():

            div = ui.div(
                ui.div(
                    ui.h2("Status of backends"),
                    output_tabulator('backends_status_tabulator'),
                    ),
                ui.div(
                    ui.h2('Current active users'),
                    output_tabulator('active_sessions_tabulator'),
                ),
                style='padding: 6px',
                )
            return div
        
    
    @render_tabulator
    def backends_status_tabulator():
        
        def get_backends_status_df(backends_status_dic: dict[str, bool]):
            cols = ['backend_name', 'status']
            data = []
            for name, status in backends_status_dic.items():
                data.append([name, status])
            return pd.DataFrame(data, columns=cols)
        
        df = get_backends_status_df(backends_status_dic=backends_status())
        
        def get_backends_status_tabulator_columns(df: pd.DataFrame):
            cols = []
            col_names_mapper = {
                'backend_name': 'Backend',
                'status': 'Alive',
            }
            for col_name in df.columns:
                user_friendly_name = col_names_mapper.get(col_name, col_name)
                dic = {
                    'title': user_friendly_name,
                    'field': col_name,
                    'hozAlign': 'left',
                    'headerFilter': True
                }
                if col_name == 'status':
                    dic |= {"formatter": "tickCross"}
                cols.append(dic)
            
            return cols
        
        tabulator_columns = get_backends_status_tabulator_columns(df=df)
        
        # columns = create_columns(self.all_list_vars['all_todo_lists'](), default_filter=True)
        return Tabulator(df, 
                            table_options=TableOptions(
                                default_editor=True,
                                height=300,
                                pagination=False,
                                pagination_add_row="table",
                                # frozen_rows=6,
                                movable_rows=True,
                                # layout="fitDataTable",
                                add_row_pos="top",
                                selectable_rows=True,
                                history=True,
                                columns=tabulator_columns,
                            )
                        )
    
    
    @render_tabulator
    def active_sessions_tabulator():
        
        def get_users_sessions_df(ws_mapper: dict[str, Workspace]):
            cols = ['username', 'sessions_nb']
            data = []
            for username, workspace in ws_mapper.items():
                sessions_nb = len(workspace.sessions_models)
                data.append([username, sessions_nb])
            return pd.DataFrame(data, columns=cols)
        
        df = get_users_sessions_df(ws_mapper=users_workspaces())
        
        def get_active_sessions_tabulator_columns(df: pd.DataFrame):
            cols = []
            col_names_mapper = {
                'username': 'Active users',
                'sessions_nb': 'Current number of active sessions',
            }
            for col_name in df.columns:
                user_friendly_name = col_names_mapper.get(col_name, col_name)
                dic = {
                    'title': user_friendly_name,
                    'field': col_name,
                    'hozAlign': 'left',
                    'headerFilter': True
                }
                cols.append(dic)
            
            return cols
        
        tabulator_columns = get_active_sessions_tabulator_columns(df=df)
        
        # columns = create_columns(self.all_list_vars['all_todo_lists'](), default_filter=True)
        return Tabulator(df, 
                            table_options=TableOptions(
                                default_editor=True,
                                height=300,
                                pagination=False,
                                pagination_add_row="table",
                                # frozen_rows=6,
                                movable_rows=True,
                                # layout="fitDataTable",
                                add_row_pos="top",
                                selectable_rows=True,
                                history=True,
                                columns=tabulator_columns,
                            )
                        )