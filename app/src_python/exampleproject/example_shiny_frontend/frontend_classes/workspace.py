from shiny import Session, ui
import asyncio
from .object_header import ObjectHeader
from .reactive_templates import ReactAttrTmpl, ReactQueueTmpl
from apscheduler.executors.pool import ProcessPoolExecutor
from pydantic import BaseModel, Field
from typing import Optional, Callable, Type, Any
from apscheduler.schedulers.asyncio import AsyncIOScheduler
from apscheduler.events import EVENT_JOB_EXECUTED, EVENT_JOB_ERROR
from pytz import utc
from apscheduler.executors.asyncio import AsyncIOExecutor
import logging
import sys
from typing import Optional, Callable
from datetime import datetime, timezone
from ..tools.errors import display_error
from ..tools.others import *
from .models import *
from boilerplate_shared_module.schemes_pydantic.frontend import ViewElements, DataFrameModel, WorkSpaceModel
from bokeh.sampledata.penguins import data as penguins_data
from bokeh.plotting import figure
from bokeh.transform import factor_cmap, factor_mark
from ..tools.class_tools import classproperty
from boilerplate_shared_module.tools.cookies import wrap_request_session_into_cookies
    

class FrontendResponse(BaseModel):
    session: Session
    workspace: "Workspace"
    notif_id: Optional[str] = None
    message: str = "OK"
    details: Optional[str] = None
    date: str = Field(default_factory=lambda: datetime.now(timezone.utc).isoformat())
    error: str | None = None
    data: Any = None
    callback: Optional[Callable] = None
    kwargs: Optional[dict] = {}

    class Config:
        arbitrary_types_allowed = True


logging.basicConfig()
logger = logging.getLogger("apscheduler")
logger.setLevel(logging.INFO)
# StreamHandler for the console
stream_handler = logging.StreamHandler(sys.stdout)
log_formatter = logging.Formatter(
    "%(asctime)s [%(processName)s: %(process)d] [%(threadName)s: %(thread)d] [%(levelname)s] %(name)s: %(message)s"
)
stream_handler.setFormatter(log_formatter)
logger.addHandler(stream_handler)


def post_process_func(response: FrontendResponse | None = None):
    if isinstance(response, FrontendResponse):
        if response.notif_id:
            msg_type = "error" if response.error else "message"
            ui.notification_show(
                response.message,
                session=response.session,
                id=response.notif_id,
                duration=None,
                close_button=True,
                type=msg_type,
            )
        if response.error is not None:
            display_error(
                session=response.session,
                title=response.message,
                error=response.error,
                error_traceback=response.details,
            )
        if response.callback:
            return response.callback(**response.kwargs)


class Workspace(ObjectHeader):
    """
    A class representing a workspace.

    This class manages data related to a user's workspace.

    Attributes
    ----------
    username : str
        The username associated with the workspace.
    """

    def __init__(self, session_model: SessionModel | None = None , **kwargs):
        """
        Initializes an instance of the Workspace class.

        Parameters
        ----------
        session_model : SessionModel
            A model that contains at least the shiny session instanciated by the server.
        
        Keyword Arguments
        -----------------
        - For the first option (creation of a new workspace using a new username):
            username : str | None
                The username associated with the workspace.
            first_name : str | None
                The first name of the user associated with the workspace.
            last_name : str | None
                The last name of the user associated with the workspace.
        - For the second option (initializing the workspace using an existing one):
            model_path : str | None
                The path of the json file describing an existing model of the workspace
        """
        # ************** Fontional attributes, DON'T TOUCH ************** #
        # ======= Templating ======= #
        super().__init__()
        
        # ======= Functional header ======= #
        self.username = kwargs.get('username', 'anonymous')
        self.first_name = kwargs.get('first_name', 'anonymous')
        self.last_name = kwargs.get('last_name', '')
        self.sessions_models: list[SessionModel] = []
        self.add_session_if_not_exists(session_model=session_model)
        
        # ======= Semaphores ======= #
        self.sync_sem = asyncio.Semaphore(1)
        
        # ======= Loops ======= #
        self.listeners_loop = asyncio.get_event_loop()
        self.jobs_loop = asyncio.new_event_loop()
        self.flush_loop = asyncio.new_event_loop()
        # self.setters_loop = asyncio.new_event_loop()

        # ======= Queues ======= #
        self.jobs_queue = ReactQueueTmpl()
        self.responses_queue = ReactQueueTmpl()

       # ======= Running jobs nb ======= #
        self.running_jobs_nb = ReactAttrTmpl(default_value=0)

        # ======= Schedulers ======= #
        self.background_scheduler = AsyncIOScheduler(
            **self.scheduler_conf(max_instances=3)
        )
        self.scheduler = AsyncIOScheduler(**self.scheduler_conf(max_instances=1))
        self.add_scheduler_listeners()
        self.start_schedulers()
        
        # ************** USER Attributes, CAN BE TOUCH ************** #
        
        # ======= Add here your workspace vars ======= #
        # Reactive vars that should not be saved
        
        # Other vars
        
        # ======= Shiny dynamic UI ======= #
        # NOTE: attributes names into the model must be the sames than widgets id
        # WARNING, these attributes must be defined into the model (file `models.py`)
        self.view_elements = ViewElements(
            example_penguins_data=DataFrameModel(data=penguins_data,
                                                    path='/runtime/data/example/penguins.csv'
                                                    ),
            example_input_bill_depth=
                [
                    None,
                    None
                ],
            )
        

    # ██████╗  █████╗ ███████╗███████╗
    # ██╔══██╗██╔══██╗██╔════╝██╔════╝
    # ██████╔╝███████║███████╗█████╗  
    # ██╔══██╗██╔══██║╚════██║██╔══╝  
    # ██████╔╝██║  ██║███████║███████╗
    # ╚═════╝ ╚═╝  ╚═╝╚══════╝╚══════╝
                                    
    
    ## Important for parsing: Mapper between pydantic models and class objects ##
    @classproperty
    def models_mapper(cls) -> dict[Type[BaseModel], Any]:
        return {
            WorkSpaceModel: cls,
            }
    
    def scheduler_conf(self, max_instances: int = 3):
        conf = {
            "executors": {
                # 'default': ThreadPoolExecutor(20),
                "default": AsyncIOExecutor(),
                "threadpool": {"type": "threadpool", "max_workers": 20},
                "processpool": ProcessPoolExecutor(),
            },
            "job_defaults": {"coalesce": False, "max_instances": max_instances},
            "timezone": utc,
        }
        return conf

    def add_scheduler_listeners(self):

        def background_job_listener(event):
            if event.exception:
                print(f"Job {event.job_id} failed: {event.exception}")
            else:
                print(f"Job {event.job_id} completed successfully")

        self.background_scheduler.add_listener(
            background_job_listener, EVENT_JOB_EXECUTED | EVENT_JOB_ERROR
        )

    def start_schedulers(self):
        self.background_scheduler.start()

    def on_shutdown(self):
        self.background_scheduler.shutdown(
            wait=False
        )  # If wait is True, it will wait for its jobs to finish their execution
    
    async def increase_counter(self, step: int = 1):
        counter_val = await self.compteur.get_objects()
        counter_val += step
        await self.compteur.set_objects(value=counter_val)
        
    
    # ███████╗███████╗███████╗███████╗██╗ ██████╗ ███╗   ██╗███████╗
    # ██╔════╝██╔════╝██╔════╝██╔════╝██║██╔═══██╗████╗  ██║██╔════╝
    # ███████╗█████╗  ███████╗███████╗██║██║   ██║██╔██╗ ██║███████╗
    # ╚════██║██╔══╝  ╚════██║╚════██║██║██║   ██║██║╚██╗██║╚════██║
    # ███████║███████╗███████║███████║██║╚██████╔╝██║ ╚████║███████║
    # ╚══════╝╚══════╝╚══════╝╚══════╝╚═╝ ╚═════╝ ╚═╝  ╚═══╝╚══════╝
                                                                  

    def get_session_model_from_shiny_sess(self, session: Session):
        for session_model in self.sessions_models:
            if session_model.shiny_sess == session:
                return session_model
        return None
    
    def add_session_if_not_exists(self, session_model: SessionModel):
        """
        Add a session to the sessions stored in the workspace.

        Parameters
        ----------
        session_model : SessionModel
            Model that contains the shiny session (and keycloak session) to add.
        """
        shiny_sess = session_model.shiny_sess
        existing_session_model = self.get_session_model_from_shiny_sess(session=shiny_sess)
        if existing_session_model is None:
            self.sessions_models.append(session_model)
        self.add_or_edit_keycloak_session(session_model=session_model)
    
    def add_or_edit_keycloak_session(self, session_model: SessionModel):
        """
        Add a session to the sessions stored in the workspace.

        Parameters
        ----------
        session : Session
            The session to add.
        """
        shiny_sess = session_model.shiny_sess
        keycloak_sess = session_model.keycloak_sess
        existing_session_model = self.get_session_model_from_shiny_sess(session=shiny_sess)
        if existing_session_model is not None:
            if keycloak_sess is not None:
                setattr(existing_session_model, 'keycloak_sess', keycloak_sess)
            

    async def remove_session(self, session: Session):
        """
        Remove a session from the store.

        Parameters
        ----------
        session : Session
            The session to remove.
        """
        existing_session_model = self.get_session_model_from_shiny_sess(session=session)
        if existing_session_model is not None:
            self.sessions_models.remove(existing_session_model)
            # TODO: Expire keycloak token
            await session.close()
    
    # ████████╗ ██████╗  ██████╗ ██╗     ███████╗
    # ╚══██╔══╝██╔═══██╗██╔═══██╗██║     ██╔════╝
    #    ██║   ██║   ██║██║   ██║██║     ███████╗
    #    ██║   ██║   ██║██║   ██║██║     ╚════██║
    #    ██║   ╚██████╔╝╚██████╔╝███████╗███████║
    #    ╚═╝    ╚═════╝  ╚═════╝ ╚══════╝╚══════╝
                                               
    
    def get_widget_value(self, widget_id: str, module_name: str | None = None):
        try:
            widget_id = self.get_view_element_name(name=widget_id, module_name=module_name)
            widget_value = getattr(self.view_elements, widget_id)
            return widget_value
        except:
            return None
    
    def get_view_element_name(self, name: str, module_name: str | None = None):
        if module_name is not None:
                name = f"{module_name}_{name}"
        return name
    
    # ███████╗██╗  ██╗ █████╗ ███╗   ███╗██████╗ ██╗     ███████╗
    # ██╔════╝╚██╗██╔╝██╔══██╗████╗ ████║██╔══██╗██║     ██╔════╝
    # █████╗   ╚███╔╝ ███████║██╔████╔██║██████╔╝██║     █████╗  
    # ██╔══╝   ██╔██╗ ██╔══██║██║╚██╔╝██║██╔═══╝ ██║     ██╔══╝  
    # ███████╗██╔╝ ██╗██║  ██║██║ ╚═╝ ██║██║     ███████╗███████╗
    # ╚══════╝╚═╝  ╚═╝╚═╝  ╚═╝╚═╝     ╚═╝╚═╝     ╚══════╝╚══════╝
                                                               
    
    def add_products_to_cart(self, products: list[str]):
        cart = self.view_elements.example_cart
        self.view_elements.example_cart = cart + products
        return self.view_elements.example_cart
    
    def get_penguins_df(self,
                        bill_depth_range: list[float],
                        sexs: list[str] | None = None,
                        species: list[str] | None = None,
                        islands: list[str] | None = None):
        # ==== Substract dataframe on corresponding criterias ==== #
        # Select the sex
        if ('ALL' in sexs) or (len(sexs) == 0):
            df = self.view_elements.example_penguins_data.data
        else:
            df = self.view_elements.example_penguins_data.data[self.view_elements.example_penguins_data.data['sex'].isin(sexs)]
        
        # Select geographic area
        if (islands is not None) and len(islands) > 0:
            df = df[df['island'].isin(islands)]
        
        # Select wanted species
        if (species is not None) and len(species) > 0:
            df = df[df['species'].isin(species)]
        
        # Select using bill depth range
        min_depth, max_depth = bill_depth_range
        df = df[(df['bill_depth_mm'] >= min_depth) & (df['bill_depth_mm'] <= max_depth)]
        
        return df
    
    def get_penguins_fig(self, bill_depth_range: list[float],
                         sexs: list[str] | None = None,
                         species: list[str] | None = None,
                         islands: list[str] | None = None ):
        
        available_species = sorted(self.view_elements.example_penguins_data.data.species.unique())
        species_markers = ['hex', 'circle_x', 'triangle']
        available_sexs = ['FEMALE', 'MALE']
        sex_colors = ['red', 'blue']

        penguins_sub_df = self.get_penguins_df(bill_depth_range=bill_depth_range,
                                               sexs=sexs,
                                               species=species,
                                               islands=islands)
        
        p = figure(title = "Penguin size", background_fill_color="#fafafa")
        p.xaxis.axis_label = 'Flipper Length (mm)'
        p.yaxis.axis_label = 'Body Mass (g)'

        p.scatter("flipper_length_mm", "body_mass_g", source=penguins_sub_df,
                legend_group="species", fill_alpha=0.4, size=12,
                marker=factor_mark('species', species_markers, available_species),
                color=factor_cmap('sex', sex_colors, available_sexs))

        p.legend.location = "top_left"
        p.legend.title = "Species"
        
        return p
    