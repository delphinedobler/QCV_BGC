# Colocation tool

## Purpose
Colocation tool allows to:
- pre-download into cache files copernicus data colocated with a full batch of in-situ observations, grouped by spatio-temporal criteria (a.k.a. medium cubes)
- download remaining missing colocated copernicus data for a specific set of observations
- output colocated data as a {"Format still to be defined - either dataset or dataframe or netCDF"}, for each observation, using colocation limits criteria (a.k.a mini-cubes)
- display colocated data

## How to
To use this tool:
1) make sure you installed the required libraries:<br>
> pip install -r requirements.txt
2) locally update you copernicus credentials (you only need it once):<br>
> import copernicusmarine <br> copernicusmarine.login()

3) update the configuration file: Colocation_cfg.py
4) run the tool:
python Colocation.py

## Configuration parameters
To be written

## inputs
To be written

## outputs
To be written


## Creating a Project from the Docker Boilerplate

The boilerplate aims to help developers initialize a project with a base setup according to their use case. Currently, it is primarily developed in Python but it is also working with R .


1. **Run Initialization Scripts:**
    The very first time :
    ```sh
    ./initfirst.sh
    make config_init
    ```

    When you nwant to for recreate the .env prototype ( because of fundamental modifications )
    ```sh
    make config_init
    ```


2. **Modify Environment Variables:**
    - Edit the `.env` files for each service.

3. **Build and Launch the Project:**
    Then, build, launch the container and enter a shel into it
    ```sh
    make dcb
    make dcud
    make dceb
    ```

4. **Start and Stop Services:**
    Depending on the structure of your project, you want to start the applications with launchers
    If there is only one main application, adapt and use :
    
    ```sh
    start-app.sh PYTHON_APP
    ```

    Other Important launchers in full stacks ( separated backends / frontend services ) :
    ```sh
    start-front.sh
    start-back.sh
    stop-front.sh
    stop-back.sh
    ```



#### Fullstack Application (Branch: dev_fullstack)

1. **Follow Initialization Instructions:**
    - Customize environment variables in `.env` files.

2. **Build and Launch the Project:**
    ```sh
    make dcb
    make services-up-fullstack
    ```

3. **Start and Stop Fullstack Application:**
    ```sh
    make dceb
    start-fullstack.sh
    stop-fullstack.sh
    ```

4. **Update Configuration:**
    - Modify `/runtime/config/config-deployment.yaml` if environment variables change.

#### Standalone Shiny Application (Branch: monolith)

1. **Initialize Services:**
    - Follow similar steps as for the `dev` branch. Don't up cache service it won't be used. Proxy service is optional.

2. **Launch the Application:**
    Go in the container :
    ```sh
    make dceb
    ```

    Then, launch the app:
    ```sh
    start-app.sh PYTHON_SHINY
    ```


### Additional Notes

#### Authentication Service

    If a Keycloak authentication service is required ( mainly fastapi backends ), ensure to place the `ServiceClient.json` file in the `dev-secrets` directory.

#### Shared Volumes, Shared Libraries

    The current boilerplate uses a shared library containing Pydantic models, API consumers, managers, and authentication middlewares. Developers may need to modify this shared library. 
    The recommended approach is to create a new git repository for the library, add it to `/libraries` folder located outside of this project, and mount a volume for this folder (environment variables `VOL_SHARED_LIB` or `DEV_EXTERNAL_LIBS` in `.env`).

    - Verify paths in `app/environment-dev.txt`.
    - Install libraries in dev mode:
        ```sh
        dev-py-dl-pckgs.py
        ```
    - You can mount log, data, and config volumes in the `.env` file to share configurations between different modules, such as a backend and frontend, or to access outputs and logs outside the container.

