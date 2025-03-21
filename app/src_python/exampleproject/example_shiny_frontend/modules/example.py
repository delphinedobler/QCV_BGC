from shiny import ui, Inputs, Outputs, Session, reactive, App, render, module
from ..frontend_classes.workspace import *
from shinywidgets import render_widget, bokeh_dependency, output_widget

MODULE_NAME = 'example'

        
@module.ui
def ui_sidebar(workspace: Workspace):

    sidebar = ui.div(
                    ui.input_selectize(id='input_shop_products',
                                    label='Select Some products to shop',
                                    choices=['Apple', 'Chicken',
                                                'Turkey', 'Orange'],
                                    selected=workspace.get_widget_value(widget_id='input_shop_products', module_name=MODULE_NAME),
                                    multiple=True),
                    ui.input_action_button('input_add_to_cart_btn',
                                        label='Add products to cart'),
                    ui.input_selectize('input_species', 
                                    'Select the wanted penguins species', 
                                    workspace.view_elements.example_penguins_data.data.species.unique().tolist(),
                                    selected=workspace.get_widget_value(widget_id='input_species', module_name=MODULE_NAME),
                                    multiple=True),
                    ui.input_selectize('input_islands', 
                                    'Select the wanted islands',
                                    workspace.view_elements.example_penguins_data.data.island.unique().tolist(),
                                    selected=workspace.get_widget_value(widget_id='input_islands', module_name=MODULE_NAME),
                                    multiple=True
                                    ),
                    ui.input_checkbox_group('input_sexs',
                                            'Select the penguins sex',
                                            choices={
                                                    'MALE': 'Male',
                                                    'FEMALE': 'Female',
                                                    'ALL': 'All',
                                                    },
                                            selected=workspace.get_widget_value(widget_id='input_sexs', module_name=MODULE_NAME),
                                            ),
                    ui.input_slider('input_bill_depth',
                                    'Select the penguins bill depth range (mm)',
                                    workspace.view_elements.example_penguins_data.data.bill_depth_mm.min(),
                                    workspace.view_elements.example_penguins_data.data.bill_depth_mm.max(),
                                    # value=[workspace.penguins_data.bill_depth_mm.min(),
                                    #     workspace.penguins_data.bill_depth_mm.max()],
                                    value=workspace.get_widget_value(widget_id='input_bill_depth', module_name=MODULE_NAME),
                                    # workspace.get_widget_kwargs()
                                    )
                )
    return sidebar

@module.ui
def ui_main():
    main = ui.div(
        ui.output_text(id='output_txt_cart'),
        bokeh_dependency(),
        output_widget(
            "penguins_plot",
        ),
    )
    return main

@module.server
def server(input: Inputs, output: Outputs, session: Session, 
           workspace: Workspace
            ):
    reactive_cart = reactive.value(workspace.get_widget_value(widget_id='cart', module_name=MODULE_NAME))

    @render.text
    def output_txt_cart():
        return f"Your cart contains following products : {reactive_cart()}"

    @reactive.effect
    @reactive.event(input.input_add_to_cart_btn)
    def add_products_to_cart_action():
        chosen_products = list(input.input_shop_products())
        # actual_cart = reactive_cart.get()
        new_products = workspace.add_products_to_cart(products=chosen_products)
        reactive_cart.set(new_products)

    @output
    @render_widget
    @reactive.event(
        input.input_species,
        input.input_islands,
        input.input_sexs,
        input.input_bill_depth,
    )
    def penguins_plot():
        fig = workspace.get_penguins_fig(
            bill_depth_range=input.input_bill_depth(),
            sexs=input.input_sexs(),
            islands=input.input_islands(),
            species=input.input_species()
        )
        return fig

    # ======== Events that update view elements in the controller ======== #
    
    @reactive.effect
    @reactive.event(input.input_shop_products)
    def update_ws_input_shop_products():
        element_name = workspace.get_view_element_name(name='input_shop_products',
                                                    module_name=MODULE_NAME)
        setattr(workspace.view_elements, element_name, list(input.input_shop_products()))
    
    @reactive.effect
    @reactive.event(input.input_species)
    def update_ws_input_species():
        element_name = workspace.get_view_element_name(name='input_species',
                                                    module_name=MODULE_NAME)
        setattr(workspace.view_elements, element_name, list(input.input_species()))


    @reactive.effect
    @reactive.event(input.input_islands)
    def update_ws_input_islands():
        element_name = workspace.get_view_element_name(name='input_islands',
                                                    module_name=MODULE_NAME)
        setattr(workspace.view_elements, element_name, list(input.input_islands()))


    @reactive.effect
    @reactive.event(input.input_sexs)
    def update_ws_input_sexs():
        element_name = workspace.get_view_element_name(name='input_sexs',
                                                    module_name=MODULE_NAME)
        setattr(workspace.view_elements, element_name, list(input.input_sexs()))

    @reactive.effect
    @reactive.event(input.input_bill_depth)
    def update_ws_input_bill_depth():
        element_name = workspace.get_view_element_name(name='input_bill_depth',
                                                    module_name=MODULE_NAME)
        setattr(workspace.view_elements, element_name, input.input_bill_depth())
