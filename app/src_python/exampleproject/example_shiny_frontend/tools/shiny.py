from shiny import Inputs
from shiny.types import SilentException


def get_input_value(input: Inputs, name: str, default=None):
        result = None
        try:
            func = getattr(input, name)
            result = func.get()
            return result
        except SilentException as e:
            return default