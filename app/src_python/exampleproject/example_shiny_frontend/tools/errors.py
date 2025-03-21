from shiny import Session, ui
import traceback


def display_error(
    session: Session,
    title: str,
    error: Exception | str,
    error_traceback: str | None = None,
):
    """
    Displays an error in the user interface.

    Parameters
    ----------
    title : str
        The title of the error message.
    error : Union[Exception, str]
        The error message or exception.
    error_traceback : str, optional
        The traceback of the error, by default None.
    """
    if not error_traceback:
        error_traceback = traceback.format_exc()
    html_str = f"<p><strong>Description:</strong> {error}</p>\n"
    html_str += "<p><strong>Details:</strong></p>\n"
    html_str += f"<p>{error_traceback}</p>"
    msg = ui.HTML(html_str)
    m = ui.modal(
        msg,
        title=ui.strong(title),
        easy_close=True,
        footer=None,
    )
    ui.modal_show(modal=m, session=session)