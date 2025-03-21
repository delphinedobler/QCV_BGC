from starlette.requests import Request
from ..handlers.authentication import *
from ..metadata.api import get_api_meta
from boilerplate_shared_module.tools.api_infos import ApiMetadata
from fastapi import Depends, HTTPException, Response, APIRouter, Header


not_auth_router = APIRouter()


@not_auth_router.get("/login")
async def route_login(
    request: Request,
    api_meta: ApiMetadata = Depends(get_api_meta)
    ):
    """
    Route for user login.

    Parameters
    ----------
    request : Request
        Incoming request.

    Returns
    -------
    RedirectResponse
        Redirects to the login page.
    """
    response = await handler_login(
        request=request,
        api_meta=api_meta
        )
    if response.http_code != 200:
        raise HTTPException(status_code=response.http_code, detail=response.model_dump())
    # Call the redirect response stored in the data of the response capsule
    return response.data

    

@not_auth_router.get("/logout")
async def route_logout(
    request: Request,
    redirect_url: str | None = None,
    api_meta: ApiMetadata = Depends(get_api_meta)
    ):
    """
    Route for user logout.

    Parameters
    ----------
    request : Request
        Incoming request.
    redirect_url : str | None
        URL where the user will be redirected after being logged out.

    Returns
    -------
    RedirectResponse
        Redirects to the logout page.
    """
    response = await handler_logout(
        request=request,
        redirect_url=redirect_url,
        api_meta=api_meta,
        )
    if response.http_code != 200:
        raise HTTPException(status_code=response.http_code, detail=response.model_dump())
    # Call the redirect response stored in the data of the response capsule
    return response.data



@not_auth_router.get("/auth", # dependencies=[Depends(valid_access_token)]
                 )
async def route_auth(
    request: Request,
    redirect_url: str | None = None,
    api_meta: ApiMetadata = Depends(get_api_meta),
    ):
    """
    Route for user authentication.

    Parameters
    ----------
    request : Request
        Incoming request.

    Returns
    -------
    RedirectResponse
        Redirects to the authentication page.
    """
    if redirect_url:
        # Stocker la redirect_url dans les headers de la réponse ou la session
        # request.headers['X-Redirect-URL'] = redirect_url  # Non standard, mais possible
        request.session['redirect_url'] = redirect_url  # En session comme fallback
    response = await handler_auth(
        request=request,
        api_meta=api_meta,
        )
    if response.http_code != 200:
        raise HTTPException(status_code=response.http_code, detail=response.model_dump())
    # Call the redirect response stored in the data of the response capsule
    return response.data


@not_auth_router.get("/refresh_token",
                 )
async def route_refresh_token(
    request: Request,
    api_meta: ApiMetadata = Depends(get_api_meta)
    ):
    """
    Route for refreshing authentication token.

    Parameters
    ----------
    request : Request
        Incoming request.

    Returns
    -------
    RedirectResponse
        Redirects to the token refresh page.
    """
    # response = await handler_refresh_token(request=request)
    # if response.http_code != 200:
    #     raise HTTPException(status_code=response.http_code, detail=response.model_dump())
    # # Call the redirect response stored in the data of the response capsule
    # return response
    # Go through the middlewrae that automatically refreshes the token
    # return Response(status_code=204)  # No Content
    # await refresh_token(request=request, authorization=authorization)
    response = await handler_refresh_token(
        request=request,
        api_meta=api_meta
        )
    if response.http_code != 200:
        raise HTTPException(status_code=response.http_code, detail=response.model_dump())
    # Call the redirect response stored in the data of the response capsule
    return response
    
@not_auth_router.get("/auth_callback")
async def route_auth_callback(
    request: Request, 
    code: str,
    state: str,
    redirect_url: str | None = None,
    api_meta: ApiMetadata = Depends(get_api_meta)
    ):
    """
    Route for authentication callback.

    Parameters
    ----------
    request : Request
        Incoming request.
    code : str
        Authentication code.
    state: str
        State of the incoming request
    redirect_url : str | None
        URL where the user will be redirected after authentication.

    Returns
    -------
    RedirectResponse
        Redirects to the authentication callback page.
    """
    response = await handler_auth_callback(
        request=request, 
        code=code,
        state=state,
        redirect_url=redirect_url,
        api_meta=api_meta
        )
    if response.http_code != 200:
        raise HTTPException(status_code=response.http_code, detail=response.model_dump())
    # Call the redirect response stored in the data of the response capsule
    return response.data



# @not_auth_router.get("/current_user")
# async def route_get_current_user(
#     request: Request,
#     api_meta: ApiMetadata = Depends(get_api_meta)
#     ):
#     """
#     Route for retrieving current user information.

#     Parameters
#     ----------
#     request : Request
#         Incoming request.

#     Returns
#     -------
#     RedirectResponse
#         Redirects to the current user information page.
#     """
#     response = await handler_get_current_user(
#         request=request,
#         api_meta=api_meta
#         )
#     if response.http_code != 200:
#         raise HTTPException(status_code=response.http_code, detail=response.model_dump())
#     # Call the redirect response stored in the data of the response capsule
#     return response.data
