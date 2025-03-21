from starlette.requests import Request
from starlette.middleware.base import BaseHTTPMiddleware
import base64
import json


class AddCookiesToSessionMiddleware(BaseHTTPMiddleware):
    async def dispatch(self, request: Request, call_next):
        # Récupérer tous les cookies de la requête
        cookies = request.cookies
        # Ajouter les cookies dans la session
        for key, value in cookies.items():
            if key != 'session':
                value_json_bytes = base64.b64decode(value)

                # Convertir les octets en chaîne JSON
                value_json = value_json_bytes.decode('utf-8')

                # Charger la chaîne JSON en objet Python (dictionnaire, liste, etc.)
                value_decoded = json.loads(value_json)
                request.session[key] = value_decoded
        
        # Continuer la requête
        response = await call_next(request)
        return response