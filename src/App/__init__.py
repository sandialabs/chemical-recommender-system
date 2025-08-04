# © 2024 National Technology & Engineering Solutions of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with NTESS, the U.S. Government retains certain rights in this software.
# SPDX-License-Identifier: BSD-3-Clause

from flask import Flask
from .state import AppState
from flask_socketio import SocketIO

socketio = SocketIO()

def create_app():
    app = Flask(__name__)

    # Initialize the state without a multiprocessing manager
    app.state = AppState()

    # Register blueprints
    from .routes.main import main_bp
    from .routes.search import search_bp
    from .routes.results import results_bp
    from .routes.contact import contact_bp
    from .routes.batch import batch_bp
    from .routes.api import api_bp

    app.register_blueprint(main_bp)
    app.register_blueprint(search_bp)
    app.register_blueprint(results_bp)
    app.register_blueprint(contact_bp)
    app.register_blueprint(batch_bp)
    app.register_blueprint(api_bp)

    socketio.init_app(app, cors_allowed_origins="*")

    return app
