from flask import Flask, session,render_template,flash,redirect,url_for
from flask_mail import Mail

from flask_cors import CORS,cross_origin

from celery import Celery



app=Flask(__name__)



SECRET_KEY = "changeme"
SESSION_TYPE = 'filesystem'
app.config.from_object(__name__)
app.config['CUSTOM_PATH']='~/Documents/microblog/app/'
CORS(app, resources={r"/*": {"origins": "*"}})



from app import routes



