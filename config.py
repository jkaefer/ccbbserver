import os

class Config(object):
	SECRET_KEY=os.environ.get('SECRET_KEY') 
	CUSTOM_PATH='~/Documents/microblog/app/'
