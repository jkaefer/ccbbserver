import os
#for Gmail smtpd api
import os.path
from googleapiclient.discovery import build
from google_auth_oauthlib.flow import InstalledAppFlow
from google.auth.transport.requests import Request
import smtplib


from flask import jsonify
from pathlib import Path
import datetime
from flask import render_template,session,flash,redirect,url_for,request,send_from_directory,make_response, current_app
from werkzeug.utils import secure_filename
#from flask_login import current_user,login_user,logout_user,login_required
#from app.models import User
from app import app
#from app import session
from app.forms import IntronSubmitForm,ExonSubmitForm,SpliceSubmitForm
#from app import db
import threading
from queue import Queue
import subprocess
import glob


import numpy as np
#import utils.tools as utils
import pandas as pd
import copy
import json


from Bio import SeqIO
import time
import pickle

from functools import update_wrapper






from celery import Celery
from kombu import Exchange, Queue
files='./static/js/files.json'

app.config.update(
    CELERY_BROKER_URL='redis://localhost:6379/1',
    CELERY_RESULT_BACKEND='redis://localhost:6379/1',
    CELERY_TASK_SERIALIZER = 'pickle',
    CELERY_RESULT_SERIALIZER = 'pickle',
    CELERY_ACCEPT_CONTENT = ['pickle']
)

def make_celery(app):
    celery = Celery(
        app.import_name,
        backend=app.config['CELERY_RESULT_BACKEND'],
        broker=app.config['CELERY_BROKER_URL'],
        #task_serializer=app.config['CELERY_SERIALIZER'],
        accept_content=['pickle'],
        
    	#task_acks_late = True,
	#worker_prefetch_multiplier = 1,
	
	#task_exchange = Exchange('default', type='direct'),
	#task_create_missing_queues = False,

	#task_queues = [
	#    Queue('intron', routing_key='intron'),
	#],

	#task_routes = {
	#    'app.routes.execIntron': {'queue': 'intron','routing_key': 'intron'}
	#}
	
    )
    celery.conf.update(app.config)

    class ContextTask(celery.Task):
        def __call__(self, *args, **kwargs):
            with app.app_context():
                return self.run(*args, **kwargs)

    celery.Task = ContextTask
    return celery
    
celery = make_celery(app)

introntask=None
splicetask=None
exontask=None
email=None
@celery.task()
def execIntron(formdict,curr_time):
	print('execIntron')
	currdir=subprocess.check_output('pwd').decode('utf-8')
	currdir=currdir[0:len(currdir)-1]

	if os.path.isdir(currdir+'/app')==True:
		os.chdir(currdir+'/app')
	#global query_id
	#global status
	
	
	query_id = str(curr_time)
	target_dir="./allJobs/Intronjob_"+query_id
	#target_dir="regsnp-intron/data/"+query_id
	
	messages = json.dumps({"main":query_id})

	print('pastmessages')
	#text=formdict['input_text']
	user=formdict['userid']
	#uploaded_file=formdict['input_file']
	inputFormat=formdict['inputFormat']
	assembly=formdict['assembly']
	#description=formdict['description']
	#email=formdict['email']
	
	print('befcmd')
	In=target_dir+"/snp_input.txt"
	Out=target_dir+"/output"
	
	
	
	
	#"regsnp_intron -f --iformat "+FORMAT+" -s ./regsnp-intron/data/regsnp_intron_web/settings.json "+INPUT+" "+OUTPUT
	cmd=["regsnp_intron -f --iformat "+inputFormat+" -s ./regsnp-intron/data/regsnp_intron_web/settings.json "+In+" "+Out]
	subprocess.check_output(cmd,shell=True)
	
	##ADDING TO THE FILES LIST############
	currjson=[]
	if os.path.exists(files):
		with open(files) as f:
			currjson=json.load(f)
	else:
		with open(files,'w') as f:
			f.write('[]')
	
	addition=json.loads('{"user":"'+str(user)+'","prg":"Intron","date":"'+str(curr_time)+'"}')
	currjson.append(addition)
	with open(files,'w') as f:
		json.dump(currjson,f)
	os.system('sed "s/.*/var files=&;/" '+str(files)+' >  ./static/js/files.js')
	os.system('cat ./static/js/files.js ./static/js/index.js > ./static/js/filesAccess.js')

@celery.task()
def execSplice(formdict,query_id):
	print('execSplice')
	currdir=subprocess.check_output('pwd').decode('utf-8')
	currdir=currdir[0:len(currdir)-1]
	if os.path.isdir(currdir+'/app')==True:
		os.chdir(currdir+'/app')
	#global query_id
	#global status
	
	
	#text=formdict['data_textbox']
	user=formdict['userid']
	#uploaded_file=formdict['data_file']
	atts=formdict['atts']
	target_dir="./allJobs/Splicejob_"+query_id
	#target_dir="regsnps-splicing/query_jobs/job_"+curr_time
	
	subprocess.run(['mkdir -p '+target_dir],shell=True)
	#subprocess.run(['mkdir('+target_dir+')'])

	
	
	#creating selected attribute file 
	with open('./allJobs/Splicejob_'+query_id+'/checked_attributes','w') as f:
	#with open('regsnps-splicing/query_jobs/job_'+curr_time+'/checked_attributes','w') as f:
		for a in atts:
			f.write(a+'\n')
	#FILE is uploaded
	uploaded_fl=''
	input_desp=''
	real_result=''
	
	target_dir='./allJobs/Splicejob_'+query_id
	uploaded_fl=target_dir+'/query'
	real_result=target_dir+'/result.html'
	input_desp='query'
	
	
	
	setsid perl ./regsnps-splicing/bg_model.pl ./allJobs/Splicejob_1/query query 1 ./allJobs/Splicejob_1/result.html
	
	cmd=['setsid perl ./regsnps-splicing/bg_model.pl '+uploaded_fl+' '+input_desp+' '+query_id+' '+real_result]
	subprocess.check_output(cmd,shell=True)
	
	
	##ADDING TO THE FILES LIST############
	currjson=[]
	if os.path.exists(files):
		with open(files) as f:
			currjson=json.load(f)
	else:
		with open(files,'w') as f:
			f.write('[]')
	
	addition=json.loads('{"user":"'+str(user)+'","prg":"Splice","date":"'+str(query_id)+'"}')
	currjson.append(addition)
	with open(files,'w') as f:
		json.dump(currjson,f)
	os.system('sed "s/.*/var files=&;/" '+str(files)+' >  ./static/js/files.js')
	os.system('cat ./static/js/files.js ./static/js/index.js > ./static/js/filesAccess.js')
	
	os.system('pwd')
	
@celery.task()
def execExon(formdict,query_id):
	print('execExon')
	currdir=subprocess.check_output('pwd').decode('utf-8')
	currdir=currdir[0:len(currdir)-1]
	
	if os.path.isdir(currdir+'/app')==True:
		os.chdir(currdir+'/app')
	
	user=formdict['userid']
	
	
	#ret=query_id
	target_dir="./allJobs/Exonjob_"+query_id
	uploaded_fl=target_dir+'/query'
	input_desp='query'

	
	cmd=['~/Downloads/jdk1.8.0_271/bin/java -Xmx2g -cp "latestRun.jar:./ExonImpact/lib_jar/*" ccbb.hrbeu.exonimpact.test.Exonimpact_for_server '+uploaded_fl+';Rscript ./ExonImpact/src/R/predict.r '+uploaded_fl+'_features.csv > '+target_dir+'/error2_'+input_desp]
	
	subprocess.check_output(cmd,shell=True)
	print('past subprocess')
	#print(subprocess.check_output('pwd',shell=True))

	
	##ADDING TO THE FILES LIST############
	currjson=[]
	if os.path.exists(files):
		with open(files) as f:
			currjson=json.load(f)
	else:
		with open(files,'w') as f:
			f.write('[]')
	#print('{"user":'+str(user)+',"prg":"Exon","date":'+str(query_id)+'}')
	#exit(1)
	addition=json.loads('{"user":"'+str(user)+'","prg":"Exon","date":"'+str(query_id)+'"}')
	currjson.append(addition)
	with open(files,'w') as f:
		json.dump(currjson,f)
	os.system('sed "s/.*/var files=&;/" '+str(files)+' >  ./static/js/files.js')
	os.system('cat ./static/js/files.js ./static/js/index.js > ./static/js/filesAccess.js')
	
	
	
	return query_id
	
@celery.task()
def add_together(a, b):
    res=subprocess.check_output('echo hello', shell=True)
    return res
    

    


@app.route('/app/<path:filename>')
def custom_static(filename):
	return send_from_directory(app.config['CUSTOM_PATH'], filename)
	
@app.route('/',methods=['GET','POST'])
#@app.route('/home',methods=['GET','POST'])
#@login_required
def index():
	
	
	#flash(session)
	#flash(session.get('username'))
	#if request.args.get['intron']:
	#	redirect(url_for('intronsub')
	#flash(session.get('loggedin'))
	#result = add_together.delay(23, 42)
	#print('Result',result.wait())
	#a=[x.split('/').pop() for x in a]
	
		
	#a=[{"fname"+str(x):x} for x in a]
	#a={"here":a}
	#b=json.dumps(a)
	#global files
	#files=jsonify(a)
	return render_template('index.html',title='Home')

'''
@app.route("/getData", methods=['GET'])
def getData():

	global files

	return files
'''
	


@app.route('/intronSub',methods=['GET','POST'])
def intronSub():
	#form = request.form
	print('intron')
	if request.method == 'POST':
		curr_time=datetime.datetime.now().strftime('%Y-%m-%d-%H-%M-%S')
		formdict={}
		text=request.form['input_text']
		formdict['userid']=request.form['userid']
		upfile=request.files['input_file']
		formdict['inputFormat']=request.form['inputFormat']
		formdict['assembly']=request.form['assembly']
		#formdict['description']=request.form['description']
		global email
		email=request.form['email']
		#formdict['email']=email
		
		query_id = str(curr_time)
		target_dir="./allJobs/Intronjob_"+query_id
		
		subprocess.run(['mkdir -p '+target_dir],shell=True)
		
		if upfile.filename=='':
			formdict['data_file']=''
			with open(target_dir+"/snp_input.txt",'w') as w:
				w.write(text)
		else:
			upfile.save(target_dir+"/"+upfile.filename)
			#removing preceding spaces else just upfile.save necessary
			with open(target_dir+"/snp_input.txt",'w') as w:
				with open(target_dir+"/"+upfile.filename,'r') as f:
					for line in f.readlines():
						w.write(line.strip(' '))
			
			formdict['data_file']=upfile.filename
		
		
		
		global introntask
		introntask=execIntron.delay(formdict,curr_time)
		
		
		
			
		return redirect(url_for('intronAftSub',query_id=query_id))
		#return redirect(url_for('intronAftSub',query_id=query_id,status2=status2.pid+2))
	return render_template('intronvariant.html')

@app.route('/intronAftSub',methods=['GET','POST'])
def intronAftSub():
	global introntask
	#print(pid)
	query_id=request.args['query_id']
	global email
	#CELERY HANDLES BLOCKING WITH WORKER SO NO NEED TO MONITOR PID###
	#pid=request.args['status']
	if introntask.ready()==True:
		
		p=None
		with open('../spec.pkl','rb') as f:
			p=pickle.load(f)
		
		#print(p)
		server = smtplib.SMTP("smtp.gmail.com", 587)
		server.ehlo()
		server.starttls()
		server.login('kaeferj1@gmail.com',p)
		msg="Your regSNPs-Intron result:\n"+str(request.host_url)+"intronRes?query_id="+str(query_id)
		submsg='Subject: {}\n\n{}'.format('regSNPs-Intron Result', msg)
		if email!='':
			server.sendmail('kaeferj1@gmail.com', str(email), submsg)
		else:
			server.sendmail('kaeferj1@gmail.com', 'kaeferj1@gmail.com', submsg)
		server.close()
		p=None
		
		return render_template('intronres.html',query_id=query_id)
	else:
		return render_template('intronsubmission.html',query_id=query_id)
	
	
@app.route('/intronRes',methods=['GET','POST'])
def intronRes():
	query_id = request.args.get('query_id')
	return render_template('intronres.html',query_id=query_id)
		


		
	
@app.route('/spliceSub',methods=['GET','POST'])
def spliceSub():
	if request.method == 'POST':
		
		#if os.path.exists('regsnps-splicing/query_jobs/')==False:
		#	subprocess.run(['mkdir -p '+target_dir],shell=True)
			#subprocess.run(['mkdir(regsnps-splicing/query_jobs/)'])
			
		query_id=datetime.datetime.now().strftime('%Y-%m-%d-%H-%M-%S')
		formdict={}
		text=request.form['input_text']
		
		upfile=request.files['input_file']
		formdict['userid']=request.form['userid']
		global email
		email=request.form['email']
		formdict['atts']=request.form.getlist('attribute')
		
		
		target_dir='./allJobs/Splicejob_'+query_id
		subprocess.run(['mkdir -p '+target_dir],shell=True)
		
		
		if upfile.filename=='':
			formdict['data_file']=''
			with open(target_dir+'/query','w') as w:
				w.write(text)
		else:
			upfile.save(target_dir+"/"+upfile.filename)
			#removing preceding spaces else just upfile.save necessary
			with open(target_dir+'/query','w') as w:
				with open(target_dir+"/"+upfile.filename,'r') as f:
					for line in f.readlines():
						w.write(line.strip(' '))
			
			formdict['data_file']=upfile.filename
			
		
		
		global splicetask
		splicetask=execSplice.delay(formdict,query_id)
		
			
			
		#print(status.check_output())
		return redirect(url_for('spliceAftSub',query_id=query_id))

	return render_template('splicevariant.html')
	
@app.route('/spliceAftSub',methods=['GET','POST'])
def spliceAftSub():
	global splicetask
	
	#print(pid)
	query_id=request.args['query_id']
	global email
	#CELERY HANDLES BLOCKING WITH WORKER SO NO NEED TO MONITOR PID###
	#pid=request.args['status']
	if splicetask.ready()==True:
		p=None
		with open('../spec.pkl','rb') as f:
			p=pickle.load(f)
		
		#print(p)
		server = smtplib.SMTP("smtp.gmail.com", 587)
		server.ehlo()
		server.starttls()
		server.login('kaeferj1@gmail.com',p)
		msg="Your regSNPs-Splicing result:\n"+str(request.host_url)+"spliceRes?query_id="+str(query_id)
		submsg='Subject: {}\n\n{}'.format('regSNPs-Splicing Result', msg)
		if email!='':
			server.sendmail('kaeferj1@gmail.com', str(email), submsg)
		else:
			server.sendmail('kaeferj1@gmail.com', 'kaeferj1@gmail.com', submsg)
		server.close()
		p=None
		return render_template('spliceres.html',query_id=query_id)
	else:
		return render_template('splicesubmission.html',query_id=query_id)
	#if status.returncode!=0:
	#	return render_template('splicesubmission.html',query_id=query_id,filename=filename)
		
	#else:
	#	return render_template('spliceres.html')
	
	
	
@app.route('/spliceRes',methods=['GET','POST'])
def spliceRes():
	query_id=request.args['query_id']
	return render_template('spliceres.html',query_id=query_id)
	
	


@app.route('/exonSub',methods=['GET','POST'])
def exonSub():
	if request.method == 'POST':
		#return render_template('exonsubmission.html',query_id='hello')
		query_id=datetime.datetime.now().strftime('%Y-%m-%d-%H-%M-%S')
		
		formdict={}
		
		formdict['userid']=request.form['userid']
		text=request.form['input_text']
		upfile=request.files['input_file']
		global email
		email=request.form['email']
		
		
		
		target_dir="./allJobs/Exonjob_"+query_id
		subprocess.run(['mkdir -p '+target_dir],shell=True)
		
		
		if upfile.filename=='':
			formdict['data_file']=''
			with open(target_dir+'/query','w') as w:
				w.write(text)
		else:
			upfile.save(target_dir+"/"+upfile.filename)
			#removing preceding spaces else just upfile.save necessary
			with open(target_dir+'/query','w') as w:
				with open(target_dir+"/"+upfile.filename,'r') as f:
					for line in f.readlines():
						w.write(line.strip(' '))
			
			formdict['data_file']=upfile.filename
		
		global exontask
		exontask=execExon.delay(formdict,query_id)
		#print(finished.ready())
		
		#result.wait()
		
		#res=execExon()
		#global query_id
		#status=1
		
		#print(status.check_output())
		#return render_template('exonvariant.html')
		return redirect(url_for('exonAftSub',query_id=query_id))
	return render_template('exonvariant.html')
		
@app.route('/exonAftSub',methods=['GET','POST'])
def exonAftSub():
	#if request.method=='POST':
	global exontask
	
	#print(pid)
	query_id=request.args['query_id']
	global email
	#CELERY HANDLES BLOCKING WITH WORKER SO NO NEED TO MONITOR PID###
	#pid=request.args['status']
	#print(exontask.ready())
	if exontask.ready()==True:
	#try:   
	#	os.kill(int(pid),0)
		#import time
		#time.sleep(2)
		print('ready')
		#print(str(request.host_url))
		p=None
		with open('../spec.pkl','rb') as f:
			p=pickle.load(f)
		
		#print(p)
		server = smtplib.SMTP("smtp.gmail.com", 587)
		server.ehlo()
		server.starttls()
		server.login('kaeferj1@gmail.com',p)
		msg="Your ExonImpact result:\n"+str(request.host_url)+"exonRes?query_id="+str(query_id)
		submsg='Subject: {}\n\n{}'.format('ExonImpact Result', msg)
		print('email',email)
		if email!='':
			server.sendmail('kaeferj1@gmail.com', str(email), submsg)
		else:
			server.sendmail('kaeferj1@gmail.com', 'kaeferj1@gmail.com', submsg)
		server.close()
		p=None
		
		return render_template('exonres.html',query_id=query_id)
	#except:
	else:
		print('not ready')
		return render_template('exonsubmission.html',query_id=query_id)
	#return render_template('exonres.html')
	#except:
		
@app.route('/exonRes',methods=['GET','POST'])
def exonRes():
	print('made it')
	#if request.method=='POST':
	#pid=request.args['status']
	#print(pid)
	query_id=request.args['query_id']
	#print(query_id)
	
	#try: 
	#	os.kill(int(pid),0)
	return render_template('exonres.html', query_id=query_id)

		



