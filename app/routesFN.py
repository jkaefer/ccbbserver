import os
#for Gmail smtpd api
import os.path
from googleapiclient.discovery import build
from google_auth_oauthlib.flow import InstalledAppFlow
from google.auth.transport.requests import Request


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

from functools import update_wrapper




#session['logged_in']=False


from celery import Celery
#print(request.host.split(':')[0])
files='./static/js/files.json'

app.config.update(
    CELERY_BROKER_URL='redis://localhost:6379',
    CELERY_RESULT_BACKEND='redis://localhost:6379'
)

def make_celery(app):
    celery = Celery(
        app.import_name,
        backend=app.config['CELERY_RESULT_BACKEND'],
        broker=app.config['CELERY_BROKER_URL']
    )
    celery.conf.update(app.config)

    class ContextTask(celery.Task):
        def __call__(self, *args, **kwargs):
            with app.app_context():
                return self.run(*args, **kwargs)

    celery.Task = ContextTask
    return celery
    
celery = make_celery(app)

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
	result = add_together.delay(23, 42)
	print('Result',result.wait())
	#a=glob.glob("/home/pcuser/Documents/microblog/app/allJobs/*")
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
	if request.method == 'POST':
		print('MADEIT')
		curr_time=datetime.datetime.now().strftime('%Y-%m-%d-%H-%M-%S')
		query_id = "query_"+str(curr_time)
		target_dir="./allJobs/Intronjob_"+query_id
		#target_dir="regsnp-intron/data/"+query_id
		
		messages = json.dumps({"main":query_id})

		
		#if form.validate_on_submit():
		#uploaded_fl=''
		text=request.form['input_text']
		user=request.form['userid']
		print('MADEIT')
		uploaded_file=request.files['input_file']
		inputFormat=request.form['inputFormat']
		print('MADEIT')
		assembly=request.form['assembly']
		description=request.form['description']
		email=request.form['email']
		print('MADEIT')
		cmd=None
		subprocess.run(['mkdir -p '+target_dir],shell=True)
		if uploaded_file.filename!='':
			
			#subprocess.run(['mkdir '+target_dir],shell=True)
			target_file = target_dir+"/snp_input.txt"
			uploaded_file.save(target_file)
			
			
			
		else:
			
			target_file = target_dir+"/snp_input.txt"			
			with open(target_file,'w') as f:
				f.write(text)
			
		cmd=["./regsnp-intron/run.sh "+target_dir+"/snp_input.txt "+target_dir+"/output "+" "+inputFormat+" "+query_id+" "+email+" > "+target_dir+"/log 2>&1 &"]
		status2=subprocess.Popen(cmd,shell=True)
		
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
			
		return redirect(url_for('intronAftSub',query_id=query_id))
		#return redirect(url_for('intronAftSub',query_id=query_id,status2=status2.pid+2))
	return render_template('intronvariant.html')

@app.route('/intronAftSub',methods=['GET','POST'])
def intronAftSub():
	query_id=request.args['query_id']
	#pid=request.args['status2']
	#print('intronpid:'+str(pid))
	
	
	#try: 
	#	absolutePath = Path('./allJobs/Splicejob_'+query_id+'/output/snp.prediction.json').resolve(strict=True)
	#	return render_template('intronres.html',query_id=query_id)
		#return render_template('splicesubmission.html',query_id=query_id,filename=filename)
	#except FileNotFoundError:
	return render_template('intronsubmission.html',query_id=query_id)
	
	
@app.route('/intronRes',methods=['GET','POST'])
def intronRes():
	query_id = request.args.get('query_id')
	return render_template('intronres.html',query_id=query_id)
		


		
	
@app.route('/spliceSub',methods=['GET','POST'])
def spliceSub():
	if request.method == 'POST':
		print('subbed')
		text=request.form['data_textbox']
		user=request.form['userid']
		uploaded_file=request.files['data_file']
		#if os.path.exists('regsnps-splicing/query_jobs/')==False:
		#	subprocess.run(['mkdir -p '+target_dir],shell=True)
			#subprocess.run(['mkdir(regsnps-splicing/query_jobs/)'])
			
		curr_time=datetime.datetime.now().strftime('%Y-%m-%d-%H-%M-%S')
		target_dir="./allJobs/Splicejob_"+curr_time
		#target_dir="regsnps-splicing/query_jobs/job_"+curr_time
  		
		subprocess.run(['mkdir -p '+target_dir],shell=True)
		#subprocess.run(['mkdir('+target_dir+')'])

		
		atts=request.form.getlist('attribute')
		#creating selected attribute file 
		with open('./allJobs/Splicejob_'+curr_time+'/checked_attributes','w') as f:
		#with open('regsnps-splicing/query_jobs/job_'+curr_time+'/checked_attributes','w') as f:
			for a in atts:
				f.write(a+'\n')
		#FILE is uploaded
		uploaded_fl=''
		input_desp=''
		real_result=''
		
		
		if uploaded_file.filename!='':
			uploaded_file.save('./allJobs/Splicejob_'+curr_time+'/'+uploaded_file.filename)
			uploaded_fl='./allJobs/Splicejob_'+curr_time+'/'+uploaded_file.filename
			real_result='./allJobs/Splicejob_'+curr_time+'/result.html'
			#uploaded_file.save('regsnps-splicing/query_jobs/job_'+curr_time+'/'+uploaded_file.filename)
			#uploaded_fl='./regsnps-splicing/query_jobs/job_'+curr_time+'/'+uploaded_file.filename
			#real_result='./regsnps-splicing/query_jobs/job_'+curr_time+'/result.html'
			input_desp=uploaded_file.filename
		else:
			uploaded_fl='./allJobs/Splicejob_'+curr_time+'/query'
			real_result='./allJobs/Splicejob_'+curr_time+'/result.html'
			#uploaded_fl='regsnps-splicing/query_jobs/job_'+curr_time+'/query'
			#real_result='./regsnps-splicing/query_jobs/job_'+curr_time+'/result.html'
			input_desp='query'
			with open(uploaded_fl,'w') as f:
				f.write(text)
			
		
		cmd=['setsid perl ./regsnps-splicing/bg_model.pl '+uploaded_fl+' '+input_desp+' '+curr_time+' '+real_result]
		status=subprocess.Popen(cmd,shell=True)
		
		
		##ADDING TO THE FILES LIST############
		currjson=[]
		if os.path.exists(files):
			with open(files) as f:
				currjson=json.load(f)
		else:
			with open(files,'w') as f:
				f.write('[]')
		
		addition=json.loads('{"user":"'+str(user)+'","prg":"Splice","date":"'+str(curr_time)+'"}')
		currjson.append(addition)
		with open(files,'w') as f:
			json.dump(currjson,f)
		os.system('sed "s/.*/var files=&;/" '+str(files)+' >  ./static/js/files.js')
		os.system('cat ./static/js/files.js ./static/js/index.js > ./static/js/filesAccess.js')
			
			
		#print(status.check_output())
		return redirect(url_for('spliceAftSub',query_id=curr_time,status=status.pid+1,filename=input_desp))

	return render_template('splicevariant.html')
	
@app.route('/spliceAftSub',methods=['GET','POST'])
def spliceAftSub():
	pid=request.args['status']
	query_id=request.args['query_id']
	filename=request.args['filename']
	print(pid)
	
	try: 
		os.kill(int(pid),0)
		return render_template('splicesubmission.html',query_id=query_id,filename=filename)
	except ProcessLookupError:
		return render_template('spliceres.html',query_id=query_id,filename=filename)
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
		print('hello')
		text=request.form['data_textbox']
		user=request.form['userid']
		uploaded_file=request.files['data_file']
		#if os.path.exists('regsnps-splicing/query_jobs/')==False:
		#	subprocess.run(['mkdir -p '+target_dir],shell=True)
			#subprocess.run(['mkdir(regsnps-splicing/query_jobs/)'])
			
		query_id=datetime.datetime.now().strftime('%Y-%m-%d-%H-%M-%S')
		target_dir="./allJobs/Exonjob_"+query_id
		#target_dir="./ExonImpact/usr_input/job_"+query_id
  		
		subprocess.run(['mkdir -p '+target_dir],shell=True)
		
		uploaded_fl=''
		input_desp=''
		real_result=''
		
		print('hello')
		if uploaded_file.filename!='':
			uploaded_file.save(target_dir+'/'+uploaded_file.filename)
			
			uploaded_fl=target_dir+'/'+uploaded_file.filename
			input_desp=uploaded_file.filename
		else:
			uploaded_fl=target_dir+'/query'
			input_desp='query'
			with open(uploaded_fl,'w') as f:
				f.write(text)
				
		#with open('input_desp.txt','w') as f:
				#f.write(input_desp)
		#print(input_desp)
		cmd=['~/Downloads/jdk1.8.0_271/bin/java -Xmx2g -cp "newRun.jar:./ExonImpact/lib_jar/*" ccbb.hrbeu.exonimpact.test.Exonimpact_for_server '+uploaded_fl+'; Rscript ./ExonImpact/src/R/predict.r '+uploaded_fl+'_features.csv > '+target_dir+'/error2_'+input_desp]
		#status=subprocess.Popen(cmd,shell=True)
		#cmd=['Rscript ./ExonImpact/src/R/predict.r '+uploaded_fl+'_features.csv > '+target_dir+'/error2_'+input_desp]
		#$error=exec("/usr/bin/Rscript R/predict.r usr_input/".$query_file_name. "_features.csv > error/error2_".$query_file_name." 2>&1",$output2);
		#$error=exec("mv ".$query_file_name." usr_input/");
		#cmd=[' setsid perl ./regsnps-splicing/bg_model.pl '+uploaded_fl+' '+input_desp+' '+curr_time+' '+real_result]
		status=subprocess.Popen(cmd,shell=True)
		
		
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
		
		#print(status.check_output())
		return redirect(url_for('exonAftSub',query_id=query_id,status=status.pid+1,filename=input_desp))
	return render_template('exonvariant.html')
		
@app.route('/exonAftSub',methods=['GET','POST'])
def exonAftSub():
	print('made it')
	#if request.method=='POST':
	pid=request.args['status']
	#print(pid)
	query_id=request.args['query_id']
	#print(query_id)
	file_name=request.args['filename']
	
	try: 
		os.kill(int(pid),0)
		return render_template('exonsubmission.html', query_id=query_id)
	except:
		return render_template('exonres.html',query_id=query_id,filename=file_name)
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
	file_name=request.args['filename']
	
	#try: 
	#	os.kill(int(pid),0)
	return render_template('exonres.html', query_id=query_id,file_name=file_name)

		



