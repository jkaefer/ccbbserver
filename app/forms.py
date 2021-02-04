from flask_wtf import FlaskForm,Form
from flask_wtf.file import FileField
from wtforms import StringField, PasswordField, BooleanField, SubmitField,SelectField,validators,SelectMultipleField
from wtforms.validators import DataRequired, ValidationError,Required,Email,EqualTo
#from app import session
#from app.models import User
#import MySQLdb
#from app import app

#db=MySQLdb.connect("localhost","root","12345","db20")




class ExonSubmitForm(Form):
	manualData=StringField('input',validators=[DataRequired()])
	#file = FileField('Browse...')
	inputFormat=SelectField(u'Format', choices = ['txt','vcf'], validators = [Required()])
	genomeAssembly=SelectField(u'Assembly', choices = ['CRCh37/hg19'], validators = [Required()])
	
	email=StringField('Email',validators=[DataRequired(),Email()])
	submit=SubmitField('Submit')
	


class IntronSubmitForm(Form):
	manualData=StringField('input',validators=[DataRequired()])
	file = FileField('Browse...')
	inputFormat=SelectField(u'Format', choices = [('txt','txt'),('vcf','vcf')], validators = [Required()])
	genomeAssembly=SelectField(u'Assembly', choices = [('CRCh37/hg19','CRCh37/hg19')], validators = [Required()])
	
	email=StringField('Email',validators=[DataRequired(),Email()])
	submit=SubmitField('Submit')
	
class SpliceSubmitForm(Form):
	manualData=StringField('input',validators=[DataRequired()])
	file = FileField('Browse...')
	inputFormat=SelectField(u'inputFormat', choices = [('txt','txt'),('vcf','vcf')], validators = [Required()])
	genomeAssembly=SelectField(u'inputFormat', choices = [('txt','txt'),('vcf','vcf')], validators = [Required()])
	
	email=StringField('Email',validators=[DataRequired(),Email()])
	submit=SubmitField('Submit')
	
	



