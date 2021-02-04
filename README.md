# ccbbserver
Makes use of regSNPs-Intron at <https://github.com/linhai86/regsnp_intron>, ExonImpact at <https://github.com/regSNPs/ExonImpact> and algorithms designed by the regSNPs-Splicing project locatable at <https://pubmed.ncbi.nlm.nih.gov/28391525/>

## Stack Diagram
![alt text](https://github.com/jkaefer/ccbbserver/blob/master/diagram.png)


## Installation
1. An Anaconda environment is needed with python 3.7 or greater and the required packages are contained within the file py37list.txt, these can be batch installed by  
```bash
conda activate py37
conda install py37list.txt
```

## Other requirements
**Perl:**
**Java version 1.8.0_271:**
**Redis Queue:**
Follow the instructions at <https://redis.io/topics/quickstart> to install. Then launch an instance in a terminal with:
```bash
redis-server
```
**Celery:**
Automatically installed with the conda enviroment. To start celery daemon for enter ccbbserver directory and run:
```bash
celery -A app.routes:celery worker --loglevel=info
```


## Usage
From ccbbserver directory ensure the flask app environment variable is set to location of ccbbserver.py:
```bash
export FLASK_APP=~/Documents/ccbbserver/ccbbserver.py
flask run 
```
visit the hosted url

## Output
Query prediction results will be stored in the allJobs directory within the ccbbserver folder







