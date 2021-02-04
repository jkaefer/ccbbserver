#!/usr/bin/env python
import argparse
import smtplib
from email.mime.multipart import MIMEMultipart
from email.mime.text import MIMEText

class EmailServer(object):
    def __init__(self, host, port, email_address, password):
        self.email_address = email_address
        self.server = self._start_server(host, port, email_address, password)

    def _start_server(self, host, port, email_address, password):
        server = smtplib.SMTP(host, port)
        server.starttls()
        server.login(email_address, password)
        return server

    def send_msg(self, to_address, msg):
        self.server.sendmail(self.email_address, to_address, msg)

class RegSNPEmail(object):
    def __init__(self, host, port, email_address, password):
        self.server = EmailServer(host, port, email_address, password)

    def send_start(self, to_address, description='', link=''):
        msg = MIMEMultipart('alternative')
        msg['Subject'] = '[regSNP-intron] Job Submitted'
        msg['From'] = self.server.email_address
        msg['To'] = to_address
        
        text = 'Dear regSNP-intron user,\nYour job has been submitted. The results will be sent to you once the job is completed.\nYou can find the current status of your job here: {0}.\nregSNP-intron team'.format(link)
        html = '''\
        <html>
          <head></head>
          <body>
            <p>Dear regSNP-intron user,</p>
            <p>Your job has been submitted. The results will be sent to you once the job is completed.</p>
            <p>You can check the current status of your job <a href="{}">here</a>.</p>
            <p>regSNP-intron team</p>
          </body>
        </html>
        '''.format(link)
        
        part1 = MIMEText(text, 'plain')
        part2 = MIMEText(html, 'html')
        
        msg.attach(part1)
        msg.attach(part2)
        self.server.send_msg(to_address, msg.as_string())

    def send_end(self, to_address, description='', link=''):
        msg = MIMEMultipart('alternative')
        msg['Subject'] = '[regSNP-intron] Job Completed'
        msg['From'] = self.server.email_address
        msg['To'] = to_address
        
        text = 'Dear regSNP-intron user,\nYour job has been completed.\nYou can find your results here: {0}.\nregSNP-intron team'.format(link)
        html = '''\
        <html>
          <head></head>
          <body>
            <p>Dear regSNP-intron user,</p>
            <p>Your job has been completed.</p>
            <p>You can find your results <a href="{0}">here</a>.</p>
            <p>regSNP-intron team</p>
          </body>
        </html>
        '''.format(link)
        
        part1 = MIMEText(text, 'plain')
        part2 = MIMEText(html, 'html')
        
        msg.attach(part1)
        msg.attach(part2)
        self.server.send_msg(to_address, msg.as_string())

def main():
    parser = argparse.ArgumentParser(description='''Sending email notification
            for regsnp-intron webserver.''')
    parser.add_argument('-s', '--server', default='smtp.gmail.com',
            help='email server')
    parser.add_argument('-p', '--port', type=int, default='587',
            help='port')
    parser.add_argument('-a', '--account', default='regsnpintron@gmail.com',
            help='account')
    parser.add_argument('-pw', '--password', default='ccbbliulab',
            help='password')
    parser.add_argument('status',
            help='start|end')
    parser.add_argument('query_id',
            help='query_id')
    parser.add_argument('email',
            help='email')
    args = parser.parse_args()
    server = RegSNPEmail(args.server, args.port, args.account, args.password)
    if args.status == 'start':
        server.send_start(args.email, link='https://regsnps-intron.ccbb.iupui.edu/submission.php?query_id={0}'.format(args.query_id))
    elif args.status == 'end':
        server.send_end(args.email, link='https://regsnps-intron.ccbb.iupui.edu/result.php?query_id={0}'.format(args.query_id))

if __name__ == '__main__':
    main()





