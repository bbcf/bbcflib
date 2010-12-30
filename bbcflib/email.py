# import smtplib
# from email.mime.text import MIMEText



# msg = MIMEText("This is a test")

# msg['Subject'] = 'Testing!'
# msg['From'] = 'ross@localhost'
# msg['Reply-To'] = 'ross@localhost'
# msg['To'] = 'ross'

# s = smtplib.SMTP()
# s.connect("localhost")
# s.sendmail('ross@localhost','ross@localhost',msg.as_string())
# s.close()


from mailer import Mailer
from mailer import Message

message = Message(From="ross@localhost",
                  To="ross@localhost",
                  charset="utf-8")
message.Subject = "Boris"
message.Body = "Boris the mad babboon runs amok!\n"

sender = Mailer('localhost')
sender.send(message)
