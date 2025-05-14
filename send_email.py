import smtplib
from email.mime.multipart import MIMEMultipart
from email.mime.text import MIMEText
from email.mime.image import MIMEImage
import os

def send_email(index):
    #创建邮件对象
    msg=MIMEMultipart()
    msg['Subject']=index['subject']
    msg['From']=index['sender']
    msg["To"] = ",".join(index['recipient'])
    #添加邮件正文
    msg.attach(MIMEText(index['message'],'plain'))
    for res in os.listdir(index['attachment_path']):
            with open(index['attachment_path']+"/"+res,'rb') as f:
                part = MIMEImage(f.read())
                part.add_header('application','attachment', filename = res)
                msg.attach(part)
    
    #连接邮件smtp服务器发送邮件
    with smtplib.SMTP_SSL('smtp.163.com', 465) as smtp:
        smtp.login(index['sender'], index['sender_code'])
        smtp.send_message(msg)
        print("邮件发送完毕")