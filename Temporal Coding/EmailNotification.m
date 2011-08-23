function EmailNotification(myaddress,mypassword,toaddress,subject,message)
% function EmailNotification(myaddress,mypassword,subject,message)
% sends an email notification (through GMail)

% myaddress = 'myaddress@gmail.com';
% mypassword = 'mypassword';

setpref('Internet','E_mail',myaddress);
setpref('Internet','SMTP_Server','smtp.gmail.com');
setpref('Internet','SMTP_Username',myaddress);
setpref('Internet','SMTP_Password',mypassword);

props = java.lang.System.getProperties;
props.setProperty('mail.smtp.auth','true');
props.setProperty('mail.smtp.socketFactory.class', ...
                  'javax.net.ssl.SSLSocketFactory');
props.setProperty('mail.smtp.socketFactory.port','465');

sendmail(toaddress, subject, message);
