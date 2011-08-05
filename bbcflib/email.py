"""
=====================
Module: bbcflib.email
=====================

This module provides an simple way to set up and send emails, possibly
configuring the subject, sender, and SMTP server via an external
configuration file.  The standard usage is to create an
``EmailReport`` object::

    r = EmailReport(sender='nobody@localhost',
                    to='boris@localhost',
                    subject='Default Subject',
                    smtp_server='localhost')

then append lines to its body::

    r.appendBody("This is some text.")

and finally call the ``send`` method to dispatch the email.::

    r.send()

The ``EmailReport`` object can also be created referring to a
``ConfigParser``, though the *to* field must still be specified.  With
a ``ConfigParser`` object ``cp``, we would create an ``EmailReport``
with::

    r = EmailReport(config=cp, to='boris@localhost')

.. autoclass:: EmailReport
"""

# Other modules #
from mailer import Mailer, Message

################################################################################
class EmailReport(object):
    """Create an email, possibly reading a configuration file.

    An ``EmailReport`` can be created either by specifying the keyword
    arguments *sender*, *to*, *subject*, and *smtp_server*, or a
    ``ConfigParser`` object as the keyword argument *config* and the
    keyword argument *to*.  The section of the ``ConfigParser`` to use
    can be set with the keyword argument *section* (and defaults to
    "emailreport").

    The keyword arguments may also be specified along with a *config*
    argument, in which case they override the value in the
    configuration file.

    The fields in the configuration file are

      * email_sender
      * email_smtp_server
      * email_default_subject

    """
    def __init__(self, sender=None, to=None, subject=None, smtp_server=None, config=None, section='emailreport'):
        if (sender==None or to==None or subject==None or smtp_server==None) and config==None:
            raise TypeError("EmailReport requires from, to, subject, and smtp_server, or config.")
        elif to==None:
            raise TypeError("Must provide a 'to' address.")
        elif config != None:
            if sender == None: self.sender = config.get(section, 'email_sender')
            else:              self.sender = sender

            self.to = to

            if subject == None: self.subject = config.get(section, 'email_default_subject')
            else:               self.subject = subject

            if smtp_server == None: self.smtp_server = config.get(section, 'email_smtp_server')
            else:                   self.smtp_server = smtp_server
        else:
            self.sender = sender
            self.to = to
            self.subject = subject
            self.smtp_server = smtp_server

        self.body = ""


    def appendBody(self, s):
        """Append string *s* to the message body."""
        self.body = self.body + s

    def dump(self):
        """Return a dictionary describing this email."""
        return {'sender': self.sender,
                'to': self.to,
                'subject': self.subject,
                'smtp_server': self.smtp_server,
                'body': self.body}

    def send(self):
        """Send the email described by this object."""
        message = Message(From=self.sender,
                          To=self.to,
                          charset="utf-8")
        message.Subject = self.subject
        message.Body = self.body
        mailer = Mailer(self.smtp_server)
        mailer.send(message)

#-----------------------------------#
# This code was written by the BBCF #
# http://bbcf.epfl.ch/              #
# webmaster.bbcf@epfl.ch            #
#-----------------------------------#
