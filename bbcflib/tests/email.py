# Built-in modules #
import socket, re, ConfigParser, cStringIO

# Internal modules #
from ..email import EmailReport

# Unitesting module #
try:
    import unittest2 as unittest
except ImportError:
    import unittest

# Nosetest flag #
__test__ = True

#--------------------------------------------------------------------------------#
def hostname_contains(pattern):
    hostname = socket.gethostbyaddr(socket.gethostname())[0]
    if re.search(pattern, hostname) == None:
        return False
    else:
        return True

#--------------------------------------------------------------------------------#
test_config_file = '''[emailreport]
email_sender=nobody@localhost
email_smtp_server=localhost
email_default_subject=Default Subject'''
def get_config_file_parser():
    file = cStringIO.StringIO()
    file.write(test_config_file)
    file.seek(0)
    config = ConfigParser.ConfigParser()
    config.readfp(file)
    return config

###################################################################################
class Test_EmailReport(unittest.TestCase):
    def setUp(self):
        self.report = EmailReport(sender='nobody@localhost', to='boris@localhost', subject='Default Subject', smtp_server='localhost')
        self.report_from_config = EmailReport(config=get_config_file_parser(), to='boris@localhost')

    def test_config_requires_to(self):
        self.assertRaises(TypeError,
                          lambda : EmailReport(config=get_config_file_parser()))


    def test_init(self):
        self.assertEqual(self.report.sender, 'nobody@localhost')
        self.assertEqual(self.report.to, 'boris@localhost')
        self.assertEqual(self.report.subject, 'Default Subject')
        self.assertEqual(self.report.smtp_server, 'localhost')

    def test_init_from_config(self):
        self.assertEqual(self.report_from_config.sender, 'nobody@localhost')
        self.assertEqual(self.report_from_config.to, 'boris@localhost')
        self.assertEqual(self.report_from_config.subject, 'Default Subject')
        self.assertEqual(self.report_from_config.smtp_server, 'localhost')

    def test_append_body(self):
        a = "This is a test of the emergency broadcasting system."
        self.report.appendBody(a)
        self.assertEqual(self.report.body, a)
        self.report.appendBody(a)
        self.assertEqual(self.report.body, a+a)

    def test_dump(self):
        self.assertEqual(self.report.dump(),
                         {'sender': 'nobody@localhost',
                          'to': 'boris@localhost',
                          'subject': 'Default Subject',
                          'smtp_server': 'localhost',
                          'body': self.report.body})

    def test_send(self):
        try:
            socket.create_connection(('localhost', 25))
        except:
            self.skipTest("Email send tests skipped because you can't send messages.")
        self.report_from_config.send()

#-----------------------------------#
# This code was written by the BBCF #
# http://bbcf.epfl.ch/              #
# webmaster.bbcf@epfl.ch            #
#-----------------------------------#
