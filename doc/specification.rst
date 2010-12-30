``bbcflib`` Specification

BBCFLib is a set of Python modules to make it easy to deploy bein workflows for general users.  It includes modules to do lookups from GenRep, fetch files from the DAF LIMS, interface with Fabrice David's Ruby on Rails web frontends, send report emails, and load configuration files for the workflow.

This specification is by no means complete or final.  Everything will doubtless get revised several times along the way, but it gives us something to aim for.

*********
Scenarios
*********

Getting a configuration
-------------------------

Jacques is writing a workflow.  He needs configuration fields "boris" and "natasha" in a section "myfields".  He creates a ConfigParser object, declaring that "natasha" defaults to "".::

    cm = Configuration()
    cm.add_section("myfields")
    cm.set('myfields','natasha','')

Then when he loads a configuration file "myconfig.txt", and reads the "boris" field from the Configuration.::

    cm.read("myconfig.txt")
    cm.get("myfields", "boris")

He tries to fetch the field "natasha", but misspells it.::

    cm.get("myfields", "natahsa")

The program throws a NoOptionError.

Fetching the path to a bowtie index from GenRep
-----------------------------------------------

Marion needs a bowtie index for her workflow.  She creates a GenRep object pointing to ``http://bbcftools.vital-it.ch/genrep/`` with::

    g = GenRep("http://bbcftools.vital-it.ch/genrep/", "/path/to/indices")
    a = g.get_assembly("MSmeg_JS623")

Fetch a FASTQ file from the DAF LIMS
------------------------------------

Jacques needs to fetch a file with reads from the DAF LIMS.  He has the facility, machine, run, and lane as "Lausanne", "R2D2", "942", and "6".  He creates a DAFLIMS object with the right URL, username, and password, and then calls fetch on it to get the file.::

    df = DAFLIMS(username="jacques", password="mypassword")
    df.fetch_file(facility="Lausanne", machine="R2D2", run=942, lane=6, write_to="/path/to/file.fastq")

This writes the file to ``/path/to/file.fastq``.

Send an email report
--------------------

Jacques needs to send an email to a user when a job finishes running.  He creates an EmailReport object::

    er = EmailReport(from="harry@nowhere.com", to="boris@nowhere.com",
                     subject="Your job finished!",
                     smtp_server="smtp.nowhere.com")
    er.appendToBody("Here's your job output.")
    er.send()

This sends off the message.

Configuring from a state
------------------------

Mado is sending email reports, fetching from the DAF LIMS and GenRep, and needs to configure these from a file, not with fields in the code itself.  In all the object creations above, she hands a single object, a Configuration, to the constructor.::

    g = GenRep(cm)
    df = DAFLIMS(cm)
    er = EmailReport(cm, to="boris@nowhere.com")

Now the fields like subject and smtp_server in EmailReport or username and password are filled in with fields from the ConfigParser.  Each defaults to a particular section name, but she could also set the section name with section=...

Fetching info from Fabrice's frontends
--------------------------------------

Marion's workflow is run from Fabrice's frontend with a run ID of 36.  Her workflow receives this as a command line argument, so she needs to fetch the information about the run.  She creates a Frontend object with that ID, which queries his frontend for metadata, groups, and samples.::

    f = Frontend(url="http://3c.vital-it.ch", run=36)
    Frontend.email # The email address to send reports to
    Frontend.groups # A dictionary of groups pointing to samples

She could also have written Frontend(cm, run=36) instead to load the URL from the Configuration.

********
Nongoals
********

*******
Details
*******

Configuration
-------------

bbcflib uses ConfigParser to handle configuration files.  The other components all allow you to pass a ConfigParser and an optional section argument in place of other constructor arguments, which will then pull the right fields from the ConfigParser.  A configuration file template for all these sections is::

    [genrep]
    genrep_url=
    genrep_root=

    [daflims]:
    daflims_username=
    daflims_password=

    [emailreport]:
    email_from_address=
    email_smtp_server=
    email_default_subject=

    [frontend]:
    frontend_url=

To add a configuration file, call read("filename") on the ConfigParser.  "filename" is the name of the INI style file to be read.

To get values from the Configuration, call the methods

get("section", "fieldname")
getboolean("section", "fieldname")
getint("section", "fieldname")
getfloat("section", "fieldname")

If no value is set, a NoOptionError is thrown.

GenRep
------

GenRep takes two arguments: a path to the GenRep JSON methods, and a path to its files.  These may both be specified by a ConfigParser and optional section (defaulting to "genrep") instead.  So the constructor is::

    def __init__(self, url=None, path=None, config=None, section='genrep'):
        # check if we have url+path or config+section.  If neither, throw TypeError

The user generally only calls get_assembly on GenRep.  It takes either a string or an integer, specifying an assembly name or ID, and returns and Assembly object.

DAF LIMS
--------

DAFLIMS takes either a username and password, or a ConfigParser and optional section argument.  Since the DAF LIMS is only usable from VITAL-IT, we hard code all the details into the object.  Then there is one method for use by the user, fetch_file, which takes facility, machine, run, lane, and an optional write_to.  If write_to is omitted, the original filename is written in the current working directory.

Email Report
------------

A single EmailReport object, either taking from, to, subject, smtp_server, or config, section (optional), and to.  Afterwards, call appendToBody to add lines to the body, and send() to send the message off.

Frontend
--------

A Frontend object takes either a url and run id or a ConfigParser and run id.  The __init__ method fetches all the information from there, and sets local fields:

email (string)
groups ({"group name": {"control": bool, "samples": {sample_id: {"facility":, "machine":, "run":, "lane": }}}})
description (string)

For testing, use the key: 9pv1x7PamOj80eXnZa14

When running your own bein execution for the frontend, set the execution's description to the key so the frontend can find it afterwards.