Tutorial
========

Here is a typical workflow that uses both ``mapseq`` and ``chipseq``. First prepare a configuration file with the following sections: first a ``Global variables`` section that defines the local environment for the pipelines::

    [Global variables]
    genrep_url='http://bbcftools.vital-it.ch/genrep/'
    bwt_root='/db/genrep'
    fastq_root='/scratch/cluster/daily/htsstation/mapseq/'
    script_path='/archive/epfl/bbcf/share'
    [[hts_chipseq]]
    url='http://htsstation.vital-it.ch/chipseq/'
    download='http://htsstation.vital-it.ch/lims/chipseq/chipseq_minilims.files/'
    [[hts_mapseq]]
    url='http://htsstation.vital-it.ch/mapseq/'
    download='http://htsstation.vital-it.ch/lims/mapseq/mapseq_minilims.files/'
    [[gdv]] 
    url='http://svitsrv25.epfl.ch/gdv'
    email='your.email@yourplace.org'
    key='pErS0na1&keY%0Ng2V'

For example, if you intend to download data from the LIMS, you need to setup your account::

    [[lims]] 
    user='limslogin'
    [[[passwd]]] 
    lgtf='xxxxxxxx'
    gva='yyyyyyy'

Similarly, if you want to receive an email upon completion of the pipeline, submit the relevant informations for the sender email::

    [[email]] 
    sender='webmaster@lab.edu'
    smtp='local.machine.edu'
    
Then comes the job description::

    [Job]
    description='config test'
    assembly_id=mm9
    email=toto@place.no
    [Options]
    read_extension=65
    input_type=0
    compute_densities=True
    discard_pcr_duplicates=True
    
Experimental conditions correspond to `groups` which are numbered, each condition may have any number of replicates (called `runs`) which are associated with their respective group via its numeric `group_id`::

    [Groups]
    [[1]]
    control=True
    name=unstimulated
    [[2]]
    control=False
    name=stimulated
    
    [Runs]
    [[1]]
    url=http://some.place.edu/my_control.fastq
    group_id=1
    [[2]]
    url=http://some.place.edu/my_test1.fastq
    group_id=2
    [[3]]
    url=http://some.place.edu/my_test2.fastq
    group_id=2

Such a configuration file can be passed as command-line argument to the scripts `run_mapseq.py <https://github.com/bbcf/bbcfutils/blob/master/Python/run_mapseq.py>`_ and `run_chipseq.py <https://github.com/bbcf/bbcfutils/blob/master/Python/run_chipseq.py>`_, e.g.::

    python run_mapseq.py -c my_config.txt -d test_lims

We next analyse how these python scripts are using these configuration and processing the files. First we import all the relevant modules::

    from bbcflib import daflims, genrep, frontend, gdv, common
    from bbcflib.mapseq import *
    from bbcflib.chipseq import *
    
Then connect to a ``MiniLIMS`` and parse the configuration file::

    M = MiniLIMS( 'test_lims' )
    (job,gl) = frontend.parseConfig( 'my_config.txt' )

This returns two dictionaries, one with the job description and one with the global variables sections. Then we fetch an assembly and define a few options::

    g_rep = genrep.GenRep( url=gl["genrep_url"], root=gl["bwt_root"], intype=job.options.get('input_type') )
    assembly = g_rep.assembly( job.assembly_id )
    dafl = dict((loc,daflims.DAFLIMS( username=gl['lims']['user'], password=pwd )) for loc,pwd in gl['lims']['passwd'].iteritems())
    job.options['ucsc_bigwig'] = job.options.get('ucsc_bigwig') or True
    job.options['gdv_project'] = job.options.get('gdv_project') or False
    via = 'lsf'

then start an execution environment in which we

* fetch the fastq files using :func:`bbcflib.mapseq.get_fastq_files`
* launch the bowtie mapping via :func:`bbcflib.mapseq.map_groups`
* generate a pdf report of the mapping statistics with :func:`bbcflib.mapseq.add_pdf_stats`
* if requested, make a density profile using :func:`bbcflib.mapseq.densities_groups`
* create the corresponfing project and tracks in :doc:`GDV <bbcflib_gdv>`.

This corresponds to the code below::

    with execution( M, description='test_mapseq' ) as ex:
        job = get_fastq_files( job, ex.working_directory, dafl )
        mapped_files = map_groups( ex, job, ex.working_directory, assembly, {'via': via} )
        pdf = add_pdf_stats( ex, mapped_files,
                             dict((k,v['name']) for k,v in job.groups.iteritems()),
                             gl['script_path'] )
        if job.options['compute_densities']:
            if not(job.options.get('read_extension')>0):
                job.options['read_extension'] = mapped_files.values()[0].values()[0]['stats']['read_length']
            density_files = densities_groups( ex, job, mapped_files, assembly.chromosomes, via=via )
            if job.options['gdv_project']:
                gdv_project = gdv.create_gdv_project( gl['gdv']['key'], gl['gdv']['email'],
                                                      job.description, hts_key, 
                                                      assembly.nr_assembly_id,
                                                      gdv_url=gl['gdv']['url'], public=True )
                add_pickle( ex, gdv_project, description='py:gdv_json' )

Finally all the output files are returned as a dictionary::

    allfiles = common.get_files( ex.id, M )

this dictionary will be organized by file type and provide a descriptive name and the actual (repository) file name, e.g.::

    {'none': {'7XgDex9cTCn8JjEk005Q': 'test.sql'}, 
    'py': {'hkwjU7nnhE0uuZostJmF': 'file_names', 'M844kgtaGpgybnq5APsb': 'test_full_bamstat', 'cRzKabyKnN0dcRHaAVsj': 'test_Poisson_threshold', 'j4EWGj2riic7Xz47hKhj': 'test_filter_bamstat'}, 
    'sql': {'7XgDex9cTCn8JjEk005Q_merged.sql': 'test_merged.sql'}, 
    'bigwig': {'UjaseL2p8Z1RnDetZ2YX': 'test_merged.bw'},
    'pdf': {'13wUAjrQEikA5hXEgTt': 'mapping_report.pdf'}, 
    'bam': {'mJP4dqP1f2K6Pw2iZ2LZ': 'test_filtered.bam', 'IRn3o49zIZ2JOOkMxAJl.bai': 'test_complete.bam.bai', 'IRn3o49zIZ2JOOkMxAJl': 'test_complete.bam', 'mJP4dqP1f2K6Pw2iZ2LZ.bai': 'test_filtered.bam.bai'}}

If you then want to continue with a ChIP-seq analysis, you can start a new execution, collect the files with :func:`bbcflib.chipseq.get_bam_wig_files` and run :func:`bbcflib.chipseq.workflow_groups`::

    with execution( M, description='test_chipseq' ) as ex:
        (mapped_files, job) = get_bam_wig_files( ex, job, 'test_lims', gl['hts_mapseq']['url'], gl['script_path'], via=via )
        chipseq_files = workflow_groups( ex, job, mapped_files, assembly.chromosomes, gl['script_path'] )


