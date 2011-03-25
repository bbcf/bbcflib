Tutorial
========

Here is a typical workflow that uses both ``mapseq`` and ``chipseq``. First import all the relevant modules::

    from bbcflib import daflims, genrep, frontend
    from bbcflib.mapseq import *
    from bbcflib.chipseq import *
    
Then connect to a ``MiniLIMS`` and define a set of jobs::

    M = MiniLIMS( 'test_lims' )
    fastq_path = "/archive/epfl/bbcf/data/fastq"
    job = {'groups': {
        'TestA': {'runs':{'1':'test_a_1.fastq', '2': 'test_a_2.fastq'},'control': False},
        'TestB': {'runs':{'1':'test_b_1.fastq'},'control': False},
	'Control': {'runs':{'1':'control_1.fastq','2':'control_2.fastq'},'control': True}
    },
    'assembly_id': "sacCer2",
    'options': {'discard_pcr_duplicates': True,
                'merge_strands': -1,
                'peak_calling': True,
                'ucsc_bigwig': False,}}
    fastq_files = dict((k,v['runs']) for k,v in job['groups'].iteritems())

Connect to ``GenRep`` to retrieve genome annotations::

    g_rep = genrep.GenRep( 'http://bbcftools.vital-it.ch/genrep/', '/scratch/frt/yearly/genrep/nr_assemblies/bowtie' )

Start a first execution to run ``bowtie`` on every fastq file, then create a pdf report of the mapping::

    with execution( M, description="Test mapping" ) as ex:
        processed = map_groups( ex, job, fastq_files, fastq_path, g_rep )
        pdf = add_pdf_stats( ex, p, dict((k,k) for k in job['groups'].keys()), '/archive/epfl/bbcf/share' )
    exid=ex.id
    print exid

From that first execution, import the resulting bamfiles and proceed to the Chipseq analysis::

    with execution( M, description="Test chipseq" ) as ex:
        (p,g)=import_mapseq_results( exid, M, ex.working_directory, job )
        g_rep_assembly = g_rep.assembly( job['assembly_id'] )
        (p,g) = workflow_groups( ex, job, p, g_rep_assembly.chromosomes, '/archive/epfl/bbcf/share' )

Finally get a dictionary of all the files produced by the analysis::

    print ex.id
    allfiles = get_files( ex.id, M )

