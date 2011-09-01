"""
=======================
Module: bbcflib.c4seq
=======================

This module provides functions to run a 4c-seq analysis from reads mapped on
a reference genome.

"""

from bein import *
from bein.util import *
from bbcflib import daflims, genrep, frontend, email, gdv, common, track
from bbcflib.mapseq import *
import sys, getopt, os, json, re
import gMiner as gm
import createlib


#-------------------------------------------#
# Functions 
#-------------------------------------------#
# call script segToFrag.awk
@program
def segToFrag(in_countsPerFragFile,regToExclude=None, script_path='./'):
	def parseStdout(p):
		print('will parse the Stdout')
		filename = unique_filename_in()
		output=open(filename,'w')
		s=''.join(p.stdout)
     	        output.write(s)
		return filename

	if regToExclude == None:
		print('will call: awk -f '+script_path+'segToFrag.awk '+in_countsPerFragFile)
		return {'arguments': ["awk","-f",script_path+'segToFrag.awk',in_countsPerFragFile],
			'return_value':None}
	else:	
		print('will call: awk -f '+script_path+'segToFrag.awk '+' -v regToExclude='+regToExclude+' '+in_countsPerFragFile)
		return {'arguments': ["awk","-f",script_path+'segToFrag.awk ',"-v","regToExclude="+regToExclude,in_countsPerFragFile],
			'return_value':None}

def call_segToFrag(*args, **kwargs):
	filename = unique_filename_in()
        kwargs["stdout"] = filename
 	future=segToFrag.nonblocking(*args, **kwargs)
	future.wait()
        return filename

# *** parse the output of call_segToFrag
def parseSegToFrag(infile):
	filename = unique_filename_in()
	output = open(filename,'w')
	with open(infile,"r") as f:
		for s in f.readlines():
			s=s.strip('\n')
			if re.search(r'IsValid',s.split('\t')[2]) and float(s.split('\t')[11])>0.0 :
				coord=((s.split('\t')[1]).split(':')[1]).split('-')
				output.write((s.split('\t')[1]).split(':')[0]+'\t'+str(int(coord[0])-1)+'\t'+coord[1]+'\t'+s.split('\t')[11]+'\n')
	output.close()
	return filename

# *** To sort a file on coordinates
@program
def sortOnCoord(infile):
	return{'arguments': ["sort","-k1,1","-k2,2n","-k3,3n",infile],
		'return_value':None }

def call_sortOnCoord(*args, **kwargs):
	filename = unique_filename_in()
        kwargs["stdout"] = filename
	future=sortOnCoord.nonblocking(*args, **kwargs)
	future.wait()
	return filename


# *** Create a dictionary with infos for each primer (from file primers.fa)
# ex: primers_dict=loadPrimers('/archive/epfl/bbcf/data/DubouleDaan/finalAnalysis/XmNGdlXjqoj6BN8Rj2Tl/primers.fa')
def loadPrimers(primersFile):
	primers={}
	with open(primersFile,'r') as f:
		for s in f.readlines():
			s=s.strip('\n')
			if re.search(r'^>',s):
				primerInfos={}
				infos=s.split('|')
				n=len(infos)-1
				name=infos[0][1:len(infos[0])]
				primerInfos['fullseq']=infos[1]
				primerInfos['baitcoord']=infos[2]
				primerInfos['primary']=infos[3]
				if re.search('Exclude',infos[len(infos)-1]):
					n=n-1
					primerInfos['regToExclude']=(infos[len(infos)-1]).split('=')[1]
				primerInfos['seqToFilter']=infos[4:n]
				if not name in primers:	
					primers[name]=primerInfos
					prevPrimer=name
			else:
				primers[name]['seq']=s
	return primers


# *** main function to compute normalised counts per fragments from a density file
# ex: resfiles=density_to_countsPerFrag(ex,mapped_files[gid][rid]['wig']['merged'],mapped_files[gid][rid]['libname'],assembly.name,reffile,regToExclude,working_dir, script_path, 'lsf')
def density_to_countsPerFrag(ex,density_file,density_name,assembly_name,reffile,regToExclude,wd,script_path, via='lsf'):
	print("will call mean_score_by_feature for t1="+density_file+"(name="+density_name+") and t2="+reffile)
        outdir=unique_filename_in()
        #ok: res=gm.run(track1=density_file,track1_name=density_name,track2=reffile,track2_name='libFile',track2_chrfile=assembly.name,operation_type='genomic_manip',manipulation='mean_score_by_feature',output_location=wd,output_name=outdir)
	
	gMiner_job = { 'track1': density_file,
                                       'track1_name':density_name,
                                       'track2':reffile,
                                       'track2_name':'libFile',
                                       'track2_chrfile':assembly_name,
                                       'operation_type':'genomic_manip',
                                       'manipulation':'mean_score_by_feature',
                                       'output_location':wd,
                                       'output_name':outdir
                      }
	print(gMiner_job)
	res = common.run_gMiner.nonblocking(ex,gMiner_job,via='lsf').wait()
	ex.add(wd+outdir+".sql",description="sql:meanScorePerFeature_"+density_name)
	countsPerFragFile=unique_filename_in()+".bed"
	with track.load(wd+outdir+".sql",'sql') as t:
		t.convert(countsPerFragFile,'bed')
	ex.add(countsPerFragFile,description="bed:meanScorePerFeature_"+density_name)
	res=call_segToFrag(ex, countsPerFragFile, regToExclude, script_path, via=via)
	ex.add(res,description="none:res_segToFrag_"+density_name+" (rough)")
	resBedGraph=parseSegToFrag(res)
	ex.add(resBedGraph,description="none:res_segToFrag_"+density_name+" (bedGraph non-sorted)")
	resBedGraph=call_sortOnCoord(ex,resBedGraph,via=via)
	headerFile=unique_filename_in();
	hfile=open(headerFile,'w')
	hfile.write('track type="bedGraph" name="'+density_name+' normalised counts per valid fragments" description="'+density_name+' normalised counts per valid fragments" visibility=full windowingFunction=maximum autoScale=off viewLimits=1:2000\n')
	hfile.close()
	sortedBedGraph=common.cat([headerFile,resBedGraph])
	ex.add(sortedBedGraph,description="bedgraph:res_segToFrag_"+density_name+" (bedGraph sorted)")	
	sortedBedGraph_sql=unique_filename_in()+".sql"
	with track.load(sortedBedGraph,'bedGraph', chrmeta=assembly_name) as t:
                t.convert(sortedBedGraph_sql,'sql')
	ex.add(sortedBedGraph_sql,description="sql:res_segToFrag_"+density_name+" (bedGraph sorted)")
	return [wd+outdir+".sql",countsPerFragFile,res,resBedGraph,sortedBedGraph,sortedBedGraph_sql]

# Main 
#-------------------------------------------#
# *** open the 4C-seq minilims and create execution
# *** 0.get/create the library 
# *** 1.when necessary, calculate the density file from the bam file (mapseq.parallel_density_sql)
# ### 2.calculate the count per fragment for each denstiy file with gFeatMiner:mean_score_by_feature to calculate)
def workflow_groups(ex, job, primers_dict, g_rep, mapseq_files, mapseq_url, script_path='', via='lsf' ):
	assembly = g_rep.assembly(job.assembly_id)
	processed={
		'lib' : {},
		'density' : {},
		'4cseq' : {}
		}
	job_groups=job.groups
	htss_mapseq = frontend.Frontend( url=mapseq_url )

	new_libs=[]
        fasta_allchr=assembly
        #fasta_allchr='/archive/epfl/bbcf/mleleu/pipeline_vMarion/pipeline_3Cseq/vWebServer_Bein/tests/test.minilims.files/8LkjoeWBoRw0mh1rLz8V' #for test: will be /scratch/cluster/monthly/htsstation/4cseq/job.id

	for gid, group in job_groups.iteritems():
		reffile=createlib.get_libForGrp(ex,group,fasta_allchr,new_libs, job.id, g_rep)
#		reffile='/archive/epfl/bbcf/data/DubouleDaan/library_Nla_30bps/library_Nla_30bps_segmentInfos.bed'
		processed['lib'][gid]=reffile
		for rid,run in group['runs'].iteritems():
                        job_mapseq=htss_mapseq.job(run['key'])
			if 'regToExclude' in primers_dict[mapseq_files[gid][rid]['libname']]:
                                regToExclude=primers_dict[mapseq_files[gid][rid]['libname']]['regToExclude']
			else:
			        regToExclude=None
                        if not job_mapseq.options.get('compute_densities') or job_mapseq.options.get('merge_strands') != 0:
				print("will call parallel_density_sql with bam:"+mapseq_files[gid][rid]['bam']+"\n")
				density_file=mapseq.parallel_density_sql( ex, mapseq_files[gid][rid]['bam'],
                        						assembly.chromosomes,
                        			                        nreads=mapseq_files[gid][rid]['stats']["total"],
                                                 			merge=0,
                                                 			convert=False,
                                                 			via=via )
				mapseq_files[gid][rid]['wig']['merged']=density_file+"merged.sql"
				print("density file:"+mapseq_files[gid][rid]['wig']['merged'])
                        else:
                                print("Will use existing density file:"+mapseq_files[gid][rid]['wig']['merged'])

			ex.add(mapseq_files[gid][rid]['wig']['merged'],description='sql:density_file_'+mapseq_files[gid][rid]['libname'])
			processed['density'][mapseq_files[gid][rid]['libname']]=mapseq_files[gid][rid]['wig']['merged']

			print("Will process to the main part of 4cseq module: calculate normalised counts per fragments from density file:"+mapseq_files[gid][rid]['wig']['merged'])
			resfiles=density_to_countsPerFrag(ex,mapseq_files[gid][rid]['wig']['merged'],mapseq_files[gid][rid]['libname'],assembly.name,reffile,regToExclude,ex.remote_working_directory+'/',script_path, via)
			processed['4cseq']=resfiles
	return processed

