"""
=======================
Module: bbcflib.c4seq
=======================

This module provides functions to run a 4c-seq analysis from reads mapped on
a reference genome.

"""

from bein import *
from bein.util import touch
from bbcflib import daflims, genrep, frontend, email, gdv, track, createlib
from bbcflib.mapseq import *
from bbcflib.common import unique_filename_in
import sys, getopt, os, json, re
import sqlite3
import gMiner as gm

grpId=1
step=0

#-------------------------------------------#
# Functions 
#-------------------------------------------#
# call script segToFrag.awk
@program
def segToFrag(in_countsPerFragFile,regToExclude=None, script_path='./'):
	''' 
		This function calls segToFrag.awk (which transforms the counts per segment to a normalised count per fragment) 
		Gives the region to exclude if any 
	'''
	if regToExclude == None:
		print('will call: awk -f '+script_path+'segToFrag.awk '+in_countsPerFragFile)
		return {'arguments': ["awk","-f",script_path+'segToFrag.awk',in_countsPerFragFile],
			'return_value':None}
	else:	
		print('will call: awk -f '+script_path+'segToFrag.awk '+' -v reg2Excl='+regToExclude+' '+in_countsPerFragFile)
		return {'arguments': ["awk","-f",script_path+'segToFrag.awk ',"-v","reg2Excl="+regToExclude,in_countsPerFragFile],
			'return_value':None}

def call_segToFrag(*args, **kwargs):
	filename = unique_filename_in()
        kwargs["stdout"] = filename
 	future=segToFrag.nonblocking(*args, **kwargs)
	future.wait()
        return filename

# *** parse the output of call_segToFrag
def parseSegToFrag(infile):
	''' Parse the output of segToFrag '''
	filename_all = unique_filename_in()
	out_all = open(filename_all,'w')
	filename = unique_filename_in()
	output = open(filename,'w')
	with open(infile,"r") as f:
		for s in f:
			s=s.strip('\n')
			if re.search(r'IsValid',s.split('\t')[2]) :
				coord=((s.split('\t')[1]).split(':')[1]).split('-')
				out_all.write((s.split('\t')[1]).split(':')[0]+'\t'+str(int(coord[0])-1)+'\t'+coord[1]+'\t'+s.split('\t')[11]+'\n')
				if float(s.split('\t')[11])>0.0:
					output.write((s.split('\t')[1]).split(':')[0]+'\t'+str(int(coord[0])-1)+'\t'+coord[1]+'\t'+s.split('\t')[11]+'\n')
	output.close()
	out_all.close()
	return [filename,filename_all]

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
	'''
		Create a dictionary with infos for each primer (from file primers.fa) 
	'''
	primers={}
	with open(primersFile,'rb') as f:
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
                                        print "primerInfos['regToExclude']: " + primerInfos['regToExclude']
				primerInfos['seqToFilter']=infos[4:n]
				if not name in primers:	
					primers[name]=primerInfos
					prevPrimer=name
			else:
				primers[name]['seq']=s
	return primers


@program
def profileCorrection(inputFile,baitCoord,name,outputFile,reportFile,script_path='./'):
	return{'arguments': ["R","--vanilla","--no-restore","--slave","-f",script_path+"profileCorrection.R","--args",inputFile,baitCoord,name,outputFile,reportFile],
                'return_value':None}

@program
def smoothFragFile(inputFile,nFragsPerWin,curName,outputFile,regToExclude=None,script_path='./'):
        if not(regToExclude): regToExclude=''
	return{'arguments': ["R","--vanilla","--no-restore","--slave","-f",script_path+"smoothData.R","--args",inputFile,nFragsPerWin,curName,outputFile,regToExclude],
        	'return_value':None}

# *** main function to compute normalised counts per fragments from a density file
# ex: resfiles=density_to_countsPerFrag(ex,mapped_files[gid][rid]['wig']['merged'],mapped_files[gid][rid]['libname'],assembly,reffile,regToExclude,working_dir, script_path, 'lsf')
def density_to_countsPerFrag(ex,density_file,density_name,assembly,reffile,regToExclude,wd,script_path, via='lsf'):
	'''
		main function to compute normalised counts per fragments from a density file 
	'''
	global grpId
	global step

	print("will call mean_score_by_feature for t1="+density_file+"(name="+density_name+") and t2="+reffile)
        output = unique_filename_in()
        from gMiner.operations.genomic_manip.scores import mean_score_by_feature
        with track.Track(density_file) as scores:
                with track.Track(reffile) as features:
                        with track.new(output,format='sql',chrmeta=assembly.chrmeta) as out:
                                for ch in scores:
                                        out.write(ch,mean_score_by_feature()(
                                                        scores.read(ch),
                                                        features.read(ch,fields=['start', 'end', 'name'])),
                                                  fields=['start', 'end', 'name', 'score'])
        connection = sqlite3.connect(output)
        cursor = self.connection.cursor()
        cursor.execute("UPDATE 'attributes' SET 'value'='%s' WHERE 'key'='%s'"%('quantitative','datatype'))
        connection.commit()
        connection.close()
        ex.add(output,description=set_file_descr("meanScorePerFeature_"+density_name+".sql",groupId=grpId,step="norm_counts_per_frag",type="sql",view="admin",gdv='1'))

	countsPerFragFile=unique_filename_in()+".bed"
	with track.load(output,'sql') as t:
		t.convert(countsPerFragFile,'bed')
	ex.add(countsPerFragFile,description=set_file_descr("meanScorePerFeature_"+density_name+".bed",groupId=grpId,step="density",type="bed",ucsc="1"))
	step += 1

	# calculate normalised score per fragments (segToFrag)
	res = call_segToFrag(ex, countsPerFragFile, regToExclude, script_path, via=via)
	ex.add(res,description=set_file_descr("res_segToFrag_"+density_name+".bedGraph",groupId=grpId,step="norm_counts_per_frag",type="bedGraph",view="admin",comment="rough"))
	[resBedGraph,resBedGraph_all]=parseSegToFrag(res)
	ex.add(resBedGraph,description=set_file_descr("res_segToFrag_"+density_name+".bedGraph",groupId=grpId,step="norm_counts_per_frag",type="bedGraph",view="admin",comment="bedGraph non-sorted"))
	ex.add(resBedGraph_all,description=set_file_descr("res_segToFrag_"+density_name+"_all_nonSorted.bedGraph",groupId=grpId,step="norm_counts_per_frag",type="bedGraph",view="admin",comment="all informative frags - null included - bedGraph non-sorted"))
	resBedGraph_all=call_sortOnCoord(ex,resBedGraph_all,via=via)
	ex.add(resBedGraph_all,description=set_file_descr("res_segToFrag_"+density_name+"_all.bedGraph",groupId=grpId,step="norm_counts_per_frag",type="bedGraph",view="admin",comment="all informative frags - null included -sorted bedGraph"))
	resBedGraph=call_sortOnCoord(ex,resBedGraph,via=via)
	headerFile=unique_filename_in();
	hfile=open(headerFile,'w')
	hfile.write('track type="bedGraph" name="'+density_name+' normalised counts per valid fragments" description="'+density_name+' normalised counts per valid fragments" visibility=full windowingFunction=maximum autoScale=off viewLimits=1:2000\n')
	hfile.close()
	sortedBedGraph=cat([headerFile,resBedGraph])
	ex.add(sortedBedGraph,description=set_file_descr("res_segToFrag_"+density_name+".bedGraph",groupId=grpId,step="norm_counts_per_frag",type="bedGraph",comment="bedGraph sorted",ucsc='1'))
	sortedBedGraph_sql=unique_filename_in()
	touch(ex,sortedBedGraph_sql)
	with track.load(sortedBedGraph,'bedGraph', chrmeta=assembly.chrmeta) as t:
                t.convert(sortedBedGraph_sql+".sql",'sql')
	ex.add(sortedBedGraph_sql+".sql",description=set_file_descr("res_segToFrag_"+density_name+".sql",groupId=grpId,step="norm_counts_per_frag",type="sql",view="admin",gdv="1",comment="bedGraph sorted"))
	step += 1
	return [output,countsPerFragFile,res,resBedGraph,sortedBedGraph,sortedBedGraph_sql]

# Main 
#-------------------------------------------#
# *** open the 4c-seq minilims and create execution
# *** 0.get/create the library 
# *** 1.when necessary, calculate the density file from the bam file (mapseq.parallel_density_sql)
# ### 2.calculate the count per fragment for each denstiy file with gFeatMiner:mean_score_by_feature to calculate)
def workflow_groups(ex, job, primers_dict, assembly, mapseq_files, mapseq_url, script_path='', via='lsf' ):
	'''
		# Main 
		#-------------------------------------------#
		# *** open the 4C-seq minilims and create execution
		# *** 0.get/create the library 
		# *** 1.when necessary, calculate the density file from the bam file (mapseq.parallel_density_sql)
		# ### 2.calculate the count per fragment for each denstiy file with gFeatMiner:mean_score_by_feature to calculate)
	'''
	global grpId
	global step
	processed={
		'lib' : {},
		'density' : {},
		'4cseq' : {}
		}
	
        job_groups=job.groups
	htss_mapseq = frontend.Frontend( url=mapseq_url )

	new_libs=[]
	
	for gid, group in job_groups.iteritems():
                grpId = gid
		reffile=createlib.get_libForGrp(ex,group,assembly,new_libs, job.id, gid)
#		reffile='/archive/epfl/bbcf/data/DubouleDaan/library_Nla_30bps/library_Nla_30bps_segmentInfos.bed'
		processed['lib'][gid]=reffile

		for rid,run in group['runs'].iteritems():
                        #job_mapseq=htss_mapseq.job(run['key'])
		
			if 'regToExclude' in primers_dict[mapseq_files[gid][rid]['libname']]:
                                regToExclude=primers_dict[mapseq_files[gid][rid]['libname']]['regToExclude']
				regToExclude=regToExclude=regToExclude.replace('\r','')
			else:
			        regToExclude=None
			print("regToExclude="+str(regToExclude))
                        if not job.options.get('compute_densities') or job.options.get('merge_strands') != 0:
				print("will call parallel_density_sql with bam:"+mapseq_files[gid][rid]['bam']+"\n")
				density_file=parallel_density_sql( ex, mapseq_files[gid][rid]['bam'],
								   assembly.chromosomes,
								   nreads=mapseq_files[gid][rid]['stats']["total"],
								   merge=0,
								   convert=False,
								   via=via )
				mapseq_files[gid][rid]['wig']['merged']=density_file+"merged.sql"
				print("name of density_file after parallel_density_sql="+density_file)
				print("density file:"+mapseq_files[gid][rid]['wig']['merged'])
                        else:
                                print("Will use existing density file:"+mapseq_files[gid][rid]['wig']['merged'])

			print("density files:")
			print(mapseq_files[gid][rid]['wig']['merged'])
			print("Will convert density file .sql to .wig")
			mapseq_wig = unique_filename_in()
			touch(ex,mapseq_wig)
			print("mapseq_wig filename will be:"+mapseq_wig)
			with track.load(mapseq_files[gid][rid]['wig']['merged'],'sql') as t:
                		t.convert(mapseq_wig+".bw",'bigWig')
			ex.add(mapseq_wig+".bw",description=set_file_descr("density_file_"+mapseq_files[gid][rid]['libname']+".bw",groupId=gid,step="density",type="bigWig",ucsc='1'))	

			ex.add(mapseq_files[gid][rid]['wig']['merged'],description=set_file_descr("density_file_"+mapseq_files[gid][rid]['libname']+".sql",groupId=gid,step="density",type="sql",view='admin',gdv="1"))
                        processed['density'][mapseq_files[gid][rid]['libname']]=mapseq_files[gid][rid]['wig']['merged']

			print("Will process to the main part of 4cseq module: calculate normalised counts per fragments from density file:"+mapseq_files[gid][rid]['wig']['merged'])
			resfiles=density_to_countsPerFrag(ex,mapseq_files[gid][rid]['wig']['merged'],mapseq_files[gid][rid]['libname'],assembly,reffile,regToExclude,ex.remote_working_directory+'/',script_path, via)
			processed['4cseq']=resfiles
			
			print("Will proceed to profile correction of file "+str(resfiles[4]))
			profileCorrectedFile=unique_filename_in()
			reportFile_profileCorrection=unique_filename_in()
			profileCorrection.nonblocking(ex,resfiles[4],primers_dict[mapseq_files[gid][rid]['libname']]['baitcoord'],mapseq_files[gid][rid]['libname'],profileCorrectedFile,reportFile_profileCorrection,script_path,via=via).wait()
		        ex.add(profileCorrectedFile,description=set_file_descr("res_segToFrag_"+mapseq_files[gid][rid]['libname']+"_profileCorrected.bedGraph",groupId=gid,step="profile_correction",type="bedGraph",comment="profile corrected data;bedGraph sorted",ucsc='1'))
			ex.add(reportFile_profileCorrection,description=set_file_descr("report_profileCorrection_"+mapseq_files[gid][rid]['libname']+".pdf",groupId=gid,step="profile_correction",type="pdf",comment="report profile correction"))
			profileCorrectedFile_sql=unique_filename_in()+".sql"
			with track.load(profileCorrectedFile,'bedGraph') as t:
				t.convert(profileCorrectedFile_sql,'sql')
			ex.add(profileCorrectedFile_sql,description=set_file_descr("res_segToFrag_"+mapseq_files[gid][rid]['libname']+"_profileCorrected.sql",groupId=gid,step="profile_correction",type="sql",comment="profile corrected data;bedGraph sorted; sql format",gdv="1",view="admin"))

			step += 1
	
			print("Will smooth data before and after profile correction")
        		nFragsPerWin=str(group['window_size'])
			print("Window size="+nFragsPerWin)
        		outputfile=unique_filename_in()
		        smoothFragFile(ex,resfiles[4],nFragsPerWin,mapseq_files[gid][rid]['libname'],outputfile,regToExclude,script_path)
			ex.add(outputfile,description=set_file_descr("res_segToFrag_"+mapseq_files[gid][rid]['libname']+"_smoothed_"+nFragsPerWin+"FragsPerWin.bedGraph",groupId=gid,step="smoothing",type="bedGraph",comment="smoothed data, before profile correction",ucsc='1'))

			smoothedFile_sql=unique_filename_in()+".sql"
			with track.load(outputfile,'bedGraph') as t:	
				t.convert(smoothedFile_sql,'sql')
			ex.add(smoothedFile_sql,description=set_file_descr("res_segToFrag_"+mapseq_files[gid][rid]['libname']+"_smoothed_"+nFragsPerWin+"FragsPerWin.sql",groupId=gid,step="smoothing",type="bedGraph",comment="smoothed data, before profile correction, sql format",gdv="1",view="admin"))		
	
        		outputfile_afterProfileCorrection=unique_filename_in()
		        smoothFragFile(ex,profileCorrectedFile,nFragsPerWin,mapseq_files[gid][rid]['libname']+"_[fromProfileCorrected]",outputfile_afterProfileCorrection,regToExclude,script_path)
			ex.add(outputfile_afterProfileCorrection,description=set_file_descr("res_segToFrag_"+mapseq_files[gid][rid]['libname']+"_profileCorrected_smoothed_"+nFragsPerWin+"FragsPerWin.bedGraph",groupId=grpId,step="smoothing",type="bedGraph",comment="smoothed data, after profile correction",ucsc='1'))

			smoothedFile_afterProfileCorrection_sql=unique_filename_in()+".sql"
			with track.load(outputfile_afterProfileCorrection,'bedGraph') as t: 
				t.convert(smoothedFile_afterProfileCorrection_sql,'sql')
			ex.add(smoothedFile_afterProfileCorrection_sql,description=set_file_descr("res_segToFrag_"+mapseq_files[gid][rid]['libname']+"_profileCorrected_smoothed_"+nFragsPerWin+"FragsPerWin.sql",groupId=grpId,step="smoothing",type="sql",comment="smoothed data, after profile correction,sql format",gdv="1",view="admin"))

		step=0
	return processed

