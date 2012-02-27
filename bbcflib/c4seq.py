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
from bbcflib.common import unique_filename_in, gzipfile
import sys, getopt, os, json, re, time
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
			if re.search(r'IsValid',s.split('\t')[2]) and not re.search(r'_and_',s.split('\t')[8]) and not re.search(r'BothRepeats',s.split('\t')[8]) and not re.search(r'notValid',s.split('\t')[8]):
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

@program
def mergeBigWig(bwfiles,resfile,assembly_name):
	return {'arguments': ['mergeBigWig.sh','-i',bwfiles,'-o',resfile,'-a',assembly_name],
                'return_value':None}


@program
def call_runDomainogram(infile,name,prefix,regCoord=None,wmaxDomainogram=500,wmax_BRICKS=50,skip=0,script_path='./'):
        time.sleep(60)
        if regCoord == None: regCoord=""
        return{'arguments': ["R","--vanilla","--no-restore","--slave","-f",script_path+"runDomainogram.R","--args",infile,name,prefix,regCoord,str(wmaxDomainogram),str(wmax_BRICKS),str(skip)],
                'return_value':prefix+".log"}


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
	time.sleep(60)	
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

	print("will call mean_score_by_feature for t1=%s (name=%s) and t2=%s" %(density_file,density_name,reffile))
        output = unique_filename_in()
        chlist = []
        from gMiner.operations.genomic_manip.scores import mean_score_by_feature
        with track.Track(density_file) as scores:
                with track.Track(reffile) as features:
                        with track.new(output,format='sql',chrmeta=assembly.chrmeta) as out:
                                for ch in scores:
                                        chlist.append(ch)
                                        out.write(ch,mean_score_by_feature()(
                                                        scores.read(ch),
                                                        features.read(ch,fields=['start', 'end', 'name'])),
                                                  fields=['start', 'end', 'name', 'score'])
	countsPerFragFile=unique_filename_in()+".bed"
	with track.load(output,'sql') as t:
		t.convert(countsPerFragFile,'bed')

#	gzipfile(ex,countsPerFragFile)
#	ex.add(countsPerFragFile+".gz",description=set_file_descr("meanScorePerFeature_"+density_name+".bed.gz",groupId=grpId,step="density",type="bed",view="admin"))
        connection = sqlite3.connect(output)
        cursor = connection.cursor()
	cursor.execute("UPDATE 'attributes' SET value='%s' WHERE key='%s'"%('quantitative','datatype'))        
	[cursor.execute("UPDATE '%s' SET name=''" %ch) for ch in chlist]
	connection.commit()
	cursor.close()
        connection.close()
        ex.add(output,description=set_file_descr("meanScorePerFeature_"+density_name+".sql",groupId=grpId,step="norm_counts_per_frag",type="sql",view="admin",gdv='1'))

	step += 1

	# calculate normalised score per fragments (segToFrag)
	res = call_segToFrag(ex, countsPerFragFile, regToExclude, script_path, via=via)
	[resBedGraph,resBedGraph_all]=parseSegToFrag(res)

	gzipfile(ex,countsPerFragFile)
	ex.add(countsPerFragFile+".gz",description=set_file_descr("meanScorePerFeature_"+density_name+".bed.gz",groupId=grpId,step="density",type="bed",view="admin"))
	gzipfile(ex,res)
	ex.add(res+".gz",description=set_file_descr("res_segToFrag_"+density_name+".bedGraph.gz",groupId=grpId,step="norm_counts_per_frag",type="bedGraph",view="admin",comment="rough"))
#	ex.add(resBedGraph,description=set_file_descr("res_segToFrag_"+density_name+".bedGraph",groupId=grpId,step="norm_counts_per_frag",type="bedGraph",view="admin",comment="bedGraph non-sorted"))
#	ex.add(resBedGraph_all,description=set_file_descr("res_segToFrag_"+density_name+"_all_nonSorted.bedGraph",groupId=grpId,step="norm_counts_per_frag",type="bedGraph",view="admin",comment="all informative frags - null included - bedGraph non-sorted"))
	resBedGraph_all=call_sortOnCoord(ex,resBedGraph_all,via=via)
	gzipfile(ex,resBedGraph_all,args=["-c"],stdout=resBedGraph_all+".gz")
	ex.add(resBedGraph_all+".gz",description=set_file_descr("res_segToFrag_"+density_name+"_all.bedGraph.gz",groupId=grpId,step="norm_counts_per_frag",type="bedGraph",view="admin",comment="all informative frags - null included -sorted bedGraph"))
	
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
	return [output,countsPerFragFile+".gz",res+".gz",resBedGraph,sortedBedGraph,sortedBedGraph_sql,resBedGraph_all]

# Main 
#-------------------------------------------#
# *** open the 4c-seq minilims and create execution
# *** 0.get/create the library 
# *** 1.when necessary, calculate the density file from the bam file (mapseq.parallel_density_sql)
# ### 2.calculate the count per fragment for each denstiy file with gFeatMiner:mean_score_by_feature to calculate)
def workflow_groups(ex, job, primers_dict, assembly, mapseq_files, mapseq_url, script_path='', logfile=None, via='lsf' ):
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

        if logfile is None: logfile = sys.stdout
	
	for gid, group in job_groups.iteritems():
                grpId = gid
		reffile=createlib.get_libForGrp(ex,group,assembly,new_libs, job.id, gid)
#		reffile='/archive/epfl/bbcf/data/DubouleDaan/library_Nla_30bps/library_Nla_30bps_segmentInfos.bed'
		processed['lib'][gid]=reffile

		bwFiles=""
		for rid,run in group['runs'].iteritems():
			if 'regToExclude' in primers_dict.get(group['name'],{}):
                                regToExclude=primers_dict[group['name']]['regToExclude']
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

                        # Add here: mergeBigWig files
			bwFiles=bwFiles+mapseq_wig+".bw,"
                        # run the rest at the grp level

		# back to grp level!
		logfile.write("density files of replicates for group "+group['name']+ "("+str(gid)+")="+bwFiles);logfile.flush()
		mergedDensityFiles=unique_filename_in()
		mergeBigWig(ex,bwFiles,mergedDensityFiles,assembly.name)
                # convert result file to sql
		mergedDensityFiles_sql=unique_filename_in()
		with track.load(mergedDensityFiles,'bigWig') as t:
			t.convert(mergedDensityFiles_sql,'sql')
		ex.add(mergedDensityFiles,description=set_file_descr("density_file_"+group['name']+"_merged.bw",groupId=gid,step="density",type="bw",ucsc="1"))
		ex.add(mergedDensityFiles_sql,description=set_file_descr("density_file_"+group['name']+"_merged.sql ",groupId=gid,step="density",type="sql",gdv="1"))

		logfile.write("Will process to the main part of 4cseq module: calculate normalised counts per fragments from density file:"+mergedDensityFiles_sql);logfile.flush()

		#resfiles=density_to_countsPerFrag(ex,mapseq_files[gid][rid]['wig']['merged'],mapseq_files[gid][rid]['libname'],assembly,reffile,regToExclude,ex.remote_working_directory+'/',script_path, via)
		resfiles=density_to_countsPerFrag(ex,mergedDensityFiles_sql,group['name'],assembly,reffile,regToExclude,ex.remote_working_directory+'/',script_path, via)
		processed['4cseq']=resfiles
		
		logfile.write("Will proceed to profile correction of file "+str(resfiles[6]));logfile.flush()
		profileCorrectedFile=unique_filename_in()
		reportFile_profileCorrection=unique_filename_in()
		profileCorrection.nonblocking(ex,resfiles[6],primers_dict[group['name']]['baitcoord'],group['name'],
                                                    profileCorrectedFile,reportFile_profileCorrection,script_path,via=via).wait()
	        ex.add(profileCorrectedFile,description=set_file_descr("res_segToFrag_"+group['name']+"_profileCorrected.bedGraph",groupId=gid,step="profile_correction",type="bedGraph",comment="profile corrected data;bedGraph sorted",ucsc='1',gdv='1'))
		ex.add(reportFile_profileCorrection,description=set_file_descr("report_profileCorrection_"+group['name']+".pdf",groupId=gid,step="profile_correction",type="pdf",comment="report profile correction"))

		step += 1
	
		logfile.write("Will smooth data before and after profile correction (winSize="+str(group['window_size'])+" fragments per window)");logfile.flush()
       		nFragsPerWin=str(group['window_size'])
       		outputfile=unique_filename_in()
	        smoothFragFile(ex,resfiles[6],nFragsPerWin,group['name'],outputfile,regToExclude,script_path)
		ex.add(outputfile,description=set_file_descr("res_segToFrag_"+group['name']+"_smoothed_"+nFragsPerWin+"FragsPerWin.bedGraph",groupId=gid,step="smoothing",type="bedGraph",comment="smoothed data, before profile correction",ucsc='1',gdv='1'))

       		outputfile_afterProfileCorrection=unique_filename_in()
	        smoothFragFile(ex,profileCorrectedFile,nFragsPerWin,group['name']+"_fromProfileCorrected",outputfile_afterProfileCorrection,regToExclude,script_path)
		ex.add(outputfile_afterProfileCorrection,description=set_file_descr("res_segToFrag_"+group['name']+"_profileCorrected_smoothed_"+nFragsPerWin+"FragsPerWin.bedGraph",groupId=grpId,step="smoothing",type="bedGraph",comment="smoothed data, after profile correction",ucsc='1',gdv='1'))


		if not('run_domainogram' in group):
			group['run_domainogram'] = False
		elif str(group['run_domainogram']).lower() in ['1','true','on','t']:
			group['run_domainogram'] = True
		else:
			group['run_domainogram'] = False

		if not('before_profile_correction' in group):
			group['before_profile_correction'] = False
		elif str(group['before_profile_correction']).lower() in ['1','true','on','t']:
			group['before_profile_correction'] = True
		else:
			group['before_profile_correction'] = False

		if group['run_domainogram']:
			if regToExclude == None: regCoord=primers_dict[group['name']]['baitcoord'] 
			else: regCoord=regToExclude
			
			if group['before_profile_correction']:
				logfile.write("Will run domainogram from informative fragments (file:"+resfiles[6]+")");logfile.flush()
				call_runDomainogram(ex,resfiles[6],group['name'],group['name'],regCoord,script_path=script_path)
			else:
				logfile.write("Will run domainogram from profile corrected data (file:"+profileCorrectedFile+")");logfile.flush()
				call_runDomainogram(ex,profileCorrectedFile,group['name'],group['name'],regCoord.split(':')[0],500,50,1,script_path=script_path)
	
		        resFiles=[]
        		startRead=0
			res_tar=tarfile.open(group['name']+"_domainogram.tar.gz", "w:gz")
		        with open(group['name']+".log",'rb') as f:
                		for s in f.readlines():
                        		s=s.strip('\n')
                        		if re.search('####resfiles####',s): startRead=1
					if startRead>0 and not re.search("RData",s) and not re.search('####resfiles####',s):
                        			resFiles.append(s)
						res_tar.add(s)
						if re.search("bedGraph$",s):
							ex.add(s,description=set_file_descr(s,groupId=1,step="domainograms",type="bedGraph",ucsc="1",gdv="1"))
				
			res_tar.close()
        		ex.add(group['name']+"_domainogram.tar.gz",description=set_file_descr(group['name']+"_domainogram.tar.gz",groupId=1,step="domainograms",type="tgz"))
	


	step=0
	return processed

