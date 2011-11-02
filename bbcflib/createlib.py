"""
======================
Module: bbcflib.createlib.py 
======================

Functions used for the creation of a library in a 4c-seq analysis.
"""

from bein import *
from bein.util import *
from bbcflib.track import Track, new
from bbcflib import genrep,frontend,common
import sys, getopt, os, json, re
import tarfile


def load_libraryParamsFile(paramsfile):
	'''
        Returns a dictionary with the parameters required for the creation of a new library
	'''
	paramslib={}
	with open(paramsfile) as f:
		for s in f.readlines():
			s=s.strip('\n')
			if re.search('Library name',s.split('=')[0]):
				paramslib['name']=s.split('=')[1]
			elif re.search('Genome name',s.split('=')[0]):
				paramslib['specie']=s.split('=')[1]
			elif re.search('Primary',s.split('=')[0]):
				paramslib['primary']=s.split('=')[1]
			elif re.search('Secondary',s.split('=')[0]):
				paramslib['secondary']=s.split('=')[1]
			elif re.search('Segment length',s.split('=')[0]):
				paramslib['length']=s.split('=')[1]
			elif re.search('Type',s.split('=')[0]):
				paramslib['type']=s.split('=')[1]
	if len(paramslib['name']) < 2:
		paramslib['name']='myLibrary'
	if len(paramslib['length']) < 1:
		paramslib['length']='30'
	if 'type' not in paramslib or len(paramslib['type']) < 2 :
		paramslib['type']='typeI'	
	return paramslib

 
@program
def call_getRestEnzymeOccAndSeq(assembly_or_fasta,prim_site,sec_site,l_seg,g_rep, remote_working_directory, l_type='typeI'):
	'''
		Will create segments and fragments files of the new library from the genome sequence (via a call to getRestEnzymeOccAndSeq.pl)
		The genome sequence (assembly_or_fasta) can be given either as a genrep assembly or as a fasta file.  
	'''
	print(assembly_or_fasta)
	print(isinstance(assembly_or_fasta,genrep.Assembly))
	if isinstance(assembly_or_fasta,genrep.Assembly):
		print(assembly_or_fasta)
		print(isinstance(assembly_or_fasta,genrep.Assembly))
		print("Will prepare fasta file")
		fasta_path=g_rep.fasta_path(assembly_or_fasta)
		tar = tarfile.open(fasta_path) 
		tar.extractall(path=remote_working_directory)
		allfiles=[]
		for finfo in tar.getmembers():
        		if not finfo.isdir():
                		allfiles.append(remote_working_directory+finfo.name)
		fasta_file=common.cat(allfiles) 
		tar.close()
	else:
		fasta_file=assembly_or_fasta
	
	print('will call getRestEnzymes with fasta_file:'+fasta_file)
	segFile = unique_filename_in()
	fragFile = unique_filename_in()
	logFile = unique_filename_in()
	outfiles=[segFile, fragFile, logFile,fasta_file]
	print("outfiles will be:"+outfiles[0])
#	script_path='/archive/epfl/bbcf/mleleu/pipeline_vMarion/pipeline_3Cseq/vWebServer_SAM/'
	if cmp(l_type,'typeI') == 0:
		print('type:typeI')
		return {'arguments': ["getRestEnzymeOccAndSeq.pl","-i",fasta_file,"-m",prim_site,"-s",sec_site,"-l",l_seg,"-o",segFile,"-f",fragFile,"-x",logFile],
                        'return_value':outfiles}
	else:
                return {'arguments': ["getRestEnzymeOccAndSeq_typeII.pl","-i",fasta_file,"-m",prim_site,"-s",sec_site,"-l",l_seg,"-o",segFile,"-f",fragFile,"-x",logFile],
                        'return_value':outfiles}


def parse_fragFile(fragfile):
	'''
		Parse fragment file to create segment info bed file and fragment bed file	
	'''
	segInfoBedFile=unique_filename_in()
	fragmentBedFile=unique_filename_in()
	segmentBedFile=unique_filename_in()
	o=open(segInfoBedFile,'w')
	obed=open(fragmentBedFile,'w')
	oseg=open(segmentBedFile,'w')
	with open(fragfile,'r') as f:
		i=0
		for s in f.readlines():
                	s=s.strip('\n')
			i=i+1
			if i>1 and not re.search('FragIsNotValid',s):
				s_split=s.split('\t')
				fragmentInfo=s_split[0]+'|chr'+s_split[1]+':'+str(int(s_split[2])+1)+'-'+s_split[3]+'|indexOfSecondRestSiteOcc='+s_split[10]+'|status='+s_split[len(s_split)-1]+'|length='+str(int(s_split[3])-int(s_split[2]))+'|0|0|0|0'
				o.write('chr'+s_split[1]+'\t'+s_split[5]+'\t'+s_split[6]+'\ttype=startSegment|'+fragmentInfo+'\n')
				o.write('chr'+s_split[1]+'\t'+s_split[8]+'\t'+s_split[9]+'\ttype=endSegment|'+fragmentInfo+'\n')
				obed.write('chr'+s_split[1]+'\t'+s_split[2]+'\t'+s_split[3]+'\tfrag'+s_split[0]+'\n')
				oseg.write('chr'+s_split[1]+'\t'+s_split[5]+'\t'+s_split[6]+'\tfrag'+s_split[0]+'_startSeq\n')
				oseg.write('chr'+s_split[1]+'\t'+s_split[8]+'\t'+s_split[9]+'\tfrag'+s_split[0]+'_endSeq\n')
	o.close()
	obed.close()
	oseg.close()
	return([segInfoBedFile,fragmentBedFile,segmentBedFile])


@program 
def call_coverageBed(file1,file2):
	"""Binds ```coverageBed`` from the 'BedTools' suite
	"""
	print('will call coverageBed -a '+file1+' -b '+file2)
	return {"arguments": ['coverageBed','-a',file1,'-b',file2], "return_value":None}

def getCoverageInRepeats(ex,infile,genomeName='mm9',via='lsf'):
	'''
	Will complete the segment info bed file with the coverage in repeats of each segment 	
	For now, works only for mm9, hg19 and dm3 	
	'''
	repeatsPath="/archive/epfl/bbcf/data/genomes/repeats/"
	repeatsFile=repeatsPath+'/'+genomeName+'/'+genomeName+'_rmsk.bed'
	print("repeatsFile="+repeatsFile)

	if not(os.path.exists(repeatsFile)):
		print("coverage in repeats not calculated as file "+repeatsFile+" does not exist.")
		return(infile)

	tmpfile = unique_filename_in()
	_ = call_coverageBed.nonblocking(ex,repeatsFile,infile,via=via,stdout=tmpfile).wait()
	resfile = unique_filename_in()
	resfile = resfile+".bed"
	o=open(resfile,'w')
	with open(tmpfile,'r') as f:
		for s in f.readlines():
			s=s.strip('\n')
			s_split=(s.split('\t')[3]).split('|')
			infos='|'.join(s_split[0:(len(s_split)-4)])
			o.write('\t'.join((s.split('\t')[0:3]))+'\t'+infos+'|'+'|'.join((s.split('\t')[4:8]))+'\n')
	o.close()	
	return resfile

def getEnzymeSeq(enzyme_id,enzymes_dict=None):
	'''
		Return the restriction site corresponding to a given enzyme id (from existing enzymes)
	'''
	if enzymes_dict == None or not isinstance(enzymes_dict,list):
		enzymes='/archive/epfl/bbcf/mleleu/pipeline_vMarion/pipeline_3Cseq/vWebServer_Bein/tests/enzymes.json'
	        g=open(enzymes)
        	enzymes_dict=json.load(g)
	for x in enzymes_dict:
		for k,v in x.iteritems():
			if v['id']==enzyme_id:
				return(v['site'])
	
	return None

def getEnzymeId(enzyme_seq,enzymes_dict=None):
	'''
		Return the enzyme id corresponding to a given restriction site sequence (from existing enzymes).
	'''
	if enzymes_dict == None or not isinstance(enzymes_dict,list):
		enzymes='/archive/epfl/bbcf/mleleu/pipeline_vMarion/pipeline_3Cseq/vWebServer_Bein/tests/enzymes.json'
		g=open(enzymes)
		enzymes_dict=json.load(g)
	for x in enzymes_dict:
		for k,v in x.iteritems():
			if v['site'] == enzyme_seq:		
				return(v['id'])
	return 0
	
def lib_exists(params,libs_dict=None,returnType="id"):
	'''
		Return id or filename corresponding to the library described in params. 	
	'''
	if libs_dict == None :	
		libs='/archive/epfl/bbcf/mleleu/pipeline_vMarion/pipeline_3Cseq/vWebServer_Bein/tests/libraries.json'
		f=open(libs)
		libs_dict = json.load(f)
	enzymes='/archive/epfl/bbcf/mleleu/pipeline_vMarion/pipeline_3Cseq/vWebServer_Bein/tests/enzymes.json' 
	g=open(enzymes)
	enzymes_dict=json.load(g)
	for lib in libs_dict: 
		print lib['library']['assembly_name']
		enz1=getEnzymeSeq(lib['library']['enzyme1_id'],enzymes_dict)
		enz2=getEnzymeSeq(lib['library']['enzyme2_id'],enzymes_dict)
		#temporary...
		if not 'type' in lib['library']: lib['library']['type']='typeI'
		if params['specie']==lib['library']['assembly_name'] and params['primary']==enz1 and params['secondary']==enz2 and int(params['length'])==int(lib['library']['segment_length']) and params['type']==lib['library']['type']:
			if returnType=="id":
				return lib['library']['id']
			else:
				return lib['library']['filename']
		
	if returnType == "id" :
		return 0
	else:
		return None

def get_libfile(id_lib):
	libs='/archive/epfl/bbcf/mleleu/pipeline_vMarion/pipeline_3Cseq/vWebServer_Bein/tests/libraries.json'
        f=open(libs)
        libs_dict = json.load(f)
	#id_lib=13
	for lib in libs_dict:
		if lib['library']['id']==int(id_lib):
			 return lib['library']['filename']
	return None

@program
def call_gunzip(path):
	"""
	unzip gzip file
	"""
	output=unique_filename()
	call = ["gunzip",path]
	return {"arguments": call, "return_value": output}

# *** main call to create the library
def createLibrary(ex,fasta_allchr,params, g_rep):
	'''
		main call to create the library
	'''	
	if len(params['primary'])<2 or len(params['secondary'])<2:
		print('Some parameters are missing, cannot create the library')
		print('primary='+params['primary']+" ; "+'secondary='+params['secondary'])
		return [None,None,None,None]
	print('Will call call_getRestEnzymeOccAndSeq')
	print(fasta_allchr)
	libfiles=call_getRestEnzymeOccAndSeq(ex,fasta_allchr,params['primary'],params['secondary'],params['length'],g_rep, ex.remote_working_directory + "/", params['type'])
	print('parse fragment file to create segment info bed file and fragment bed file\n')
	bedfiles=parse_fragFile(libfiles[1])
	print('calculate coverage in repeats for segments')
	resfile=getCoverageInRepeats(ex,bedfiles[0],params['specie'],via='local')
	print('convert bed segments infos file to sql format')
	resfile_sql=resfile+".sql"
	with Track(resfile,format="bed",chrmeta=params['specie']) as Lib:
		with new(resfile_sql, 'sql',chrmeta=params['specie']) as out:
			for ch in Lib:
				out.write(ch,Lib.read(ch,fields=['start', 'end', 'name']),fields=['start', 'end', 'name'])
	
	infos_lib={'assembly_name':params['specie'],'enzyme1_id':getEnzymeId(params['primary']),'enzyme2_id':getEnzymeId(params['secondary']),'segment_length':params['length'],'type':params['type'],'filename':resfile}
	return([libfiles,bedfiles,resfile,infos_lib,resfile_sql])


def get_libForGrp(ex,group,fasta_or_assembly,new_libraries, job_id, g_rep):
	#wd_archive="/archive/epfl/bbcf/mleleu/pipeline_vMarion/pipeline_3Cseq/vWebServer_Bein/" #temporary: will be /scratch/cluster/monthly/htsstation/4cseq/job.id
	lib_dir = "/scratch/cluster/monthly/htsstation/4cseq/" + str(job_id) + "/"
	print "Group:\n"
	print group
	if 'library_param_file' in group and group['library_param_file'] != "" :
		library_filename = lib_dir + 'group_' + group['name'] + "_paramsFileLibrary.txt"
		paramslib=load_libraryParamsFile(library_filename);
		lib_id=lib_exists(paramslib)
		ex_libfile=lib_exists(paramslib,new_libraries,returnType="filename")
		print("lib_id="+str(lib_id))
		if lib_id == 0 and ex_libfile == None :
			print("will call createlib.createLibrary with:"+str(fasta_or_assembly)+" and "+ library_filename)
			libfiles=createLibrary(ex,fasta_or_assembly,paramslib, g_rep);
			reffile=libfiles[4]
			ex.add(libfiles[2],description=common.set_file_descr("new_library.bed",group=group['name'],step=0,type="bed"))
			ex.add(reffile,description=common.set_file_descr("new_library.sql",group=group['name'],step=0,type="sql",view=admin))
			new_libraries.append({'library':libfiles[3]})	
		elif lib_id > 0 :
			print("This library already exists (id="+str(lib_id)+")")
			reffile=get_libfile(lib_id)
		else:
			print("This library has just been created ("+ex_libfile+")")
			reffile=ex_libfile
	elif 'library_id' in group and group['library_id'] > 0:
		reffile=get_libfile(group['library_id'])
		if reffile==None:
			raise TypeError("No valid parameter passed for the library.")
		elif not os.path.exists(reffile):
			reffile=reffile+'.sql'
		else:
			raise TypeError("library file ("+reffile+") is not valid")
	elif 'library_file_url' in group and group['library_file_url'] != "" :
		reffile=group['library_file_url']
	else:
		reffile=None

	print("reffile="+reffile)				 
	return reffile

