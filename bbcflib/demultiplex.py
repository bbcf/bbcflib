"""
===========================
Module: bbcflib.demultiplex
===========================

"""
from bein import *
from bein.util import *
import sys, getopt, os

# call 
#python demultiplex.py -i "/archive/epfl/bbcf/mleleu/pipeline_vMarion/pipeline_3Cseq/vWebServer_Bein/tests/export_part.fa" -p "/archive/epfl/bbcf/mleleu/pipeline_vMarion/pipeline_3Cseq/vWebServer_Bein/tests/primers.fa" -x 30 -n 1 -s 108 -l 50

grpId=1
step=0

@program
def exportToFasta(exportFile,n=1,x=22):
	faFile=unique_filename_in()
	scriptPath="/archive/epfl/bbcf/mleleu/pipeline_vMarion/pipeline_3Cseq/vWebServer_Bein/"
	return{'arguments': ["python",scriptPath+"exportToFasta.py","-i",exportFile,"-n",str(n),"-x",str(x)],
               'return_value':faFile}

@program
def fastqToFasta(fqFile,n=1,x=22):
	faFile=unique_filename_in()
	scriptPath="/archive/epfl/bbcf/mleleu/pipeline_vMarion/pipeline_3Cseq/vWebServer_Bein/"
        return{'arguments': ["python",scriptPath+"fastqToFasta.py","-i",fqFile,"-o",faFile,"-n",str(n),"-x",str(x)],
               'return_value':faFile}

def split_exonerate(filename,n=1,x=22,l=30):	    
	correction={}
	files={}
	filenames={}
	with open(filename,"r") as f:
		for s in f:
			s=s.strip('\n')
			if re.search(r'vulgar',s):
				key = s.split(' ')[5].split('|')[0]
				if not key in correction:
					correction[key]=len(s.split('|')[1])-len(s.split('|')[3])+1
					filenames[key]=unique_filename_in()					
					files[key]=open(filenames[key],"w")
				k=int(s.split(' ')[3])-int(s.split(' ')[7])+correction[key]
	#			print("key="+str(key)+";correction="+str(correction[key])+";k="+str(k))
				seq=s.split(' ')[1].split('_')[1][k:l+k]
				qual=s.split(' ')[1].split('_')[2][k:l+k]
				files[key].write("@"+s.split(' ')[1].split('_')[0]+"_"+seq+"_"+qual+"\n"+seq+"\n+\n"+qual+"\n")
	for f in files.itervalues():
		f.close()
	return filenames



@program
def raw_exonerate(fastaFile,dbFile,minScore=77,n=1,x=22,l=30):
	return{'arguments': ["exonerate","--model","affine:local","-o","-4","-e","-12","-s",str(minScore),fastaFile,dbFile],
	       'return_value': None} 

def exonerate(ex,infile, dbFile, minScore=77,n=1,x=22,l=30,via="local"):
#	filename = unique_filename_in()
	#kwargs["stdout"] = filename
	global grpId
	global step

        print "Will prepare fasta file from input file\n"
	subfiles=split_file(ex,infile,n_lines=8000000)
	faSubFiles=[]
	
        if re.search(r'\.fa$',infile) or re.search(r'\.fasta$',infile):
                faSubFiles=subfiles
        elif re.search(r'export',infile):
                futures=[exportToFasta.nonblocking(ex,exportFile=sf,n=n,x=x,via=via) for sf in subfiles]
		faSubFiles=[f.wait() for f in futures]
        elif re.search(r'\.fastq$',infile) or re.search(r'\.fq$',infile):
		futures=[fastqToFasta.nonblocking(ex,sf,n=n,x=x,via=via) for sf in subfiles]
                faSubFiles=[f.wait() for f in futures]

	for f in faSubFiles:
		ex.add(f,description="input.fasta (part) [group:" + str(grpId) + ",step:"+ str(step) + ",type:fa,view:admin]")
	
	print("Will call raw_exonerate for each fasta files")	
#	faSubFiles=split_file(ex,faFile,n_lines=500000)
	#print faSubFiles

	futures = []
	res = []
	resExonerate = []
	for sf in faSubFiles:
		subResFile = unique_filename_in() 
		futures.append(raw_exonerate.nonblocking(ex,sf,dbFile,minScore=minScore,n=n,x=x,l=l,via=via,stdout=subResFile))
		resExonerate.append(subResFile)
	#print resExonerate
	for f in futures: f.wait()
	for f in resExonerate:
		ex.add(f,description="exonerate (part) [group:" + str(grpId) + ",step:" + str(step) + ",type:txt,view:admin]")
		res.append(split_exonerate(f,n=n,x=x,l=l))
	step += 1
	#print res
 
	return res


def demultiplex(ex,fastaFile,dbFile,minScore=77,n=1,x=22,l=30,via="local"):
	resFiles={}
	print("in demultiplex, call exonerate\n")
	exonerated = exonerate(ex, fastaFile, dbFile, minScore=int(minScore),n=n,x=x,l=int(l),via=via)
	#print exonerated
#	ex.add(exonerated,description="fileExonerate")
#	split_files = split_exonerate(exonerated,n=n, x=x, l=l)
	print("demultiplex done\n")
	for sr in exonerated:
		for k,v in sr.iteritems():
               		#print("k="+k)
                        if k in resFiles:
				f=open(resFiles[k],"a+")
				f.write(open(v).read())
				f.close()
			else:
                        	resFiles[k]=v
	return resFiles


def getFileFromURL(file_loc,od="", suffix=None):
	suffix = (suffix and ("." + suffix)) or ''
	resfile=os.path.join(od,unique_filename_in()) + suffix
        if file_loc.startswith("http://") or file_loc.startswith("https://") or file_loc.startswith("ftp://"):
		urllib.urlretrieve( file_loc, resfile )
        elif os.path.exists(file_loc):
                shutil.copy( file_loc, resfile )
	else:
		raise ValueError("Couldn't find this file anywhere: %s." %file_loc)

	return resfile


def load_paramsFile(paramsfile):
	'''
		Return a dictionary with the parameters required for exonerate
	'''
	params={}
	with open(paramsfile) as f:
		for s in f.readlines():
			if not re.search('^#',s):
				s=s.strip('\n')
				(k,v)=s.split('=')
				if re.search('Search the primer from base i',k):
					params['n']=v
				if re.search('Search the primer in the next n bps of the reads',k):
					params['x']=v
				if re.search('Minimum score for Exonerate',k):
					params['s']=v
				if re.search('Generate fastq output files',k):
					if re.search('Y',v) or re.search('y',v):
						params['q']=True
					else:
						params['q']=False
				if re.search('Length of the reads to align',k):
					params['l']=v
	return params


@program
def count_lines(filename):
   def parse_output(p):
        m = re.search(r'^\s*(\d+)\s+' + filename, ''.join(p.stdout))
        return int(m.groups()[0])
   return {'arguments': ["wc","-l",filename],
           'return_value': parse_output}

def workflow_groups(ex, job, scriptPath):
	processed = {}
	job_groups=job.groups
	resFiles={}
	global grpId
	global step
	for gid, group in job_groups.iteritems():
		for rid,run in group['runs'].iteritems():
			print(group)
			print(run)
			lib_dir="/scratch/cluster/monthly/htsstation/demultiplexing/" + str(job.id) + "/"
			suffix = run['url'].split('.')[-1]
			print suffix
			infile=getFileFromURL(run['url'],lib_dir, suffix)
			primersFile = lib_dir + 'group_' + group['name'] + "_primer_file.fa"	
			ex.add(primersFile,description="group_" + group['name'] + "_primer_file.fa [group:"+ str(grpId) +",step:"+ str(step) + ",type:fa]" )
			paramsFile = lib_dir + 'group_' + group['name'] + "_param_file.txt"	
			ex.add(paramsFile,description="group_" + group['name'] + "_param_file.txt [group:"+ str(grpId) +",step:" + str(step) + ",type:txt]")
			params=load_paramsFile(paramsFile)
			print(params)	
			#demultiplex(ex,infile,opts['-p'],int(opts['-s']),opts['-n'],opts['-x'],opts['-l'],via="lsf")
			resExonerate = demultiplex(ex,infile,primersFile,params['s'],params['n'],params['x'],params['l'],via='lsf')
			
	        	filteredFastq={}
	        	for k,f in resExonerate.iteritems():
        		        ex.add(f,description="k:"+k+".fastq [group:" + str(grpId) + ",step:"+ str(step) + ",type:fastq,view:admin]")

			step += 1

		        print "Will filter the sequences\n"
		        filteredFastq=filterSeq(ex,resExonerate,primersFile)
			
		        for k,f in filteredFastq.iteritems():
        #                        resFiles[k]=[]
		                ex.add(f,description="k:"+k+"_filtered.fastq [group:" + str(grpId) + ",step:" + str(step) + ",type:fastq,view:admin]")
	#			if k in resFiles:
	#				resFiles[k].append(f)

			step += 1
			reportFile=unique_filename_in()
#		!! STILL HAVE TO GENERATE THE REPORT...	
		# for each fastq file, call count_lines(fastq)/4
	#		ex.add(reportFile,description="pdf:report_demultiplexing")
	#		resExonerate.append(reportFile)

		grpId += 1
		step = 0

	return resFiles 


@program
def call_createReport(numbersFile,reportFile):
	Rout=unique_filename_in()
	return{'arguments': ['R','CMD','BATCH','--vanilla','--no-restore','"--args',numbersFile, numbersFile ,reportFile,'"',scriptPath+"plotGraphsDemultiplexing.R",Rout],
        	'return_value':None}


def createReport(ex):
	countsFile=unique_filename_in()
	out=open(countsFile,'w')
	out.write('Primer\tTotal_number_of_reads\tnFiltered\tnValidReads\n')
	allfiles = common.get_files( ex.id, M )
	allcounts={}
#	for f in allfiles['fastq']
#		curPrimer=f[0:f.index(':')]
#		if curPrimer not in allcounts: allcounts[curPrimer]=[0,0,0]
#		if re.search(r'filtered',f):
#			allcounts[curPrimer][1]=count_lines(f)
#		else:
#			allcounts[curPrimer][0]=count_lines(f)
	
#	resFile=unique_filename_in()
#	call_createReport(countsFile,resFile)
#bsub -J "${runName}_generateReportDemultiplexing" -o "${wd}${runName}_generatReportDemultiplexing.log" -e "${wd}${runName}_generateReportDemultiplexing.err" "${RPath}R CMD BATCH --vanilla --no-restore '--args ${reportStep1} ${reportStep2} ${plotsDemultiplexing}' ${scriptsDirectory}plotGraphsDemultiplexing.R ${wd}${runName}_generateReportDemultiplexing.Rout"
	

def filterSeq(ex,fastqFiles,primersFile):
#seqToFilter=`awk -v primer=${curPrimer} '{if($0 ~ /^>/){n=split($0,a,"|");curPrimer=a[1];gsub(">","",curPrimer); if(curPrimer == primer){seqToFilter="";for(i=5;i<n;i++){if(a[i] !~ /Exclude=/){seqToFilter=seqToFilter""a[i]","}} if(a[n] !~ /Exclude=/){seqToFilter=seqToFilter""a[n]}else{gsub(/,$/,"",seqToFilter)};print seqToFilter}}}' ${primersFile}`
	seqToFilter={}
	filenames={}
	global grpId
	global step 

	with open(primersFile,"r") as f:
		for s in f:	
			if re.search(r'^>',s):
				key=s.split('|')[0].replace(">","")
				filenames[key]=unique_filename_in() 
				seqToFilter[key]=open(filenames[key],"w")
				for i in range(4,len(s.split('|'))):
					if not cmp(s.split('|')[i],'.') == 0 or not cmp(s.split('|')[i],'Exclude'):
						seqToFilter[key].write(">seq"+str(i)+"\n"+s.split('|')[i]+"\n")
	
	print "Will build bowtie indexes for seqToFilter\n" 
	indexSeqToFilter={}
	for k,f in seqToFilter.iteritems():
		f.close()
	
	indexFiles={}	
	for k,f in filenames.iteritems():
                ex.add(f,description=k+":"+k+"_seqToFilter.fa [group:"+ str(grpId) + ",step:" + str(step) + ",type:fa,view:admin]")
		#indexFiles[k]=add_bowtie_index(ex,f,description=k+'_bowtie index',stdout="../out")
		indexFiles[k]=bowtie_build.nonblocking(ex,f,via='lsf')

	print "Wait for index\n"
	bowtie_index={}
	unalignedFiles={}
	alignedFiles={}
	futures={}

	for k,f in filenames.iteritems():	
		if k in fastqFiles:	
			#bowtie_index[k]=indexFiles[k]
			bowtie_index[k]=indexFiles[k].wait()
			unalignedFiles[k]=unique_filename_in()
			#print("primer="+k+"=>index="+bowtie_index[k])
			print "Will filter reads (call bowtie)\n"
			bwtarg=["-a","-q","-n","2","-l","20","--un",unalignedFiles[k]]
			futures[k]=bowtie.nonblocking( ex, bowtie_index[k], fastqFiles[k],  
                                     bwtarg, via='lsf')
		else:
			print("key "+k+" is not defined in fastqFiles")
			#print fastqFiles

	for k,f in filenames.iteritems():
		alignedFiles[k]=futures[k].wait()

	return unalignedFiles



