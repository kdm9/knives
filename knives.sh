#!/bin/bash

# REQUIREMENTS:
# - sickle, scythe, seqqs and pairs
# - gzip
# - awk
# - which
# - date

#############################
##   Print a help message  ##
#############################
function knives_help() {
	cat <<EOF
USAGE:
    knives.sh [options] -a CONTAM -i IN1 -o OUT1 -I IN2 -O OUT2 -u UNPAIRED -p PREFIX -r PRIOR

OPTIONS:

Mandatory parameters:
    -a CONTAM	Fasta file containing adaptor sequences, for scythe.
    -i IN1	Input R1 file.
    -I IN2	Input R2 file.
    -o OUT1	Output R1 file.
    -O OUT2	Output R2 file.
    -u UNPAIRED	Unpaired ouput file.
    -p PREFIX	Prefix of ouput files, including seqqs reports and log files.
    -r PRIOR	Prior probabilty of adaptor contamination. See scythe docs.

Optional parameters:
    -q QUAL	Minimum quality score for windowed quality trimming
                [default: 20]
    -t TYPE	Quality score encoding. One of phred, sanger, solexa or
                illumina. See the help message of scythe for an explanation.
                [default: sanger]
    -l LEN	Minimum length of sequence remaining after trimming.
                [default: 40]
    -s SKIP	Skip steps. Comma seperated list of program names. Must be one
                or more of: scythe, sickle, seqqs.
    -S		Don't log stderr of programs.
    -n		Remove any sequence with Ns.
    -Q QCDIR	Directory to store QC reports. [Default: .]
    -L LOGDIR	Directory to store logs. [Default: .]
    -K SICKLE	Location of the sickle executable.
    -C SCYTHE	Location of the scythe executable.
    -P PAIRS	Location of the pairs executable.
    -E SEQQS	Location of the seqqs executable.
EOF
}


#############################
##     Parse arguments     ##
#############################
now="$(date +%Y-%m-%d_%H-%M-%S)"
log_stderr=1
filter_ns=0
logdir="."
qcdir="."
sickle="$(which sickle)"
scythe="$(which scythe)"
pairs="$(which pairs)"
seqqs="$(which seqqs)"
qualtype="sanger"
qual="20"
minlen="40"
skip=""
while getopts "a:i:I:o:O:u:p:r:K:C:P:E:q:t:l:s:L:Sn" opt
do
	case $opt in
		a)
			contam="$OPTARG"
			;;
		i)
			infile1="$OPTARG"
			;;
		I)
			infile2="$OPTARG"
			;;
		o)
			outfile1="$OPTARG"
			;;
		O)
			outfile2="$OPTARG"
			;;
		u)
			outfileU="$OPTARG"
			;;
		p)
			prefix="$OPTARG"
			;;
		r)
			prior="$OPTARG"
			;;
		K)
			sickle="$OPTARG"
			;;
		C)
			scythe="$OPTARG"
			;;
		P)
			pairs="$OPTARG"
			;;
		E)
			seqqs="$OPTARG"
			;;
		q)
			qual="$OPTARG"
			;;
		t)
			qualtype="$OPTARG"
			;;
		l)
			minlen="$OPTARG"
			;;
		s)
			skip="$OPTARG"
			;;
		Q)
			qcdir="$OPTARG"
			if [ ! -d $qcdir ]
			then
				echo "[main] ERROR: QC dir must exist (which '$qcdir' does not)." >&2
				knives_help >&2
				exit -1
			fi
			;;
		L)
			logdir="$OPTARG"
			if [ ! -d $logdir ]
			then
				echo "[main] ERROR: Log dir must exist (which '$logdir' does not)." >&2
				knives_help >&2
				exit -1
			fi
			;;
		S)
			log_stderr=0
			;;
		n)
			filter_ns=1
			;;
		\?)
			echo "[main] ERROR: invalid option '-$OPTARG'" >&2
			knives_help >&2
			exit -1
			;;
	esac
done

if [ ! -e "$contam" ]
then
	echo "[main] ERROR: Adaptor file must exist (which '$contam' doesn't)." >&2
	knives_help >&2
	exit -1
fi

if [ ! -e "$infile1" ]
then
	echo "[main] ERROR: R1 input file must exist (which '$infile1' doesn't)." >&2
	knives_help >&2
	exit -1
fi

if [ ! -e "$infile2" ]
then
	echo "[main] ERROR: R2 input file must exist (which '$infile2' doesn't)." >&2
	knives_help >&2
	exit -1
fi

if [ -z "$outfile1" ] || [ ! -d "$(dirname ${outfile1})" ]
then
	echo "[main] ERROR: R1 output file must be createable (which '$outfile1' isn't)." >&2
	knives_help >&2
	exit -1
fi

if [ -z "$outfile2" ] || [ ! -d "$(dirname ${outfile2})" ]
then
	echo "[main] ERROR: R2 output file must be createable (which '$outfile2' isn't)." >&2
	knives_help >&2
	exit -1
fi

if [ -z "$outfileU" ] || [ ! -d "$(dirname ${outfileU})" ]
then
	echo "[main] ERROR: Unpaired read output file must be createable (which '$outfileU' isn't)." >&2
	knives_help >&2
	exit -1
fi

if [ -z "$prefix" ]
then
	echo "[main] ERROR: prefix must be provided and non-empty." >&2
	knives_help >&2
	exit -1
fi

if [ -z "$prior" ] || [ -z "$(echo ${prior} | awk '{if (0.0 < $1 && $1 < 1.0) print $1;}')" ]
then
	echo "[main] ERROR: prior must be provided and between 0 and 1." >&2
	knives_help >&2
	exit -1
fi

#############################
##     Setup Logging       ##
#############################
# log for commands
cmd_log="${logdir}/${prefix}_${now}_cmds.log"
# log file for the stderr of each process
if [ $log_stderr -eq 1 ]
then
	stderr_log="${logdir}/${prefix}_${now}_stderr.log"
else
	stderr_log="/dev/null"
fi
# log for a run summary
summary_log="${logdir}/${prefix}_${now}_summary.log"

#############################
##    assemble pipeline    ##
#############################
fifo_dir="$(mktemp -d)"
fifocmd="mkfifo"

# Interleave FQ files
il_fifo="${fifo_dir}/il.fq"
$fifocmd $il_fifo
il_cmd="${pairs} join ${infile1} ${infile2}"
echo "$il_cmd >${il_fifo}" >>${cmd_log}
$il_cmd >${il_fifo} 2>>${stderr_log} &

# Run seqqs on interleaved files
seqqs_il_fifo="${fifo_dir}/seqqs_il.fq"
$fifocmd $seqqs_il_fifo
seqqs_il_cmd="${seqqs} -e -i -p ${qcdir}/${prefix}_initial ${il_fifo}"
echo "$seqqs_il_cmd >${seqqs_il_fifo}" >>${cmd_log}
$seqqs_il_cmd >${seqqs_il_fifo} 2>>${stderr_log} &

# Run scythe
scythe_fifo="${fifo_dir}/scythe.fq"
$fifocmd $scythe_fifo
scythe_cmd="${scythe} -a ${contam} -p ${prior} -q ${qualtype} ${seqqs_il_fifo}"
echo "$scythe_cmd >${scythe_fifo}" >>${cmd_log}
$scythe_cmd >${scythe_fifo} 2>>${stderr_log} &

# Run seqqs on scythe'd files
seqqs_na_fifo="${fifo_dir}/seqqs_na.fq"
$fifocmd $seqqs_na_fifo
seqqs_na_cmd="${seqqs} -e -i -p ${qcdir}/${prefix}_noadapt ${scythe_fifo}"
echo "$seqqs_na_cmd >${seqqs_na_fifo}" >>${cmd_log}
$seqqs_na_cmd >${seqqs_na_fifo} 2>>${stderr_log} &

# Run sickle
sickle_fifo="${fifo_dir}/sickle.fq"
$fifocmd $sickle_fifo
sickle_singles_fifo="${fifo_dir}/sickle_singles.fq"
$fifocmd $sickle_singles_fifo
sickle_cmd="${sickle} pe -q ${qual} -l 1 -t ${qualtype} -c ${seqqs_na_fifo} -m ${sickle_fifo} -s ${sickle_singles_fifo}"
echo "$sickle_cmd" >>${cmd_log}
$sickle_cmd 2>>${stderr_log} &

# Run seqqs on sickle'd files
seqqs_qt_fifo="${fifo_dir}/seqqs_qt.fq"
$fifocmd $seqqs_qt_fifo
seqqs_qt_cmd="${seqqs} -e -i -p ${qcdir}/${prefix}_noadapt ${sickle_fifo}"
echo "$seqqs_qt_cmd >${seqqs_qt_fifo}" >>${cmd_log}
$seqqs_qt_cmd >${seqqs_qt_fifo} 2>>${stderr_log} &

# Deinterleave FQ files
r1o_fifo="${fifo_dir}/r1o.fq"
r2o_fifo="${fifo_dir}/r2o.fq"
rUo_fifo="${fifo_dir}/rUo.fq"
$fifocmd $r1o_fifo
$fifocmd $r2o_fifo
$fifocmd $rUo_fifo
#de_il_cmd="${pairs} split -1 ${outfile1} -2 ${outfile2} -u ${outfileU} ${seqqs_qt_fifo}"
de_il_cmd="${pairs} split -1 ${r1o_fifo} -2 ${r2o_fifo} -u ${rUo_fifo} ${seqqs_qt_fifo}"
echo "$de_il_cmd" >>${cmd_log}
$de_il_cmd 2>>${stderr_log} &

# gzip outputs
gzip="gzip -c"
r1_cmd="$gzip < ${r1o_fifo}"
echo "$r1_cmd >${outfile1}" >>$cmd_log
$gzip < ${r1o_fifo} >${outfile1} 2>>${stderr_log} &

r2_cmd="$gzip < ${r2o_fifo}"
echo "$r2_cmd >${outfile2}" >>$cmd_log
$gzip < ${r2o_fifo} >${outfile2} 2>>${stderr_log} &

rU_cmd="$gzip < ${rUo_fifo}" 
echo "$rU_cmd >>${outfileU}" >>$cmd_log
$gzip < ${rUo_fifo} >>${outfileU} 2>>${stderr_log} &

rU_cmd="$gzip < ${rUo_fifo}" 
echo "$rU_cmd >${outfileU}" >>$cmd_log
$gzip < ${sickle_singles_fifo} >>${outfileU} 2>>${stderr_log} &

#############################
##         Clean Up        ##
#############################
wait $(jobs -l | awk '{printf("%s ", $2)}')
rm -rf $fifo_dir && echo "[main] Removed temp files. All done now!"
