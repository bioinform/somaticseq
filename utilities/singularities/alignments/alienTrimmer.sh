#!/bin/bash

set -e

OPTS=`getopt -o o: --long output-dir:,in-fq1:,in-fq2:,out-fq1:,out-fq2:,singleton:,adapters:,min-length:,out-script:,action -n 'alienTrimmer.sh'  -- "$@"`

if [ $? != 0 ] ; then echo "Failed parsing options." >&2 ; exit 1 ; fi

echo "$OPTS"
eval set -- "$OPTS"

MYDIR="$( cd "$( dirname "$0" )" && pwd )"

adapter_file='/opt/Trimmomatic/adapters/TruSeq3-PE-2.fa'
min_length=36
standalone=1
timestamp=$( date +"%Y-%m-%d_%H-%M-%S-%N" )
action=echo

while true; do
	case "$1" in
		-o | --output-dir )
			case "$2" in
				"") shift 2 ;;
				*)  outdir=$2 ; shift 2 ;;
			esac ;;

		--in-fq1 )
			case "$2" in
				"") shift 2 ;;
				*)  in_fq1=$2 ; shift 2 ;;
			esac ;;

		--in-fq2 )
			case "$2" in
				"") shift 2 ;;
				*)  in_fq2=$2 ; shift 2 ;;
			esac ;;

		--out-fq1 )
			case "$2" in
				"") shift 2 ;;
				*)  out_fq1=$2 ; shift 2 ;;
			esac ;;

		--out-fq2 )
			case "$2" in
				"") shift 2 ;;
				*)  out_fq2=$2 ; shift 2 ;;
			esac ;;

		--singleton )
			case "$2" in
				"") shift 2 ;;
				*)  singleton_file=$2 ; shift 2 ;;
			esac ;;

		--adapters )
			case "$2" in
				"") shift 2 ;;
				*)  adapter_file=$2 ; shift 2 ;;
			esac ;;
		
		--min-length )
			case "$2" in
				"") shift 2 ;;
				*)  min_length=$2 ; shift 2 ;;
			esac ;;

		--out-script )
			case "$2" in
				"") shift 2 ;;
				*)  min_length=$2 ; shift 2 ;;
			esac ;;

		--action )
			case "$2" in
				"") shift 2 ;;
				*)  action=$2 ; shift 2 ;;
			esac ;;

		-- ) shift; break ;;

		* ) break ;;
	esac
done



logdir=${outdir}/logs
mkdir -p ${logdir}

if [[ ${out_script_name} ]]
then
    out_script="${out_script_name}"
else
    out_script="${logdir}/alienTrim.${timestamp}.cmd"
fi

if [[ $standalone ]]
then
    echo "#!/bin/bash" > $out_script
    echo "" >> $out_script
    echo "#$ -o ${logdir}" >> $out_script
    echo "#$ -e ${logdir}" >> $out_script
    echo "#$ -S /bin/bash" >> $out_script
    echo '#$ -l h_vmem=36G' >> $out_script
    echo 'set -e' >> $out_script
fi


echo "" >> $out_script

if [ "${in_fq1##*.}" = "gz" ]; then
	fq1=`basename ${in_fq1}`
	fq1="${fq1%.gz}"
	echo "gunzip -c ${in_fq1} > ${outdir}/${fq1}" >> $out_script
	in_fq1="${outdir}/${fq1}"
	intermediate_files="${outdir}/${fq1} $intermediate_files"
fi


if [ "${in_fq2##*.}" = "gz" ]; then
	fq2=`basename ${in_fq2}`
	fq2="${fq2%.gz}"
	echo "gunzip -c ${in_fq2} > ${outdir}/${fq2}" >> $out_script
	in_fq2="${outdir}/${fq2}"
	intermediate_files="${outdir}/${fq2} $intermediate_files"
fi

echo "" >> $out_script

echo "singularity exec --bind /:/mnt   docker://lethalfang/alientrimmer:0.4.0 \\" >> $out_script
echo "/opt/AlienTrimmer_0.4.0/src/AlienTrimmer \\" >> $out_script
echo "-if /mnt/${in_fq1} -ir /mnt/${in_fq2} \\" >> $out_script
echo "-c ${adapter_file} \\" >> $out_script
echo "-of /mnt/${outdir}/${out_fq1} -or /mnt/${outdir}/${out_fq2} \\" >> $out_script
echo "-os /mnt/${outdir}/${singleton_file} \\" >> $out_script
echo "-l ${min_length}" >> $out_script

echo '' >> $out_script

echo "bgzip ${outdir}/${out_fq1}" >> $out_script
echo "bgzip ${outdir}/${out_fq2}" >> $out_script
echo "bgzip ${outdir}/${singleton_file}" >> $out_script

echo '' >> $out_script
echo "for file in $intermediate_files" >> $out_script
echo '    do rm -v $file' >> $out_script
echo "done" >> $out_script


${action} $out_script
