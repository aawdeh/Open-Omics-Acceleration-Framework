source config
source miniconda3/bin/activate distbwa

###########################################################
# Get number of ranks from lscpu
###########################################################

echo "localhost" > hostfile
num_nodes=1
lscpu > compute_config

num_cpus_per_node=$(cat compute_config | grep -E '^CPU\(s\)' | awk  '{print $2}')
num_socket=$(cat compute_config | grep -E '^Socket'| awk  '{print $2}')
num_numa=$(cat compute_config | grep '^NUMA node(s)' | awk '{print $3}')
num_cpus_all_node=`expr ${num_cpus_per_node} \* ${num_nodes}`
threads_per_core=$(cat compute_config | grep -E '^Thread' | awk  '{print $4}')
echo "#############################################"
echo "Total number of CPUs: $num_cpus_all_node"
echo "Number of sockets: "$num_socket
echo "Number of NUMA domains: "$num_numa

num_physical_cores_all_nodes=`expr ${num_cpus_all_node} / ${threads_per_core}`
num_physical_cores_per_nodes=`expr ${num_cpus_per_node} / ${threads_per_core}`
num_physical_cores_per_socket=`expr ${num_physical_cores_all_nodes} / ${num_socket}`
num_physical_cores_per_numa=`expr ${num_physical_cores_all_nodes} / ${num_numa}`
echo "Num physical cores per nodes: "$num_physical_cores_per_nodes
echo "Num physical cores per socket: "$num_physical_cores_per_socket
echo "Num physical cores per numa: "$num_physical_cores_per_numa

th=`expr ${num_physical_cores_per_numa} / 2`  #${num_physical_cores_per_numa}  ##20
if [ $th -le 10 ]
then
    th=${num_physical_cores_per_numa}
fi

while [ $num_physical_cores_per_nodes -gt $th ]
do
    num_physical_cores_per_nodes=`expr $num_physical_cores_per_nodes / 2`
done

num_physical_cores_per_rank=$num_physical_cores_per_nodes
total_num_ranks=`expr ${num_physical_cores_all_nodes} / ${num_physical_cores_per_rank}`

ranks_per_node=`expr ${total_num_ranks} / ${num_nodes}`
echo "Cores per node: "$num_physical_cores_per_nodes
echo "Total number of ranks: "${total_num_ranks}
echo "Ranks per node: "${ranks_per_node}
echo "#############################################"

INDIR=$INPUT_DIR
OUTDIR=$OUTPUT_DIR

#* ranks: Number of mpi process that we want the pipeline to run on
#* threads/shards: parameters to different tools in the pipeline, calculated as below
ppn=${ranks_per_node}

Sockets=$(cat compute_config | grep -E '^Socket\(s\)' | awk  '{print $2}')   #2
Cores=$(cat compute_config | grep -E '^Core\(s\)' | awk  '{print $4}')  #56
Thread=$(cat compute_config | grep -E '^Thread' | awk  '{print $4}')  #2
a=$(( $(( ${Cores}*${Thread}*${Sockets} / $ppn )) - 2*${Thread} ))   #24 (Four threads are removed for IO)

b=$(( $(( ${Cores}*${Sockets} )) / $ppn ))   #14

if [ $a -lt 1 ]
then
    echo 'Number of cpus are less to run the pipeline.'
    exit 0
fi

N=${total_num_ranks}
PPN=${ranks_per_node}
CPUS=$a
THREADS=$a
BINDING=socket

###########################################################
# Set parameters 
###########################################################

## parameters
READ1=${R1[@]}
READ2=${R2[@]}

echo "reads:"
echo "READ1 $READ1"
echo "READ2 $READ2"
echo "R1PREFIX $R1PREFIX"
echo "R2PREFIX $R2PREFIX"

WHITELIST=$WHITELIST
REF=$REF
sample_id=""
outfile=""

## STAR Parameters that might change
PARAMS1="--soloType Droplet \
--soloStrand $star_strand_mode"
echo $PARAMS1

## Default values
PARAMS2="--soloUMIlen $soloUMIlen \
--soloCBlen $soloCBlen \
--soloFeatures $soloFeatures \
--clipAdapterType CellRanger4 \
--outFilterScoreMin 30 \
--soloCBmatchWLtype 1MM_multi_Nbase_pseudocounts \
--soloUMIdedup 1MM_Directional_UMItools \
--outSAMtype BAM SortedByCoordinate \
--outSAMattributes UB UR UY CR CB CY NH GX GN sF \
--soloBarcodeReadLength 0 \
--soloCellReadStats Standard \
--soloMultiMappers $soloMultiMappers \
--outStd BAM_SortedByCoordinate" \

echo $PARAMS2
## --outSAMattributes UB UR UY CR CB CY NH GX GN sF \
     
[[ -n $SAMPLE_ID ]] && sample_id="--sample_id $SAMPLE_ID" && echo "sample_id: $sample_id"
[[ -n $OUTFILE ]] && outfile="--outfile $OUTFILE" && echo "outfile: $outfile"

echo "Input directory: $INDIR"
echo "Output directory: $OUTDIR"
mkdir -p ${OUTDIR}

echo Starting run with $N ranks, $CPUS threads,$THREADS threads, $PPN ppn.
# -in -sindex are required only once for indexing.
# Todo : Make index creation parameterized.

###########################################################
# Error messages
###########################################################
# Check if number of ranks equals number of splits
echo $R1; 
echo $R2; 

R1_LEN=`echo $R1 | tr ' ' '\n' | wc -l`
R2_LEN=`echo $R2 | tr ' ' '\n' | wc -l`

if [ "$N" != "$R1_LEN" ]; then
    echo "Error: Number of ranks ("$N") does not equal number of splits ("$R1_LEN"). Program failed."
    exit 1
fi

# Check if number of R1 and R3 fastq files is equal
if [ "$R1_LEN" != "$R2_LEN" ]; then
    echo "Error: Number of R1 fastq files doesnt equal number of R2 files. Program failed."
    exit 1
fi

###########################################################
# Call dist_star.py 
###########################################################

exec=dist_star.py
mpiexec -bootstrap ssh -bind-to $BINDING -map-by $BINDING --hostfile hostfile -n $N -ppn $PPN python -u \
$exec --input $INDIR --output  $OUTDIR $TEMPDIR $REFDIR --index $REF --read1 $READ1 --read2 $READ2 \
--cpus $CPUS --threads $THREADS --keep_unmapped \
--whitelist $WHITELIST --reference_genome $REF --params1 "$PARAMS1" --params2 "$PARAMS2" \
${outfile} ${istart} ${sample_id} ${output_format} --r1prefix $R1PREFIX --r2prefix $R2PREFIX --suffix $SUFFIX \
2>&1 | tee ${OUTDIR}log.txt




