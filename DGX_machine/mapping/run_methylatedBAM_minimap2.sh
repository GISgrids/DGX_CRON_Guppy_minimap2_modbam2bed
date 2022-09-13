#!/bin/bash

THREAD=$1
THREAD_LESS=$2 

##############################################################################################
RUN_DIR="/raid/scratch/lizh/methyl_test/mapping"
BAM_TMP="${RUN_DIR}/bam"
BAM_INPUT="/mnt/hexagon-scratch/lizh/methy_test/guppy6_result"
BAM_BACKUP="/mnt/hexagon-scratch/lizh/methy_test/bam_test"

FASTA="/mnt/hexagon-scratch/lizh/methy_test/grch38/GRCh38.primary_assembly.genome.fa"
MMI_INDEX="/mnt/hexagon-scratch/lizh/methy_test/grch38/grch38.mmi"

MONGO_EXEC="/raid/scratch/lizh/mongo/bin/mongo"
DB_LOG="${RUN_DIR}/mongo.log"

# Run admin - max run in a row and lock file name
MAX_RUN=100
LOCK_DIR="minimap2_cron.lock"

#### BACKUP ORIGINAL variable
# OUT_DIR="/mnt/hexagon-scratch/lizh/methy_test/bam" 
# MONGO_EXEC=$4
# BAM_DIR="/mnt/hexagon-scratch/lizh/methy_test/guppy6_result"
# LOGFILE="${OUT_DIR}/minimap2_cron.log"
# ERROR="${OUT_DIR}/minimap2_cron.ERROR"
# BACKUP_DIR="/mnt/hexagon-scratch/lizh/methy_test/bam"

#############################################################################################

cd ${RUN_DIR}

echo -e `date +"%F\t%T\t"`"[TEST]"; 

lockdir="${RUN_DIR}/${LOCK_DIR}"

# Create temp bam directory on DGX
[ ! -d ${BAM_TMP} ] && mkdir -p ${BAM_TMP} && echo "Created ${BAM_TMP}" 

declare -i i=0
while [ $i -lt ${MAX_RUN} ]
do
	if mkdir -- "$lockdir" ; then
		# printf >&2 'successfully acquired lock\n' >>${ERROR}
		
		EVAL='db.ARG_methylation_MUX.findOne({"Status": "DONE_BASECALL"})._id ; '
		echo "querySql: $EVAL" >> ${DB_LOG}
		DOC_ID=$(echo "$EVAL" | $MONGO_EXEC -u promappuser -p "q^a#W3LTjq?t^kwB" 54.251.49.171:8080/promethion --authenticationDatabase promethion --quiet)
		SUCCESS=$?
		
		if [ ${SUCCESS} -eq 0 ]; then 

			EVAL='db.ARG_methylation_MUX.findOne('"${DOC_ID}"').mux_full_path ; '
			echo "querySql: $EVAL" >> ${DB_LOG}
			MUX_FULL_PATH=$(echo "$EVAL" | $MONGO_EXEC -u promappuser -p "q^a#W3LTjq?t^kwB" 54.251.49.171:8080/promethion --authenticationDatabase promethion --quiet)
		
			# eg ON002-DNA-R00141/WHB12183-T1/20210208_0141_1G_PAG51949_51a655a0/fast5_fail/PAG51949_fail_7683c71e_99.fast5
			IFS='/' read -ra ARRAY <<< "${MUX_FULL_PATH}"
			RUN=${ARRAY[-3]}
			ELM=${ARRAY[-2]}
			INDEX=${ARRAY[-1]}
			echo -e "Processing: ${RUN}/${ELM}/${INDEX}" | gawk '{ print strftime("[%Y-%m-%d %H:%M:%S]"), $0 }' 
			
			# Update MUX_ID Status -> "mapping"
			UPDATE='db.ARG_methylation_MUX.updateOne( {"_id": '"${DOC_ID}"' } , { $set:{"Status": "mapping"} } ) ; '
			echo "UpdateSql: $UPDATE" >> ${DB_LOG}
			result=$(echo "$UPDATE" | $MONGO_EXEC -u promappuser -p "q^a#W3LTjq?t^kwB" 54.251.49.171:8080/promethion --authenticationDatabase promethion --quiet)
			echo "outcome: $result" >> ${DB_LOG}
			
			# cat QC passed reads 
			echo -e `date +"%F\t%T\t"`"[START] cat bam ${RUN}"; 
			INPUT_PREFIX="${BAM_TMP}/${RUN}_${INDEX}"
			samtools cat -@${THREAD} -o ${INPUT_PREFIX}.bam ${BAM_INPUT}/${MUX_FULL_PATH}/Guppy_6.2.1/pass/*.bam 
			echo -e `date +"%F\t%T\t"`"[END] cat bam ${RUN}"; 

			# cat QC failed reads
			# samtools cat -o ${OUT_DIR}/${RUN}.fail.bam \
			# ${BAM_DIR}/${MUX_1}/*/*/Guppy_6.2.1/fail/*.bam ${BAM_DIR}/${MUX_2}/*/*/Guppy_6.2.1/fail/*.bam

			echo -e `date +"%F\t%T\t"`"[START] minimap2 ${RUN}_${INDEX}.bam";
			#CMD="samtools fastq --threads ${THREAD_LESS} -T Mm,Ml -0 - ${INPUT_PREFIX}.bam | minimap2 -ax map-ont -y -t ${THREAD} ${MMI_INDEX} - | samtools sort --threads ${THREAD_LESS} -o ${INPUT_PREFIX}.minimap2.sorted.bam - "
			#echo -e ${CMD}
			samtools fastq --threads ${THREAD_LESS} -T Mm,Ml -0 - ${INPUT_PREFIX}.bam | minimap2 -ax map-ont -y -t ${THREAD} ${MMI_INDEX} - | samtools sort --threads ${THREAD_LESS} -o ${INPUT_PREFIX}.minimap2.sorted.bam -
			SUCCESS_MINIMAP=$?
			echo -e `date +"%F\t%T\t"`"[END] minimap2 ${RUN}";

			# Clean-up
			if [ ${SUCCESS_MINIMAP} -eq 0 ]; then 
				
				SUCCESS_MERGE=1
				
				# remove concate bam file
				rm ${INPUT}
				
				# Check if BAM file exist in output dir
				OUT_BAM="${BAM_BACKUP}/${RUN}.minimap2.sorted.bam"
				if [ -e ${OUT_BAM} ]; then
					# exist!...Merge!!!
					echo -e `date +"%F\t%T\t"`"[START] Merge ${RUN}";
					samtools merge --threads ${THREAD} --write-index -o ${BAM_TMP}/${RUN}.merged.minimap2.sorted.bam ${OUT_BAM} ${BAM_TMP}/${RUN}_${INDEX}.minimap2.sorted.bam
					SUCCESS_MERGE=$?
					echo -e `date +"%F\t%T\t"`"[END] Merge ${RUN}";
					
					if [ ${SUCCESS_MERGE} -eq 0 ]; then
						echo -e `date +"%F\t%T\t"`"[REMOVE] ${OUT_BAM} and ${BAM_TMP}/${RUN}_${INDEX}.minimap2.sorted.bam";
						rm ${OUT_BAM}
						rm ${BAM_TMP}/${RUN}_${INDEX}.minimap2.sorted.bam
						
						echo -e `date +"%F\t%T\t"`"[RENAME] to ${BAM_TMP}/${RUN}.minimap2.sorted.bam";
						mv ${BAM_TMP}/${RUN}.merged.minimap2.sorted.bam ${BAM_TMP}/${RUN}.minimap2.sorted.bam
					fi
				else
					echo -e `date +"%F\t%T\t"`"[RENAME] to ${BAM_TMP}/${RUN}.minimap2.sorted.bam";
					mv -n ${BAM_TMP}/${RUN}_${INDEX}.minimap2.sorted.bam ${BAM_TMP}/${RUN}.minimap2.sorted.bam
				fi 
				
				# Index result
				echo -e `date +"%F\t%T\t"`"[INDEX] ${BAM_TMP}/${RUN}.minimap2.sorted.bam";
				samtools index -@ ${THREAD_LESS} "${BAM_TMP}/${RUN}.minimap2.sorted.bam"
				
				# Run modbam2bed 
				echo -e `date +"%F\t%T\t"`"[START] modbam2bed ${BAM_TMP}/${RUN}.minimap2.sorted.bam";
				modbam2bed -e -m 5mC -t ${THREAD} -p "${BAM_TMP}/${RUN}" --cpg --chg --chh ${FASTA} ${BAM_TMP}/${RUN}.minimap2.sorted.bam
				SUCCESS_MODBAM2BED=$?
				echo -e `date +"%F\t%T\t"`"[END] modbam2bed ${BAM_TMP}/${RUN}.minimap2.sorted.bam";
				if [ ${SUCCESS_MODBAM2BED} -eq 0 ]; then 
					# Rmove from hexagon if merge happened
					if [ ${SUCCESS_MERGE} -eq 0 ]; then 
						echo -e `date +"%F\t%T\t"`"[REMOVE] modbam2bed ${BAM_BACKUP}/${RUN}";
						rm ${BAM_BACKUP}/${RUN}.*.bed.gz
					fi
					
					echo -e `date +"%F\t%T\t"`"[GZIP_MODBAM] ${BAM_TMP}/${RUN}";
					pigz -p ${THREAD} ${BAM_TMP}/${RUN}.chg.bed
					pigz -p ${THREAD} ${BAM_TMP}/${RUN}.chh.bed 
					pigz -p ${THREAD} ${BAM_TMP}/${RUN}.cpg.bed
				fi

				# Run mosdepth
				echo -e `date +"%F\t%T\t"`"[START] mosdepth ${BAM_TMP}/${RUN}.minimap2.sorted.bam";
				mosdepth -t ${THREAD} --no-per-base "${BAM_TMP}/${RUN}" "${BAM_TMP}/${RUN}.minimap2.sorted.bam"
				if [ ${SUCCESS_MERGE} -eq 0 ]; then 
					echo -e `date +"%F\t%T\t"`"[REMOVE] mosdepth ${BAM_BACKUP}/${RUN}";
					rm "${BAM_BACKUP}/${RUN}.mosdepth.*"
				fi
				
				# Update MUX_ID Status -> "MODBAM2BED"
				UPDATE='db.ARG_methylation_MUX.updateOne( {"_id": '"${DOC_ID}"' } , { $set:{"Status": "MODBAM2BED"} } ) ; '
				echo "UpdateSql: $UPDATE" >> ${DB_LOG}
				result=$(echo "$UPDATE" | $MONGO_EXEC -u promappuser -p "q^a#W3LTjq?t^kwB" 54.251.49.171:8080/promethion --authenticationDatabase promethion --quiet)
				echo "outcome: $result" >> ${DB_LOG}
				
				# Rsync to hexagon 
				echo -e `date +"%F\t%T\t"`"[START] rsync to ${BAM_BACKUP}";
				rsync -avu --remove-source-files ${BAM_TMP} ${BAM_BACKUP}
				echo -e `date +"%F\t%T\t"`"[END] rsync to ${BAM_BACKUP}";
				
				# Remove lock
				rm -rf -- "$lockdir"; echo "remove lock2"
			fi
		fi 
		trap 'rm -rf -- "$lockdir"; echo "remove lock2"' 0    # remove directory when script finishes
		i+=1

	else
		printf >&2 "Run $i cannot acquire lock, giving up on %s\n $lockdir"
		exit 0

	fi
done

