#!/bin/bash

#----------------------------------------------------
# FOLDERS VARIABLE
# mounted network folder
SOURCE="/mnt/hexagon-scratch/Instrument_Data_Share/JJ_ARG_Project_methylation"
# DGX local processing folder
DESTINATION="/raid/scratch/lizh/methyl_test"
# Backup locale
S3_dir="/mnt/hexagon-scratch/lizh/methy_test/guppy6_result"
# Switch to local drive temporarily due to hexagon no space
# S3_dir="/raid/scratch/lizh/methyl_test/basecall/"
# Log files
LOGFILE="${DESTINATION}/Guppy6_cron.log"
ERROR="${DESTINATION}/Guppy6_cron.ERROR"

# MongoDb 
MONGO_EXEC="/raid/scratch/lizh/mongo/bin/mongo"
DB_LOG="${DESTINATION}/mongo.log"

DOCKER="679ddb617c63"
GUPPY_VERSION="guppy_6.2.1"

# Max Guppy run in a row
MAX_RUN=10
#-----------------------------------------------------

lockdir=${DESTINATION}/tmp/guppy6_cron.lock
declare -i i=0
while [ $i -lt ${MAX_RUN} ]
do 
	echo "Run $i" | gawk '{ print strftime("[%Y-%m-%d %H:%M:%S]"), $0 }' >>${ERROR}
	if mkdir -- "$lockdir" ; then
		printf >&2 'successfully acquired lock\n' >>${ERROR}
		
		RESUME=""
		# Check for MUX_ID with Status "To_Basecall" from mongoDB
		EVAL='db.ARG_methylation_MUX.findOne({"Status": "RERUN_BC"})._id ; '
		echo "querySql: $EVAL" >> ${DB_LOG}
		DOC_ID=$(echo "$EVAL" | $MONGO_EXEC -u promappuser -p "q^a#W3LTjq?t^kwB" 54.251.49.171:8080/promethion --authenticationDatabase promethion --quiet)
		if [[ ${DOC_ID} =~ 630* ]]; then
			RESUME="--resume"
		else
			# Check for MUX_ID with Status "To_Basecall" from mongoDB
			EVAL='db.ARG_methylation_MUX.findOne({"Status": "To_Basecall"})._id ; '
			echo "querySql: $EVAL" >> ${DB_LOG}
			DOC_ID=$(echo "$EVAL" | $MONGO_EXEC -u promappuser -p "q^a#W3LTjq?t^kwB" 54.251.49.171:8080/promethion --authenticationDatabase promethion --quiet)
		fi 
		
		EVAL='db.ARG_methylation_MUX.findOne('"${DOC_ID}"').mux_full_path ; '
		echo "querySql: $EVAL" >> ${DB_LOG}
		MUX_FULL_PATH=$(echo "$EVAL" | $MONGO_EXEC -u promappuser -p "q^a#W3LTjq?t^kwB" 54.251.49.171:8080/promethion --authenticationDatabase promethion --quiet)
		# echo "${DOC_ID} , ${MUX_FULL_PATH}"
		SUCCESS=$?

		if [ "${SUCCESS}" = "0" ]; then
			# eg ON002-DNA-R00141/WHB12183-T1/20210208_0141_1G_PAG51949_51a655a0/fast5_fail/PAG51949_fail_7683c71e_99.fast5
			IFS='/' read -ra ARRAY <<< "${MUX_FULL_PATH}"
			RUN=${ARRAY[-3]}
			ELM=${ARRAY[-2]}
			INDEX=${ARRAY[-1]}

			echo "MUX Full Path is ${RUN}/${ELM}/${INDEX}" | gawk '{ print strftime("[%Y-%m-%d %H:%M:%S]"), $0 }' >>${ERROR}

			# Update MUX_ID Status -> "Basecalling"
			UPDATE='db.ARG_methylation_MUX.updateOne( {"_id": '"${DOC_ID}"' } , { $set:{"Status": "Basecalling"} } ) ; '
			echo "UpdateSql: $UPDATE" >> ${DB_LOG}
			result=$(echo "$UPDATE" | $MONGO_EXEC -u promappuser -p "q^a#W3LTjq?t^kwB" 54.251.49.171:8080/promethion --authenticationDatabase promethion --quiet)
			echo "outcome: $result" >> ${DB_LOG}
			
			cd ${DESTINATION}
			#[ ! -d basecall ] && mkdir -p basecall && echo "Created basecall"
			
			# Sync fast5 from Hexagon to scratch
			SYNC_SOURCE=${SOURCE}/data/${RUN}/${ELM}/${INDEX}/
			SYNC_DEST=${DESTINATION}/data/${RUN}/${ELM}/${INDEX}/
			echo "Processing ${RUN}" | gawk '{ print strftime("[%Y-%m-%d %H:%M:%S]"), $0 }' >>${ERROR}
			echo "START Sync data from ${SYNC_SOURCE} to ${SYNC_DEST}" | gawk '{ print strftime("[%Y-%m-%d %H:%M:%S]"), $0 }' >>${ERROR}
			# NOT --dry-run
			[ ! -d ${SYNC_DEST} ] && mkdir -p ${SYNC_DEST} && echo "Created ${SYNC_DEST}" >>${ERROR}
			run-one rsync -avu -q --include="*/" --include="*.fast5" --exclude="*" ${SYNC_SOURCE} ${SYNC_DEST} >> ${LOGFILE} 2>>${ERROR}
			echo "END Sync data " | gawk '{ print strftime("[%Y-%m-%d %H:%M:%S]"), $0 }' >>${ERROR}

			BASECALL_FOLDER="/mnt/scratch/data/${RUN}/${ELM}/${INDEX}"
			BASECALL_OUTPUT="/mnt/scratch/data/${RUN}/${ELM}/${INDEX}/Guppy_6.2.1"
				
			echo "Running ${GUPPY_VERSION} on ${BASECALL_FOLDER}" | gawk '{ print strftime("[%Y-%m-%d %H:%M:%S]"), $0 }' >>${ERROR}
			/usr/bin/docker exec ${DOCKER} guppy_basecaller -i ${BASECALL_FOLDER} -s ${BASECALL_OUTPUT} \
			-c dna_r10.4.1_e8.2_400bps_modbases_5mc_cg_hac_prom.cfg \
			-x "cuda:0,1,2,4,7:8G cuda:3,5,6:20G" \
			--num_callers 78 \
			--disable_pings --compress_fastq --recursive --min_qscore 7 \
			--bam_out --index ${RESUME} \
			>> ${LOGFILE} 2>>${ERROR}
			SUCCESS=$?
			
			# --chunks_per_runner 1024 --gpu_runners_per_device 36 --chunk_size 2000 \
			# "cuda:0:18G cuda:1,2,3,4,5,7:20G cuda:6:10G"

			#SUCCESS=0
			if [[ ${SUCCESS} -eq 0 ]]; then
				# Remove Guppy lockfile
				echo "Guppy COMPLETED, removing file lock" | gawk '{ print strftime("[%Y-%m-%d %H:%M:%S]"), $0 }' >>${ERROR}
				rm -rf -- "$lockdir"
				
				# Sync basecall result back to hexagon				
				S3_source=${SYNC_DEST}/Guppy_6.2.1/
				S3_dest=${S3_dir}/${RUN}/${ELM}/${INDEX}/Guppy_6.2.1/
				[ ! -d ${S3_dest} ] && mkdir -p ${S3_dest} && echo "Created ${S3_dest}"
				echo "START Rsync to Hexagon : ${S3_dest}" | gawk '{ print strftime("[%Y-%m-%d %H:%M:%S]"), $0 }' >>${ERROR}
				#aws s3 sync ${S3_source} ${S3_dest} 
				run-one rsync -avu -q ${S3_source} ${S3_dest} >> ${LOGFILE} 2>>${ERROR}
				# run-one rsync -avu -q ${S3_source} ${SYNC_SOURCE} >> ${LOGFILE} 2>>${ERROR}
				echo "END Rsync " | gawk '{ print strftime("[%Y-%m-%d %H:%M:%S]"), $0 }' >>${ERROR}
				
				# Update MongoDB
				UPDATE='db.ARG_methylation_MUX.updateOne( {"_id": '"${DOC_ID}"' } , { $set:{"Status": "DONE_BASECALL"} } ) ; '
				echo "UpdateSql: $UPDATE" >> ${DB_LOG}
				result=$(echo "$UPDATE" | $MONGO_EXEC -u promappuser -p "q^a#W3LTjq?t^kwB" 54.251.49.171:8080/promethion --authenticationDatabase promethion --quiet)
				echo "outcome: $result" >> ${DB_LOG}
				
				# Remove folder
				rm -rf ${DESTINATION}/data/${RUN}/${ELM}/${INDEX}
				find ${DESTINATION}/data/${RUN} -depth -exec rmdir {} \;  
				
			else
				echo "Guppy FAILED for MUX ${RUN}" | gawk '{ print strftime("[%Y-%m-%d %H:%M:%S]"), $0 }' >>${ERROR}
				# Update MongoDB
				UPDATE='db.ARG_methylation_MUX.updateOne( {"_id": '"${DOC_ID}"' } , { $set:{"Status": "RERUN_BC"} } ) ; '
				echo "UpdateSql: $UPDATE" >> ${DB_LOG}
				result=$(echo "$UPDATE" | $MONGO_EXEC -u promappuser -p "q^a#W3LTjq?t^kwB" 54.251.49.171:8080/promethion --authenticationDatabase promethion --quiet)
				echo "outcome: $result" >> ${DB_LOG}
				exit 2
			fi
			
		else
			echo "MUX Full Path is empty" | gawk '{ print strftime("[%Y-%m-%d %H:%M:%S]"), $0 }' >>${ERROR}
			exit 1
		fi 
		trap 'rm -rf -- "$lockdir"; echo "remove lock2"' 0    # remove directory when script finishes
		i+=1
	else
		printf >&2 "Run $i cannot acquire lock, giving up on %s\n $lockdir"
		exit 0
	fi
done
