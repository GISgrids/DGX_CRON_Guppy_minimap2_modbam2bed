#!/bin/bash
#Purpose = sync all other miscellaneous files from Promethion machine to mnt/seq network folder 
#Start
# run this and check error and email error
# */30 * * * * run-one rsync -avu --include="*/" --include="*.fast5" --include="*.fastq.gz" --exclude="*" \
# /data/TEST_CRON/ /mnt/seq/gridion/MinION_testing_data/ >> /data/TEST_CRON/rsync.log 2>&1

#Variables
SOURCE="/data/CRON_HEXAGON"
PROC="${SOURCE}/COMPLETE_PROJECT.file"
LOGFILE="${SOURCE}/rsync_COMPLETED.log"
ERROR="${SOURCE}/rsync_COMPLETED.ERROR"
DESTINATION="/mnt/hexagon-scratch/Instrument_Data_Share/JJ_ARG_Project_methylation"
TEST="/test_run_low_data_volume"

# MongoDb 
MONGO_EXEC="/data/TEST_CRON/mongo/bin/mongo"
DB_LOG="/data/CRON_HEXAGON/mongo.log"

# Email
SENDER="GIS.Nanopore@gmail.com"
RECIPIENT="li_zhihui@gis.a-star.edu.sg"
SUBJECT=$(echo "Promethion Project sync FAILED")
BODY="Promethion Error check log file $LOGFILE"

find ${SOURCE}/data -iname "final_summary_*.txt" | while read SUMMARY_FILE ; do
  	# echo "${SUMMARY_FILE} "
	IFS='/' read -ra ARRAY <<< "${SUMMARY_FILE}"
	RUN=${ARRAY[-4]}
	ELM=${ARRAY[-3]}
	INDEX=${ARRAY[-2]}
	FILE=${ARRAY[-1]}
	MUX_FULL_PATH=${RUN}/${ELM}/${INDEX}
	
	# Summary file's fast5 count
	FAST5_COUNT_SUMMARY=$(grep "fast5_files_in_final_dest" ${SUMMARY_FILE} | cut -d "=" -f 2 - )
	
	echo -e `date +"%F\t%T\t"`"Processing: ${FILE} ( RUN: ${RUN}, ELM: ${ELM}, INDEX: ${INDEX}) " >>${ERROR}
	
	# rsync ONLY 1 MUX folder to hexagon
	SYNC_SOURCE=${SOURCE}/data/${MUX_FULL_PATH}/
	SYNC_DEST=${DESTINATION}/data/${MUX_FULL_PATH}/
	run-one rsync -avu --remove-source-files --include="*/" --exclude="*.fast5" --exclude="*.fastq.gz" ${SYNC_SOURCE} ${SYNC_DEST} >> ${LOGFILE} 2>>${ERROR}
	SUCCESS=$?
	
	#Email Notification
	if [ "${SUCCESS}" = "0" ]; 
	then
		# Count number of fast5 in hexagon
		FAST5_COUNT=$(find ${SYNC_DEST} -iname "*.fast5" | wc -l)
		
		# Status of MUX 
		# 1) To_Basecall : DGX/Parabrick to basecall 
		# 2) Fast5_not_equal : Summary count vs hexagon counted fast5 do not tally, needs to troubleshoot
		STATUS=""
		if [ ${FAST5_COUNT} -eq ${FAST5_COUNT_SUMMARY} ]; then
			STATUS="To_Basecall"
		else
			STATUS="Fast5_not_equal"
		fi
		
		EVAL='db.ARG_methylation_MUX.insert({"mux_id": "'"${RUN}"'" , "elm_id": "'"${ELM}"'" , "index_id": "'"${INDEX}"'" , "mux_full_path": "'"${MUX_FULL_PATH}"'" , "fast5Count_summary": "'"${FAST5_COUNT_SUMMARY}"'" , "fast5Count_hexagon": "'"${FAST5_COUNT}"'" , "Status": "'"${STATUS}"'" });'
		echo "insertSql: $EVAL" >> ${DB_LOG}
		result=$(echo $EVAL | ${MONGO_EXEC} -u promappuser -p "q^a#W3LTjq?t^kwB" 54.251.49.171:8080/promethion --authenticationDatabase promethion --quiet)
		echo "outcome: $result" >> ${DB_LOG}
		echo -e `date +"%F\t%T\t"`"COMPLETED: ${FILE} ( RUN: ${RUN}, ELM: ${ELM}, INDEX: ${INDEX}) " >>${ERROR}
	else
		#echo -e `date +"%F\t%T\t"`"FAILED: ${FILE} ( RUN: ${RUN}, ELM: ${ELM}, INDEX: ${INDEX}) " >>${ERROR} 
		#echo "$BODY" | mail -s "$SUBJECT" -r "$SENDER" "$RECIPIENT" 
		#echo "Email Sent."
		exit 1
	fi
done
#END