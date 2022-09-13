#!/bin/bash
# Purpose = sync files from Promethion machine to hexagon-scratch 

# Variables
SOURCE="/data/CRON_HEXAGON"
DESTINATION="/mnt/hexagon-scratch/Instrument_Data_Share/JJ_ARG_Project_methylation"
LOGFILE="/data/CRON_HEXAGON/rsync.log"
ERROR="/data/CRON_HEXAGON/rsync.error.log"
SENDER="GIS.Nanopore@gmail.com"
RECIPIENT="li_zhihui@gis.a-star.edu.sg"
SUBJECT="Promethion Data sync FAILED"
BODY="Promethion Error check log file $LOGFILE"

# echo -e `date +"%F\t%T\t"`"[START] Rsync ${SOURCE} to ${DESTINATION}" >>${ERROR}
#Do the deed 
run-one rsync -avu --remove-source-files --include="*/" --include="*.fast5" --include="*.fastq.gz" --exclude="*" ${SOURCE}/data ${DESTINATION} >> $LOGFILE 2>>${ERROR}
SUCCESS=$?
echo $SUCCESS

#Email Notification
if [ "${SUCCESS}" = "0" ]; then
# echo -e `date +"%F\t%T\t"`"[SUCCESS] Rsync" >>${ERROR}
cd ${SOURCE}
grep ".fast5$" ${LOGFILE} | sort --unique > ${SOURCE}/FAST5.LIST
mv -f ${SOURCE}/FAST5.LIST ${DESTINATION} 
else
echo -e `date +"%F\t%T\t"`"[FAIL] Rsync \n" >>${ERROR}
echo "$BODY" | mail -s "$SUBJECT" -r "$SENDER" -a "${ERROR}" "$RECIPIENT"
echo "Email Sent."
exit 1
fi
#END