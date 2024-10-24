ROCKSDB_PATH=$HOME/EZC/rocksdb
LOG_PATH=$HOME/EZC/log_
RAW_ZNS=/nvme0n1
RAW_ZNS_PATH=/sys/block/${RAW_ZNS}/queue/scheduler
RESULT_DIR_PATH=$HOME/EZC/data
ZNS_SCHEDULER=/sys/block/nvme0n1/queue/scheduler

## zc_scheme
GREEDY=0
CBZC1=1
CBZC2=2
CBZC3=3

## Dataset
#MED=12582912 # 12
#MED=83886080 #80<<20 # 80MB*1024KB=80GB
#MED=73400320 #70GB
#MED=9437184 #9GB
MED=67108864 #64GB

##Tuning Point
T=80
T_COMPACTION=4
T_FLUSH=4
T_SUBCOMPACTION=8
ZC_KICKS=20
UNTIL=20

# sudo rm -rf ${RESULT_DIR_PATH}
# mkdir ${RESULT_DIR_PATH}


for zc_scheme_name in GREEDY CBZC1 CBZC2 CBZC3; do
    zc_scheme=${!zc_scheme_name}

    for run in 1 2 3; do
        if [ ! -d ${RESULT_DIR_PATH} ]; then
            echo "NO ${RESULT_DIR_PATH}"
            mkdir -p ${RESULT_DIR_PATH}
        fi

        echo "mq-deadline" | sudo tee ${ZNS_SCHEDULER}
        sudo rm -rf ${LOG_PATH}
        mkdir ${LOG_PATH}
        sudo ${ROCKSDB_PATH}/plugin/zenfs/util/zenfs mkfs --force --enable_gc --zbd=${RAW_ZNS} --aux_path=${LOG_PATH}
        sleep 3

        result_file=${RESULT_DIR_PATH}/result_${zc_scheme_name}_run_${run}
        echo $result_file
        
        sudo ${ROCKSDB_PATH}/db_bench \
            -num=$MED \
            -benchmarks="fillrandom,stats" \
            --fs_uri=zenfs://dev:nvme0n1 \
            -statistics \
            -value_size=1024 \
            -file_opening_threads=4 \
            -max_background_compactions=${T_COMPACTION} \
            -max_background_flushes=${T_FLUSH} \
            -subcompactions=${T_SUBCOMPACTION} \
            -histogram \
            -tuning_point=$T \
            -reset_scheme=0 \
            -partial_reset_scheme=1 \
            -zc=${ZC_KICKS} \
            -until=${UNTIL} \
            -allocation_scheme=0 \
            -zc_scheme=${zc_scheme} \
            -compaction_scheme=0 \
            -input_aware_scheme=0 \
            -max_compaction_kick=0 \
            > ${result_file}

        echo "Result saved to ${result_file}"
    done
    sleep 5
done



# gdb 명령어를 작성할 임시 스크립트 생성
#echo "run -num=$MED -benchmarks=\"fillrandom,stats\" --fs_uri=zenfs://dev:nvme0n1 -statistics -value_size=1024 \
#-file_opening_threads=4 -max_background_flushes=4 -max_background_compactions=4 -histogram > ${RESULT_DIR_PATH}/tmp" > gdb_commands.txt

# gdb 실행
#sudo gdb -x gdb_commands.txt --args ${ROCKSDB_PATH}/db_bench

# gdb 명령어 스크립트 삭제
#rm gdb_commands.txt
