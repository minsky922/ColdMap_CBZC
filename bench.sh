ROCKSDB_PATH=$HOME/EZC/rocksdb
LOG_PATH=$HOME/EZC/log_
RAW_ZNS=/nvme0n1
RAW_ZNS_PATH=/sys/block/${RAW_ZNS}/queue/scheduler
RESULT_DIR_PATH=$HOME/EZC/data

## Algorithm
EAGER=0
LAZY=1

## Dataset
#MED=12582912 # 12
#MED=83886080 #80<<20 # 80MB*1024KB=80GB
#MED=73400320 #70GB
#MED=9437184 #9GB
MED=67108864 #64GB

##Tuning Point
T=80


sudo rm -rf ${RESULT_DIR_PATH}
mkdir ${RESULT_DIR_PATH}

if [ ! -d ${RESULT_DIR_PATH} ]
then
    echo "NO ${RESULT_DIR_PATH}"
    mkdir -p ${RESULT_DIR_PATH}
fi

sudo ${ROCKSDB_PATH}/db_bench -num=$MED -benchmarks="fillrandom,stats" --fs_uri=zenfs://dev:nvme0n1 -statistics -value_size=1024 \
-file_opening_threads=4 -max_background_flushes=4 -max_background_compactions=4 -histogram -reset_scheme=$ALGORITHM -reset_at_foreground=true -tuning_point=$T > ${RESULT_DIR_PATH}/tmp

# gdb 명령어를 작성할 임시 스크립트 생성
#echo "run -num=$MED -benchmarks=\"fillrandom,stats\" --fs_uri=zenfs://dev:nvme0n1 -statistics -value_size=1024 \
#-file_opening_threads=4 -max_background_flushes=4 -max_background_compactions=4 -histogram > ${RESULT_DIR_PATH}/tmp" > gdb_commands.txt

# gdb 실행
#sudo gdb -x gdb_commands.txt --args ${ROCKSDB_PATH}/db_bench

# gdb 명령어 스크립트 삭제
#rm gdb_commands.txt
