ROCKSDB_PATH=$HOME/EZC/rocksdb
ZENFS_PATH=$HOME/EZC/rocksdb/plugin/zenfs

# cd ${ZENFS_PATH} && git pull
git pull https://minsky922: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@github.com/minsky922/EZC

cd ${ROCKSDB_PATH} && sudo DEBUG_LEVEL=0 ROCKSDB_PLUGINS=zenfs make -j16 db_bench install

if [ $? -ne 0 ]; then
    exit
fi
cd ${ROCKSDB_PATH}/plugin/zenfs/util && make clean && make
