# EZC

지금까지 **EZC** 레포지토리에 **rocksdb**와 **zenfs**를 서브트리로 추가하고 푸시한 과정을 요약하면 다음과 같습니다.

### 1. **EZC 리포지토리에 rocksdb 추가**
- 먼저 **rocksdb**를 **EZC** 리포지토리에 서브트리로 추가.
- `rocksdb`를 원격 저장소로 푸시:
  ```bash
  git subtree add --prefix=rocksdb https://github.com/minsky922/myrocksdb.git kw --squash
  git push origin main
  ```

### 2. **rocksdb 내 plugin/zenfs 디렉토리에 zenfs 추가**
- **rocksdb** 내 `plugin/zenfs` 디렉토리에 **zenfs**를 서브트리로 추가:
  ```bash
  git subtree add --prefix=rocksdb/plugin/zenfs https://github.com/minsky922/myzenfs.git master --squash
  ```

### 3. **변경 사항 커밋 및 원격 저장소로 푸시**
- **zenfs**가 추가된 변경 사항을 커밋하고, 이를 원격 저장소로 푸시:
  ```bash
  git commit -m "Add zenfs as a subtree in rocksdb/plugin/zenfs"
  git push origin main
  ```

### 최종 결과:
- **EZC** 리포지토리에는 **rocksdb**와 **zenfs**가 각각 서브트리로 추가되어 있으며, **rocksdb** 안에 **plugin/zenfs**가 포함된 구조로 관리되고 있습니다.
