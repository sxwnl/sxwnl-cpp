# Sxwnl-cpp  
寿星天文历(万年历)c++版  
Sxwnl for cpp version
# Build
## 1.Use cmake  
```bash
# build or update bin
./autobuild.sh

# clean obj files and build
./rebuild.sh

```
## 2.Use make  
```bash
# clean obj files
make clean
# for windows
make cleanw

### Build test0
make test0
### Build test1
make test1
### build all
make

```

# Run test  
```bash
# cmake build
./build/test0
./build/test1
# make build
./test0
./test1
```

# Bugs
TODO

# 特别说明
本次移植大多数工作在手机上完成，使用[C4droid](https://blog.qaiu.top/archives/c4droid)开发工具，所以代码格式不太规范，但不影响使用，经过了后续多次迭代精度已经和js版完全一致，计算效率远高于js版。
