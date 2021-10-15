g++ -std=c++11 -Wno-deprecated -Wno-unknown-pragmas -D_REENTRANT -DNDEBUG -DSC_MEMORY_ALIGNMENT=1 -DSC_LINUX -Wno-deprecated -fPIC -I ../../sc-310/supercollider/include/plugin_interface/ -I ../../sc-310/supercollider/include/common/ -I ../../sc-310/supercollider/include/server/ -O2 -c SpecTrans.cpp 

g++ -shared -g -o SpecTrans.so SpecTrans.o

cp SpecTrans.sc ~/.local/share/SuperCollider/Extensions/
cp SpecTrans.so ~/.local/share/SuperCollider/Extensions/
