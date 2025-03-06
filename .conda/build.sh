set -x
wd=$(pwd)
#echo $PREFIX >&2
#ls $wd >&2

# preparing INCLUDE files
mkdir -p $PREFIX/include
cp -r $wd/cereal/include/cereal $PREFIX/include
cp -r $wd/Eigen $PREFIX/include

# build binaries
make
make install PREFIX=$PREFIX
#ls $PREFIX/bin >&2
