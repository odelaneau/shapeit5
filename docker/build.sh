#Clean up
rm ../static_bins/*
rm resources/*


#Compile phase 1
cd ../phase_common/
make clean
make -j static_exe
cp bin/phase_common_static ../static_bins/.

#Compile phase 2
cd ../phase_rare/
make clean
make -j static_exe
cp bin/phase_rare_static ../static_bins/.

#Compile switch
cd ../switch/
make clean
make -j static_exe
cp bin/switch_static ../static_bins/.

#Compile ligate
cd ../ligate/
make clean
make -j static_exe
cp bin/ligate_static ../static_bins/.

#Compile xcftools
cd ../xcftools/
make clean
make -j static_exe
cp bin/xcftools_static ../static_bins/.

#Buld docker image
LAB=shapeit5_$(git log -1 --format=%cd --date=short)\_$(git rev-parse --short HEAD)

cd ../docker/
mkdir -p resources
cp ../static_bins/SHAPEIT5* resources/.

docker build -t $LAB -f Dockerfile .
docker save $LAB | gzip -c > $LAB\.tar.gz
