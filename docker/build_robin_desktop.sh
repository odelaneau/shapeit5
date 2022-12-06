#Clean up
rm ../static_bins/SHAPEIT5_*
rm resources/SHAPEIT5_*


#Compile phase 1
cd ../phase_common/
make clean
make -j static_exe_robin_desktop
cp bin/SHAPEIT5_phase_common_static ../static_bins/.

#Compile phase 2
cd ../phase_rare/
make clean
make -j static_exe_robin_desktop
cp bin/SHAPEIT5_phase_rare_static ../static_bins/.

#Compile switch
cd ../switch/
make clean
make -j static_exe_robin_desktop
cp bin/SHAPEIT5_switch_static ../static_bins/.

#Compile ligate
cd ../ligate/
make clean
make -j static_exe_robin_desktop
cp bin/SHAPEIT5_ligate_static ../static_bins/.

#Compile convert
cd ../convert/
make clean
make -j static_exe_robin_desktop
cp bin/SHAPEIT5_convert_static ../static_bins/.

#Buld docker image
LAB=shapeit5_$(git log -1 --format=%cd --date=short)\_$(git rev-parse --short HEAD)

cd ../docker/
mkdir -p resources
cp ../static_bins/SHAPEIT5* resources/.

docker build -t $LAB -f Dockerfile .
docker save $LAB | gzip -c > $LAB\.tar.gz
