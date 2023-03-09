#RELEASE
RLS=v1.0.0

#Clean up
mkdir -p ../static_bins/
rm ../static_bins/SHAPEIT5_*
rm resources/SHAPEIT5_*


#Compile phase 1
cd ../phase_common/
make clean
make -j static_exe_robin_desktop
cp bin/SHAPEIT5_phase_common_static_${RLS} ../static_bins/.

#Compile phase 2
cd ../phase_rare/
make clean
make -j static_exe_robin_desktop
cp bin/SHAPEIT5_phase_rare_static_${RLS} ../static_bins/.

#Compile switch
cd ../switch/
make clean
make -j static_exe_robin_desktop
cp bin/SHAPEIT5_switch_static_${RLS} ../static_bins/.

#Compile ligate
cd ../ligate/
make clean
make -j static_exe_robin_desktop
cp bin/SHAPEIT5_ligate_static_${RLS} ../static_bins/.

#Compile convert
cd ../convert/
make clean
make -j static_exe_robin_desktop
cp bin/SHAPEIT5_convert_static_${RLS} ../static_bins/.

#Buld docker image
LAB=shapeit5_${RLS}

cd ../docker/
mkdir -p resources
cp ../static_bins/SHAPEIT5* resources/.

sudo docker build -t $LAB -f Dockerfile .
sudo docker save $LAB | gzip -c > $LAB\.tar.gz
