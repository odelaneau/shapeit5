#Compile phase 1
cd ../phase1/
make clean
make -j static_exe
cp bin/SHAPEIT5_phase1_static ../static_bins/.

#Compile phase 2
cd ../phase2/
make clean
make -j static_exe
cp bin/SHAPEIT5_phase2_static ../static_bins/.

#Compile switch
cd ../switch/
make clean
make -j static_exe
cp bin/SHAPEIT5_switch_static ../static_bins/.

#Compile ligate
cd ../ligate/
make clean
make -j static_exe
cp bin/SHAPEIT5_ligate_static ../static_bins/.

#Buld docker image
cd ../docker/
cp ../static_bins/SHAPEIT5* ressources/.
docker build -t shapeit5_0.0.1 -f Dockerfile .
docker save shapeit5_0.0.1 | gzip -c > shapeit5_0.0.1.tar.gz
