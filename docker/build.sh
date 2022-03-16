#Compile phase 1
cd ../phase1/
make clean
make -j static_exe
cp bin/SHAPEIT5_phase1_static ../static_bins/.

cd ../phase2/
make clean
make -j static_exe
cp bin/SHAPEIT5_phase2_static ../static_bins/.

cd ../switch/
make clean
make -j static_exe
cp bin/SHAPEIT5_switch_static ../static_bins/.

cd ../ligate/
make clean
make -j static_exe
cp bin/SHAPEIT5_ligate_static ../static_bins/.

cd ../docker/
docker build -t shapeit5:0.0.1 -f Dockerfile .
docker save shapeit5:0.0.1 | gzip -c > shapeit5:0.0.1.tar.gz
