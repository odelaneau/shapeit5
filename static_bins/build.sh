#Compile phase 1
cd ../phase_common/
make clean
make -j static_exe
mv bin/SHAPEIT5_phase_common_static ../static_bins/.

#Compile phase 2
cd ../phase_rare/
make clean
make -j static_exe
mv bin/SHAPEIT5_phase_rare_static ../static_bins/.

#Compile switch
cd ../switch/
make clean
make -j static_exe
mv bin/SHAPEIT5_switch_static ../static_bins/.

#Compile ligate
cd ../ligate/
make clean
make -j static_exe
mv bin/SHAPEIT5_ligate_static ../static_bins/.
