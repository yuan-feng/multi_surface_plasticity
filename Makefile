vm:
	g++ -o multi_ys_VM test_multi_ys_shear.cpp vonMises_multi_surface.cpp MultiYieldSurfaceMaterial.cpp -std=c++11 -Wall
	python test_vonMises_multi_shear.py
	
dp:
	g++ -o multi_ys_DP test_multi_ys_DP_shear.cpp DruckerPrager_multi_yield_surface.cpp MultiYieldSurfaceMaterial.cpp -std=c++11 -Wall
	python test_GUI_MYS_DP.py



rmc:
	g++ -o test_MYS_RMC_equal_shear test_MYS_RMC_equal_shear.cpp RoundedMohrCoulomb_multi_surface.cpp MultiYieldSurfaceMaterial.cpp -std=c++11 -Wall
	python test_GUI_MYS_RMC.py	

# compile: 
# 	g++ -o multi_ys_DP test_multi_ys_DP_shear.cpp DruckerPrager_multi_yield_surface.cpp -std=c++11 -Wall

# plot:
# 	python plot.py

# GUI:
# 	python test_GUI_MYS_DP.py
	
clean:
	-rm -f  multi_ys
	-rm -f  strain_stress.txt
	-rm -f  results.pdf