
#    Just run the application
#Mesh information can be provided to c3po using the setting/mesh.json file
#The name of the *.csv file need to be specified in the system/c3po.json file
#The example is run using  4 cores, the csv interface will reconstruct the mesh from the information provided by the user
#and the domain is divided in orizontal slices
$C3PO_SRC_DIR/interface_CSV/c3po_csv
#mpirun -np 8 $C3PO_SRC_DIR/interface_CSV/c3po_csv
./Allclean.sh

