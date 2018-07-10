m4 -P system/blockMeshDict.m4 > system/blockMeshDict 
blockMesh
rotateMesh "(0 0 1)" "(1 0 0)"
renumberMesh -overwrite
checkMesh
