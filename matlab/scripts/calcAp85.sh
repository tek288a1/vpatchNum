#!/bin/bash
echo "#### The beginning ############################"
date
echo "###############################################"
echo ""
matlab -nodesktop -nodisplay -r "calcA_fnc(256,0.85,1,1024); exit" | tee diary1.log 
echo ""
echo "#### The end ##################################"
date
echo "###############################################"
