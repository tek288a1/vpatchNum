#!/bin/bash
echo "#### The beginning ############################"
date
echo "###############################################"
echo ""
matlab -nodesktop -nodisplay -r "calcA_fnc(256,0.95,1,2048); exit" | tee diary3.log 
echo ""
echo "#### The end ##################################"
date
echo "###############################################"
