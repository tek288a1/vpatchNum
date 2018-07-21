#!/bin/bash
echo "#### The beginning ############################"
date
echo "###############################################"
echo ""
matlab -nodesktop -nodisplay -r "calcA_fnc(256, 0.90, 1, 1024); exit" | tee diary2.log 
echo ""
echo "#### The end ##################################"
date
echo "###############################################"
