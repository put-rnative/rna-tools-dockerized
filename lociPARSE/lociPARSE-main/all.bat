@echo off
echo "========Running lociPARSE========="
echo ""
echo "========Generating Features======="
python Scripts/Feature.py
echo "Done"

echo "=========Running Inference========"
python Scripts/prediction.py
echo ""
echo "=======Inference complete.========"