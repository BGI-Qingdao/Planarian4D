#awk '{print $2}' ../input_data/Aorder.txt > new_gene.txt
bash gen_data.sh $1; 
python draw.line.py $1;
python draw.heatmap.py $1;
#python hcluster.py
python sort_hcluster.py $1;

