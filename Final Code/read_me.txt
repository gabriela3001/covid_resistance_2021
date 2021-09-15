In this folder is the code necessary to reproduce results from the paper "Evolution of resistance to COVID-19
vaccination with dynamic lockdown" by Gabriela Lobinska, Ady Pauzner, Arne Traulsen, Yitzhak Pilpel, Martin A.
Nowak (2021). 

In scripts_basic_model/basic_model and in scripts_model_extensions/* are the scripts needed for each model 
extension simulation. To run simulations on a computer cluster, modify the script to contain your desired
target directory for the results, and run the wrapper provided. 

To parse the results, use the provided Jupyter notebooks in parsing_model_extensions/