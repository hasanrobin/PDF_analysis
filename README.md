
# PDF_analysis

#!/bin/bash
#-------------We can make a bash script with the commands mentioned below -------------------------#

# Loading modules and activating conda environment 
module load anaconda/3.7
source /etc/profile.d/modules.sh
source activate hist

###### 1---- Instruction for PDF fit computation using the file geofit_zenedo.py -------
"""
Command format (example): python file_name.py in_year  end_year in_month  out_month var_name PDF_name
and, use the target PDF_name from the list: "skew", "weib", "expweib", "laplace"
"""
###---> Option1 : An example bash script for running  all  analyses(variables) together inside the ZEUS server 
for i in T2M D2M U10M V10M;
do
    bsub -q p_short -P R000 -R "span[ptile=36]" python geofit_zenedo.py 2006 2020 1 12 $i skew; 
    #echo $i    
done 

###---> Option2: An example of single variable analysis from the CLI/terminal using the below command
#  var_name list = "T2M", "D2M", "MSl", "U10M", "V10M"  

bsub -q p_short -P R000 -R  python geofit_zenedo.py 2006 2020 1 12 T2M skew; 


###### 2---- Instruction for PDF moments computation using the file moments_pdf.py -------
"""
Command format (example):  python file_name.py var_name
"""
###---> Option 1: Using the BSUB command , and for available variables (T2M, D2M, U10M, V10M, MSL)
bsub -q p_short -P R000 python moments_pdf.py T2M   # Change/input target variable name only

###---> Option 2: Under the CLI/terminal interface 
python moments_pdf.py T2M                            # Change/input target variable name only
