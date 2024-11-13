# LongitudinalProf-MSMS-Method
## Preamble

This repository contains the code for the paper titled "Longitudinal Profiling of MS/MS Fragments Improves Accuracy of Metabolite Identification". The project focuses on developing and validating methods for accurate metabolite identification using longitudinal profiling of MS/MS fragments.


## Code Availability
The code used in this study have been made accessible in this repository to ensure transparency and reproducibility of the findings. The following files are included:

1. 001-008.ipynb: This series of files includes key data processing steps, ranging from MS1 and MS/MS spectrum extraction, spectrum filtering, invalid fragment removal, fragment rationality assessment, CQMU database construction, isomer similarity calculations, to metabolite identification. These processes encompass data preprocessing, optimization, and selection, establishing a robust data foundation for subsequent analysis and metabolite identification.
2. mzbatch.py: This file contains a library of commonly used methods for mass spectrometry data processing, offering modules for file operations, data validation, data generation, data calculation, and data extraction. This file is designed for quick access to essential functions within Python-based mass spectrometry data analysis projects, facilitating efficient data management and analysis.
3. requirements.txt: This file lists all Python packages and their versions needed for this project.

## Usage
To utilize the LongitudinalProf-MSMS-Method and replicate the research results, please follow these steps:

1. **Clone or Download the Repository**

   ```bash
   git clone https://github.com/yourusername/LongitudinalProf-MSMS-Method.git
   cd LongitudinalProf-MSMS-Method

2. **Install Required Dependencies**

   ```bash
   conda create -n longprof python=3.8.13
   conda activate longprof
   pip install -r requirements.txt

3. **Run the Code**
Run every cell of "001-008.ipynb" in sequence in Jupyter Notebook. Make sure that the data set is stored in the same path as specified in the script.  

   
## Citation
The citation of this work would be greatly appreciated if you find it valuable or choose to build upon it.


## Contact
For any inquiries or questions regarding this research or its code, please feel free to contact Jun-Yan Liu at jyliu@cqmu.edu.cn.


## Copyright License
Permission is hereby granted, free of charge, to any person obtaining a copy of this project and associated documentation files (the "project"), to deal in the project without restriction, including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the project, and to permit persons to whom the project is furnished to do so, subject to the following conditions:

The data of this study, including the dataset and associated findings, are intended for experimental reference only and should not be directly interpreted or used to guide human medicine. The authors do not assume any responsibility or liability for any consequences arising from the use or application of the data. The project is provided "as is", without warranty of any kind, express or implied, including but not limited to the warranties of merchantability, fitness for a particular purpose, and non-infringement. in no event shall the authors or copyright holders be liable for any claim, damages, or other liability, whether in an action of contract, tort, or otherwise, arising from, out of, or in connection with the project or the use or other dealings in the project.

The above copyright notice and this permission notice shall be included in all copies or substantial portions of the project.
