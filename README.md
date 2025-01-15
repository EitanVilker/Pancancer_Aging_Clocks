# Pancancer_Aging_Clocks

## Environment Setup

### Option 1: Using Anaconda
1. Install Anaconda: [Anaconda Installation Guide](https://www.anaconda.com/products/distribution#download-section).
2. Clone this repository:
   ```bash
   git clone https://github.com/EitanVilker/Pancancer_Aging_Clocks.git
   cd Pancancer_Aging_Clocks
   ```
3. Create the environment using the provided environment.yml file:
   ```bash
   conda env create -f environment.yml
   ```
5. Activate the environment
   ```bash
   conda activate r_env
   ```

### Option 2: Manual
1. Install R (version 4.4.0) from [CRAN](https://cran.r-project.org/).
2. Clone this repository:
   ```bash
   git clone https://github.com/EitanVilker/Pancancer_Aging_Clocks.git
   cd Pancancer_Aging_Clocks
   ```
3. Run the provided `requirements.R` script to install the necessary libraries:
   ```R
   source("requirements.R")
   ```
