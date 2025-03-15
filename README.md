# Social_ALAN_Data
This repository contains data and analysis code for a study on social effects under Artificial Light at Night (ALAN).

## Repository Structure

Social_ALAN_Data/
├── Codes/
│   ├── Behavior_Code.R        # R script for behavioral data analysis
│   └── codes.R                # Additional R analysis scripts
├── Data/
│   ├── Raw/
│   │   ├── Behavior/         # Raw behavioral data files
│   │   │   ├── meta*.csv     # Metadata files
│   │   │   └── Raw Data*.zip # Raw behavioral measurements
│   │   ├── Gene_Expression/  # Gene expression measurement data
│   │   └── Melatonin/        # Melatonin measurement data
│   └── Social_ALAN.csv       # Processed dataset
└── README.md                 # Project documentation

```

## Description

This repository contains raw data and analysis scripts. The project includes both behavioral measurements, gene expression, and melatonin data.

## Data
- **Gene Expression Data**: Located in `Data/Raw/Gene_Expression/`, containing gene expression measurements
- **Behavior Data**: Located in `Data/Raw/Behavior/`, including raw measurements and metadata
- **Melatonin Data**: Located in `Data/Raw/Melatonin/`, containing experimental measurements
- **Processed Data**: `Social_ALAN.csv` contains the processed dataset

## Analysis Code

The analysis is performed using R scripts located in the `Codes/` directory:
- `Behavior_Code.R`: Main script for behavioral data analysis
- `codes.R`: Additional analysis routines