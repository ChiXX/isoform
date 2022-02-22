# Development and comparison of feature-selection pipelines for predicting patient outcome from RNA-seq data in breast cancer

This is a 60-credit degree project. 

The aims of this project:

* To predict patient outcomes based on gene and transcript expression.
* To explore and improve methods for feature selection which extract the informative genes/transcripts associated with breast cancer relapse.

With the aim to select informative features associated with breast cancer patient outcomes, a customized ensemble feature selection model based on penalized Cox proportional hazards (PH) model was developed and compared with the univariate method and single-run Cox PH model. This method can be employed to select features in unbalanced, ultra-high-dimensional, time-to-event data.

## Introduction

### Breast cancer

![https://gco.iarc.fr/today/online-analysis-pie](mdFigs/breast.png)

Breast cancer overtaking lung cancer became the most common cancer worldwide. In 2020, 2.26 million new cases and 685 000 deaths were caused by breast cancer, which was the fifth leading cause of cancer death.

### SCAN-B

[The Sweden Cancerome Analysis Network–Breast (SCAN-B) Initiative](https://pubmed.ncbi.nlm.nih.gov/25722745/)  was initially proposed in 2009 and began enrolling patients in 2010. Within SCAN-B, breast tumors from hospitals across a wide geography of Sweden are routinely being processed and RNA-sequenced. Up to now, more than 17 thousand patients have been enrolled in SCAN-B, and more than 11 thousand tumors have been RNA-sequenced, which provides a rich resource for breast cancer research. To our current knowledge, SCAN-B is the largest RNA-seq breast cancer study in the world!

### Survival analysis and time to event data

Survival analysis is a type of regression problem which tries to establish a connection between covariates and the time of an event. Here the definition of the start point is usually the date of primary treatment. But in SCAN-B, the start point is the date of diagnosis. As for the endpoint, according to the specific question, it can be **overall survival**  (OS) which includes all kinds of death events, **relapse-free survival** (RFS), which includes death and relapse event, and **recurrence-free interval** (RFi), which includes recurrence events only. 

A common issue in survival analysis is that the data is censored.

<img src="mdFigs/censored_data.png" alt="https://scikit-survival.readthedocs.io/en/stable/user_guide/understanding_predictions.html" style="zoom:50%;" />

The events for patients B and D are recorded. But for patients A, C, and E, the only available information is they are event-free up to their last follow-up. So they are censored.

## MATERIALS AND USAGE

### Materials

The breast cancer data is obtained form SCAN-B, which contains **2874 ER+/HER2- samples**, which are used to predict patient outcomes.

|                 | **RFi event** | **RFi event** | **RFi event** | **RFi event** |
| :-------------: | :-----------: | :-----------: | :-----------: | :-----------: |
|                 |     Train     |  Validation   |     Train     |  Validation   |
|   With event    |      116      |      27       |      342      |      95       |
|  Without event  |     1139      |      339      |     1957      |      480      |
| Censoring ratio |     90.8%     |     92.6%     |     85.1%     |     83.5%     |

### Usage

All the scripts are stored in `scripts` and managed with *Snakemake*.  It is recommend to use cluster to carry out the project.

How to start:

```bash
# A conda enviorment is available in '/env/ballgown'
# Unpack environment into directory `ballgown`
cd env
mkdir -p ballgown
tar -xzf ballgown.tar.gz -C ballgown
# Activate the environment.
$ source ballgown/bin/activate
```

Workflow design

* *de novo* assembly

![image-20220222180811657](mdFigs/workflow1.png)



