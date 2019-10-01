---
title: "v3.4"
author: "Steven N. Hart, Ph.D"
date: '9/30/2019'
output:
  html_document:
    keep_md: true
---





Plot all predictions


Set some useful functions




Load all functional data










```
## , ,  = BRCA1
## 
##                         
##                          Damaging Neutral
##   Fernandes et al., 2019        0      16
##   Findlay et al., 2018         55     259
##   Guidugli et al., 2014         0       0
##   Guidugli et al., 2018         0       0
##   Hart et al, 2018              1       2
##   Lee et al., 2010             14      13
##   Lindor et al., 2012           0       0
##   Starita et al., 2015         11    1044
##   This paper                    0       0
##   Woods et al., 2016            3      36
## 
## , ,  = BRCA2
## 
##                         
##                          Damaging Neutral
##   Fernandes et al., 2019        0       0
##   Findlay et al., 2018          0       0
##   Guidugli et al., 2014         4       0
##   Guidugli et al., 2018        38      48
##   Hart et al, 2018             16      47
##   Lee et al., 2010              0       0
##   Lindor et al., 2012           9      18
##   Starita et al., 2015          0       0
##   This paper                    7      15
##   Woods et al., 2016            0       0
```

```
##           
##            BRCA1 BRCA2
##   Damaging    84    74
##   Neutral   1370   128
```

More custom functions



Define columns to TRAIN on


Run the training

```
## 
## H2O is not running yet, starting it now...
## 
## Note:  In case of errors look at the following log files:
##     /tmp/RtmpH7IIge/h2o_m087494_started_from_r.out
##     /tmp/RtmpH7IIge/h2o_m087494_started_from_r.err
## 
## 
## Starting H2O JVM and connecting: . Connection successful!
## 
## R is connected to the H2O cluster: 
##     H2O cluster uptime:         1 seconds 449 milliseconds 
##     H2O cluster timezone:       CST6CDT 
##     H2O data parsing timezone:  UTC 
##     H2O cluster version:        3.22.1.1 
##     H2O cluster version age:    9 months and 2 days !!! 
##     H2O cluster name:           H2O_started_from_R_m087494_whe707 
##     H2O cluster total nodes:    1 
##     H2O cluster total memory:   28.44 GB 
##     H2O cluster total cores:    80 
##     H2O cluster allowed cores:  20 
##     H2O cluster healthy:        TRUE 
##     H2O Connection ip:          localhost 
##     H2O Connection port:        54321 
##     H2O Connection proxy:       NA 
##     H2O Internal Security:      FALSE 
##     H2O API Extensions:         XGBoost, Algos, AutoML, Core V3, Core V4 
##     R Version:                  R version 3.5.2 (2018-12-20) 
## 
```


Make functions to determine cutoff for other missense predictors



Compute for BRCA1



![](BRCA_ML_files/figure-html/merge-brca1-1.png)<!-- -->

Compute for BRCA2
![](BRCA_ML_files/figure-html/compute-brca2-1.png)<!-- -->


Make a combined image
![](BRCA_ML_files/figure-html/combine-p1-p2-1.png)<!-- -->

Make PR Curves
![](BRCA_ML_files/figure-html/make-pr-curves-1.png)<!-- -->![](BRCA_ML_files/figure-html/make-pr-curves-2.png)<!-- -->![](BRCA_ML_files/figure-html/make-pr-curves-3.png)<!-- -->
All together now
![](BRCA_ML_files/figure-html/plot-pr-1.png)<!-- -->

Make pair plot
![](BRCA_ML_files/figure-html/pair-plot-1.png)<!-- -->

Print AUPROC

```
##      threshold sensitivity specificity   tp tn fp  fn        mcc                  model      PRauc  Gene
## 1   0.50463650  0.96093750  0.86486486  123 64 10   5 0.83897033                BRCA-ML 0.95635057 BRCA2
## 2   0.09124986  0.94817518  0.92857143 1299 78  6  71 0.67449018                BRCA-ML 0.88857558 BRCA1
## 3   0.29937650  0.85294118  0.83050847   87 49 10  15 0.67290966               BayesDel 0.82886006 BRCA2
## 4   0.28744524  0.76530612  0.84745763   75 50  9  23 0.59501928              MCapScore 0.00000000 BRCA2
## 5   0.66550000  0.66666667  0.93220339   68 55  4  34 0.58035996             RevelScore 0.77389379 BRCA2
## 6   0.80550000  0.80392157  0.77966102   82 46 13  20 0.57173482             Vest3Score 0.68749167 BRCA2
## 7   0.73920000  0.71568627  0.86440678   73 51  8  29 0.55903147            MetalrScore 0.76499182 BRCA2
## 8   0.76750000  0.81578947  0.72222222   62 39 15  14 0.53952565           MutpredScore 0.00000000 BRCA2
## 9   0.36070000  0.60784314  0.93220339   62 55  4  40 0.52908262           MetasvmScore 0.77965613 BRCA2
## 10 -2.43500000  0.89215686  0.57627119   91 34 25  11 0.50296033            FathmmScore 0.23772313 BRCA2
## 11  0.98750000  0.55882353  0.93220339   57 55  4  45 0.48771541     Polyphen2HvarScore 0.60575001 BRCA2
## 12  0.81084154  0.76470588  0.69491525   78 41 18  24 0.45137114               EigenRaw 0.70829596 BRCA2
## 13  0.73882692  0.69607843  0.72881356   71 43 16  31 0.41079878             EigenPcRaw 0.66636386 BRCA2
## 14  0.85400000  0.94463668  0.77777778  819 14  4  48 0.39951811           MutpredScore 0.00000000 BRCA1
## 15  0.00000050  0.64705882  0.76271186   66 45 14  36 0.39489304               LrtScore 0.27371018 BRCA2
## 16  0.99999999  0.75490196  0.61016949   77 36 23  25 0.36260939        GenocanyonScore 0.53286356 BRCA2
## 17  0.94256000  0.45098039  0.89830508   46 53  6  56 0.35990816   FathmmMklCodingScore 0.50315656 BRCA2
## 18  0.16000000  0.77828467  0.89552239  853 60  7 243 0.35769074         AlignGVGDPrior 0.22992242 BRCA1
## 19  0.47500000  0.56862745  0.79661017   58 47 12  44 0.35500288         AlignGVGDPrior 0.52009029 BRCA2
## 20  0.00050000  0.54901961  0.81355932   56 48 11  46 0.35442743              SiftScore 0.25855403 BRCA2
## 21 -1.69500000  0.63725490  0.72881356   65 43 16  37 0.35277752           ProveanScore 0.24861842 BRCA2
## 22  0.82250000  0.77737226  0.88059701  852 59  8 244 0.34928325             Vest3Score 0.28788479 BRCA1
## 23  0.99992050  0.48039216  0.86440678   49 51  8  53 0.34740648    MutationtasterScore 0.49785279 BRCA2
## 24  5.27148350  0.35294118  0.94915254   36 56  3  66 0.33974695                CaddRaw 0.51075749 BRCA2
## 25  0.71700000  0.80200730  0.80597015  879 54 13 217 0.33509124             RevelScore 0.23897527 BRCA1
## 26  0.33291750  0.76003650  0.86567164  833 58  9 263 0.32614139               BayesDel 0.21665934 BRCA1
## 27  0.00950000  0.78102190  0.80597015  856 54 13 240 0.31469609              SiftScore 0.03195874 BRCA1
## 28  0.50056819  0.76459854  0.82089552  838 55 12 258 0.30759802             EigenPcRaw 0.18944944 BRCA1
## 29  0.48704278  0.73540146  0.86567164  806 58  9 290 0.30584511               EigenRaw 0.19559553 BRCA1
## 30  0.99850000  0.76733577  0.80597015  841 54 13 255 0.30242725     Polyphen2HvarScore 0.21799701 BRCA1
## 31  0.62648930  0.73572744  0.85074627  799 57 10 287 0.29987992              MCapScore 0.00000000 BRCA1
## 32 -3.82500000  0.68339416  0.83582090  749 56 11 347 0.25423225           ProveanScore 0.03316140 BRCA1
## 33  0.00030500  0.72262774  0.77611940  792 52 15 304 0.25215137               LrtScore 0.03506381 BRCA1
## 34  0.71682950  0.64689781  0.85074627  709 57 10 387 0.23867426    Gm12878FitconsScore 0.10616048 BRCA1
## 35  0.99950000  0.58850365  0.91044776  645 61  6 451 0.23419384     Polyphen2HdivScore 0.11587771 BRCA1
## 36  0.66725550  0.34313725  0.86440678   35 51  8  67 0.22602741     H1HescFitconsScore 0.32518476 BRCA2
## 37  0.71947300  0.70529197  0.74626866  773 50 17 323 0.22541914 IntegratedFitconsScore 0.12207350 BRCA1
## 38  0.83088500  0.53102190  0.92537313  582 62  5 514 0.21269305   FathmmMklCodingScore 0.10619715 BRCA1
## 39  0.98464723  0.20588235  0.94915254   21 56  3  81 0.20974405              DannScore 0.43303285 BRCA2
## 40  0.64960000  0.69981752  0.71641791  767 48 19 329 0.20720478           MetasvmScore 0.14109121 BRCA1
## 41  0.01255058  0.53010949  0.89552239  581 60  7 515 0.19836006        GenocanyonScore 0.07165259 BRCA1
## 42  5.73541700  0.77098540  0.56716418  845 38 29 251 0.18232528                CaddRaw 0.21501036 BRCA1
## 43  0.77505000  0.61952555  0.71641791  679 48 19 417 0.15979187            MetalrScore 0.15233433 BRCA1
## 44  0.71723900  0.68795620  0.61194030  754 41 26 342 0.14868508     H1HescFitconsScore 0.13435074 BRCA1
## 45  0.67874900  0.25490196  0.86440678   26 51  8  76 0.14085006 IntegratedFitconsScore 0.33394710 BRCA2
## 46  0.99172894  0.40054745  0.80597015  439 54 13 657 0.09871786              DannScore 0.10790486 BRCA1
## 47 -4.01000000  0.32481752  0.85074627  356 57 10 740 0.08808621            FathmmScore 0.06508961 BRCA1
## 48  0.99987900  0.99726277  0.01492537 1093  1 66   3 0.04850747    MutationtasterScore 0.05808745 BRCA1
## 49  0.66945200  0.64705882  0.38983051   66 23 36  36 0.03688933      HuvecFitconsScore 0.34071703 BRCA2
## 50  0.65331700  0.69799270  0.37313433  765 25 42 331 0.03595966      HuvecFitconsScore 0.08734389 BRCA1
## 51  0.71682950  0.04901961  0.96610169    5 57  2  97 0.03572777    Gm12878FitconsScore 0.39562637 BRCA2
## 52        -Inf  1.00000000  0.00000000  102  0 59   0        NaN     Polyphen2HdivScore 0.50874078 BRCA2
## 53         Inf  0.00000000  1.00000000    0 59  0 102        NaN     Polyphen2HdivScore 0.50874078 BRCA2
```

![](BRCA_ML_files/figure-html/print-summary-1.png)<!-- -->

Save table S2.



