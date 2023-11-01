compare_clusters
================
Natalie Davidson
8/10/2021

# Intro

In this notebook I compare the clusters from the published output of the
Way pipeline and my run. Since the clusters are mapped to TCGA clusters,
we can use the cluster IDs to compare them. To do this, I treated the
cluster IDs as expected (published) and predicted (my run) labels. This
enables me to calculate balanced accuracy between the published results
and my pipeline results. We do this across all datasets and in Kmeans
and NMF.

This notebook contains a lot of detailed results, but the important
results are at the bottom in the section “Balanced Accuracy Summary”.
The summary contains the balanced accuracy across all datasets and sizes
of K (2-4). Balanced accuracy is used when the cluster sizes differ –
which in this case they do slightly.

It is calculated as (Sensitivity + Specificity) / 2

Sensitivity (true positive rate) = TP / (TP + FN)

Specificity (true negative rate) = TN / (TN + FT)

# First read the files and process

``` r
new_file = file.path(proj_dir, 
                    "/data/way_pipeline_results_10removed_NeoRemoved/2.Clustering_DiffExprs-Tables-ClusterMembership/FullClusterMembership.csv")

new_df = data.frame(fread(new_file))
colnames(new_df)[1] = "ID"

old_file = file.path(proj_dir, 
                    "/data/way_publication_results/033514_tables3.csv")

old_df = data.frame(fread(old_file))
old_df$Dataset[which(old_df$Dataset == "Mayo")] = "mayo.eset"

# only take samples that are in both
keep_ids = intersect(new_df$ID, old_df$ID)
new_df = subset(new_df, ID %in% keep_ids)
old_df = subset(old_df, ID %in% keep_ids)

# order them so we can directly compare
new_df = new_df[order(new_df$ID), ]
old_df = old_df[order(old_df$ID), ]
stopifnot(all(new_df$ID == old_df$ID))
```

# Plot confusion matrices

``` r
get_conf_matrix <- function(curr_dataset, curr_clust, new_df, old_df){
    

    expected_val = subset(old_df, Dataset == curr_dataset)
    predicted_val = subset(new_df, Dataset == curr_dataset)
    stopifnot(all(expected_val$ID == predicted_val$ID))
    
    expected_val = expected_val[,curr_clust]
    predicted_val = predicted_val[,curr_clust]
    
    conf_res <- confusionMatrix(data=as.factor(predicted_val), 
                                reference = as.factor(expected_val))
    return(conf_res)

}
```

### TCGA Kmeans

``` r
print("K=2")
```

    ## [1] "K=2"

``` r
conf_res2 <- get_conf_matrix("TCGA", "ClusterK2_kmeans", new_df, old_df)
print(conf_res2)
```

    ## Confusion Matrix and Statistics
    ## 
    ##           Reference
    ## Prediction   1   2
    ##          1 285  28
    ##          2   2 184
    ##                                           
    ##                Accuracy : 0.9399          
    ##                  95% CI : (0.9153, 0.9591)
    ##     No Information Rate : 0.5752          
    ##     P-Value [Acc > NIR] : < 2.2e-16       
    ##                                           
    ##                   Kappa : 0.875           
    ##                                           
    ##  Mcnemar's Test P-Value : 5.01e-06        
    ##                                           
    ##             Sensitivity : 0.9930          
    ##             Specificity : 0.8679          
    ##          Pos Pred Value : 0.9105          
    ##          Neg Pred Value : 0.9892          
    ##              Prevalence : 0.5752          
    ##          Detection Rate : 0.5711          
    ##    Detection Prevalence : 0.6273          
    ##       Balanced Accuracy : 0.9305          
    ##                                           
    ##        'Positive' Class : 1               
    ## 

``` r
print("K=3")
```

    ## [1] "K=3"

``` r
conf_res3 <- get_conf_matrix("TCGA", "ClusterK3_kmeans", new_df, old_df)
print(conf_res3)
```

    ## Confusion Matrix and Statistics
    ## 
    ##           Reference
    ## Prediction   1   2   3
    ##          1 135   2   3
    ##          2   0 132   1
    ##          3   2   5 219
    ## 
    ## Overall Statistics
    ##                                           
    ##                Accuracy : 0.9739          
    ##                  95% CI : (0.9559, 0.9861)
    ##     No Information Rate : 0.4469          
    ##     P-Value [Acc > NIR] : <2e-16          
    ##                                           
    ##                   Kappa : 0.9597          
    ##                                           
    ##  Mcnemar's Test P-Value : 0.1818          
    ## 
    ## Statistics by Class:
    ## 
    ##                      Class: 1 Class: 2 Class: 3
    ## Sensitivity            0.9854   0.9496   0.9821
    ## Specificity            0.9862   0.9972   0.9746
    ## Pos Pred Value         0.9643   0.9925   0.9690
    ## Neg Pred Value         0.9944   0.9809   0.9853
    ## Prevalence             0.2745   0.2786   0.4469
    ## Detection Rate         0.2705   0.2645   0.4389
    ## Detection Prevalence   0.2806   0.2665   0.4529
    ## Balanced Accuracy      0.9858   0.9734   0.9784

``` r
print("K=4")
```

    ## [1] "K=4"

``` r
conf_res4 <- get_conf_matrix("TCGA", "ClusterK4_kmeans", new_df, old_df)
print(conf_res4)
```

    ## Confusion Matrix and Statistics
    ## 
    ##           Reference
    ## Prediction   1   2   3   4
    ##          1 119   1   1   2
    ##          2   0  94   0   2
    ##          3   9   0  97  14
    ##          4   2   4   2 152
    ## 
    ## Overall Statistics
    ##                                           
    ##                Accuracy : 0.9259          
    ##                  95% CI : (0.8992, 0.9473)
    ##     No Information Rate : 0.3407          
    ##     P-Value [Acc > NIR] : < 2.2e-16       
    ##                                           
    ##                   Kappa : 0.8998          
    ##                                           
    ##  Mcnemar's Test P-Value : NA              
    ## 
    ## Statistics by Class:
    ## 
    ##                      Class: 1 Class: 2 Class: 3 Class: 4
    ## Sensitivity            0.9154   0.9495   0.9700   0.8941
    ## Specificity            0.9892   0.9950   0.9424   0.9757
    ## Pos Pred Value         0.9675   0.9792   0.8083   0.9500
    ## Neg Pred Value         0.9707   0.9876   0.9921   0.9469
    ## Prevalence             0.2605   0.1984   0.2004   0.3407
    ## Detection Rate         0.2385   0.1884   0.1944   0.3046
    ## Detection Prevalence   0.2465   0.1924   0.2405   0.3206
    ## Balanced Accuracy      0.9523   0.9722   0.9562   0.9349

``` r
BA_tcga_kmeans = c(conf_res2$byClass[11], conf_res3$byClass[,11], conf_res4$byClass[,11])
BA_tcga_kmeans
```

    ## Balanced Accuracy          Class: 1          Class: 2          Class: 3 
    ##         0.9304779         0.9857947         0.9734313         0.9783502 
    ##          Class: 1          Class: 2          Class: 3          Class: 4 
    ##         0.9522723         0.9722475         0.9561779         0.9349008

### TCGA NMF

``` r
print("K=2")
```

    ## [1] "K=2"

``` r
conf_res2 <- get_conf_matrix("TCGA", "ClusterK2_NMF", new_df, old_df)
print(conf_res2)
```

    ## Confusion Matrix and Statistics
    ## 
    ##           Reference
    ## Prediction   1   2
    ##          1 245   2
    ##          2  15 237
    ##                                        
    ##                Accuracy : 0.9659       
    ##                  95% CI : (0.946, 0.98)
    ##     No Information Rate : 0.521        
    ##     P-Value [Acc > NIR] : < 2.2e-16    
    ##                                        
    ##                   Kappa : 0.9319       
    ##                                        
    ##  Mcnemar's Test P-Value : 0.003609     
    ##                                        
    ##             Sensitivity : 0.9423       
    ##             Specificity : 0.9916       
    ##          Pos Pred Value : 0.9919       
    ##          Neg Pred Value : 0.9405       
    ##              Prevalence : 0.5210       
    ##          Detection Rate : 0.4910       
    ##    Detection Prevalence : 0.4950       
    ##       Balanced Accuracy : 0.9670       
    ##                                        
    ##        'Positive' Class : 1            
    ## 

``` r
print("K=3")
```

    ## [1] "K=3"

``` r
conf_res3 <- get_conf_matrix("TCGA", "ClusterK3_NMF", new_df, old_df)
print(conf_res3)
```

    ## Confusion Matrix and Statistics
    ## 
    ##           Reference
    ## Prediction   1   2   3
    ##          1 139   7   1
    ##          2   0 150   1
    ##          3   7  10 184
    ## 
    ## Overall Statistics
    ##                                           
    ##                Accuracy : 0.9479          
    ##                  95% CI : (0.9246, 0.9657)
    ##     No Information Rate : 0.3727          
    ##     P-Value [Acc > NIR] : < 2.2e-16       
    ##                                           
    ##                   Kappa : 0.9213          
    ##                                           
    ##  Mcnemar's Test P-Value : 0.0002917       
    ## 
    ## Statistics by Class:
    ## 
    ##                      Class: 1 Class: 2 Class: 3
    ## Sensitivity            0.9521   0.8982   0.9892
    ## Specificity            0.9773   0.9970   0.9457
    ## Pos Pred Value         0.9456   0.9934   0.9154
    ## Neg Pred Value         0.9801   0.9511   0.9933
    ## Prevalence             0.2926   0.3347   0.3727
    ## Detection Rate         0.2786   0.3006   0.3687
    ## Detection Prevalence   0.2946   0.3026   0.4028
    ## Balanced Accuracy      0.9647   0.9476   0.9675

``` r
print("K=4")
```

    ## [1] "K=4"

``` r
conf_res4 <- get_conf_matrix("TCGA", "ClusterK4_NMF", new_df, old_df)
print(conf_res4)
```

    ## Confusion Matrix and Statistics
    ## 
    ##           Reference
    ## Prediction  1  2  3  4
    ##          1 83  0 27  6
    ##          2  2 90  5  0
    ##          3  0  8 72 96
    ##          4 15 43  0 52
    ## 
    ## Overall Statistics
    ##                                           
    ##                Accuracy : 0.5952          
    ##                  95% CI : (0.5507, 0.6386)
    ##     No Information Rate : 0.3086          
    ##     P-Value [Acc > NIR] : < 2.2e-16       
    ##                                           
    ##                   Kappa : 0.4652          
    ##                                           
    ##  Mcnemar's Test P-Value : < 2.2e-16       
    ## 
    ## Statistics by Class:
    ## 
    ##                      Class: 1 Class: 2 Class: 3 Class: 4
    ## Sensitivity            0.8300   0.6383   0.6923   0.3377
    ## Specificity            0.9173   0.9804   0.7367   0.8319
    ## Pos Pred Value         0.7155   0.9278   0.4091   0.4727
    ## Neg Pred Value         0.9556   0.8731   0.9009   0.7378
    ## Prevalence             0.2004   0.2826   0.2084   0.3086
    ## Detection Rate         0.1663   0.1804   0.1443   0.1042
    ## Detection Prevalence   0.2325   0.1944   0.3527   0.2204
    ## Balanced Accuracy      0.8736   0.8094   0.7145   0.5848

``` r
BA_tcga_nmf = c(conf_res2$byClass[11], conf_res3$byClass[,11], conf_res4$byClass[,11])
BA_tcga_nmf
```

    ## Balanced Accuracy          Class: 1          Class: 2          Class: 3 
    ##         0.9669697         0.9646960         0.9475958         0.9674671 
    ##          Class: 1          Class: 2          Class: 3          Class: 4 
    ##         0.8736466         0.8093724         0.7145083         0.5847732

### Mayo Kmeans

``` r
print("K=2")
```

    ## [1] "K=2"

``` r
conf_res2 <- get_conf_matrix("mayo.eset", "ClusterK2_kmeans", new_df, old_df)
print(conf_res2)
```

    ## Confusion Matrix and Statistics
    ## 
    ##           Reference
    ## Prediction   1   2
    ##          1 189   1
    ##          2  10 177
    ##                                           
    ##                Accuracy : 0.9708          
    ##                  95% CI : (0.9484, 0.9853)
    ##     No Information Rate : 0.5279          
    ##     P-Value [Acc > NIR] : < 2e-16         
    ##                                           
    ##                   Kappa : 0.9416          
    ##                                           
    ##  Mcnemar's Test P-Value : 0.01586         
    ##                                           
    ##             Sensitivity : 0.9497          
    ##             Specificity : 0.9944          
    ##          Pos Pred Value : 0.9947          
    ##          Neg Pred Value : 0.9465          
    ##              Prevalence : 0.5279          
    ##          Detection Rate : 0.5013          
    ##    Detection Prevalence : 0.5040          
    ##       Balanced Accuracy : 0.9721          
    ##                                           
    ##        'Positive' Class : 1               
    ## 

``` r
print("K=3")
```

    ## [1] "K=3"

``` r
conf_res3 <- get_conf_matrix("mayo.eset", "ClusterK3_kmeans", new_df, old_df)
print(conf_res3)
```

    ## Confusion Matrix and Statistics
    ## 
    ##           Reference
    ## Prediction   1   2   3
    ##          1 138   3   6
    ##          2   1  87   0
    ##          3   0  16 126
    ## 
    ## Overall Statistics
    ##                                           
    ##                Accuracy : 0.931           
    ##                  95% CI : (0.9006, 0.9545)
    ##     No Information Rate : 0.3687          
    ##     P-Value [Acc > NIR] : < 2.2e-16       
    ##                                           
    ##                   Kappa : 0.8953          
    ##                                           
    ##  Mcnemar's Test P-Value : 4.038e-05       
    ## 
    ## Statistics by Class:
    ## 
    ##                      Class: 1 Class: 2 Class: 3
    ## Sensitivity            0.9928   0.8208   0.9545
    ## Specificity            0.9622   0.9963   0.9347
    ## Pos Pred Value         0.9388   0.9886   0.8873
    ## Neg Pred Value         0.9957   0.9343   0.9745
    ## Prevalence             0.3687   0.2812   0.3501
    ## Detection Rate         0.3660   0.2308   0.3342
    ## Detection Prevalence   0.3899   0.2334   0.3767
    ## Balanced Accuracy      0.9775   0.9085   0.9446

``` r
print("K=4")
```

    ## [1] "K=4"

``` r
conf_res4 <- get_conf_matrix("mayo.eset", "ClusterK4_kmeans", new_df, old_df)
print(conf_res4)
```

    ## Confusion Matrix and Statistics
    ## 
    ##           Reference
    ## Prediction   1   2   3   4
    ##          1 103   1   1   1
    ##          2   2  75   1   1
    ##          3   4   0  90   1
    ##          4   0   3   1  93
    ## 
    ## Overall Statistics
    ##                                          
    ##                Accuracy : 0.9576         
    ##                  95% CI : (0.932, 0.9756)
    ##     No Information Rate : 0.2891         
    ##     P-Value [Acc > NIR] : <2e-16         
    ##                                          
    ##                   Kappa : 0.9432         
    ##                                          
    ##  Mcnemar's Test P-Value : 0.5268         
    ## 
    ## Statistics by Class:
    ## 
    ##                      Class: 1 Class: 2 Class: 3 Class: 4
    ## Sensitivity            0.9450   0.9494   0.9677   0.9688
    ## Specificity            0.9888   0.9866   0.9824   0.9858
    ## Pos Pred Value         0.9717   0.9494   0.9474   0.9588
    ## Neg Pred Value         0.9779   0.9866   0.9894   0.9893
    ## Prevalence             0.2891   0.2095   0.2467   0.2546
    ## Detection Rate         0.2732   0.1989   0.2387   0.2467
    ## Detection Prevalence   0.2812   0.2095   0.2520   0.2573
    ## Balanced Accuracy      0.9669   0.9680   0.9751   0.9773

``` r
BA_mayo_kmeans = c(conf_res2$byClass[11], conf_res3$byClass[,11], conf_res4$byClass[,11])
BA_mayo_kmeans
```

    ## Balanced Accuracy          Class: 1          Class: 2          Class: 3 
    ##         0.9720654         0.9774953         0.9085323         0.9446197 
    ##          Class: 1          Class: 2          Class: 3          Class: 4 
    ##         0.9668800         0.9679721         0.9750682         0.9772576

### Mayo NMF

``` r
print("K=2")
```

    ## [1] "K=2"

``` r
conf_res2 <- get_conf_matrix("mayo.eset", "ClusterK2_NMF", new_df, old_df)
print(conf_res2)
```

    ## Confusion Matrix and Statistics
    ## 
    ##           Reference
    ## Prediction   1   2
    ##          1 147   0
    ##          2  32 198
    ##                                           
    ##                Accuracy : 0.9151          
    ##                  95% CI : (0.8823, 0.9412)
    ##     No Information Rate : 0.5252          
    ##     P-Value [Acc > NIR] : < 2.2e-16       
    ##                                           
    ##                   Kappa : 0.8283          
    ##                                           
    ##  Mcnemar's Test P-Value : 4.251e-08       
    ##                                           
    ##             Sensitivity : 0.8212          
    ##             Specificity : 1.0000          
    ##          Pos Pred Value : 1.0000          
    ##          Neg Pred Value : 0.8609          
    ##              Prevalence : 0.4748          
    ##          Detection Rate : 0.3899          
    ##    Detection Prevalence : 0.3899          
    ##       Balanced Accuracy : 0.9106          
    ##                                           
    ##        'Positive' Class : 1               
    ## 

``` r
print("K=3")
```

    ## [1] "K=3"

``` r
conf_res3 <- get_conf_matrix("mayo.eset", "ClusterK3_NMF", new_df, old_df)
print(conf_res3)
```

    ## Confusion Matrix and Statistics
    ## 
    ##           Reference
    ## Prediction  1  2  3
    ##          1 93  2 37
    ##          2  5 97  1
    ##          3  0 45 97
    ## 
    ## Overall Statistics
    ##                                          
    ##                Accuracy : 0.7613         
    ##                  95% CI : (0.715, 0.8034)
    ##     No Information Rate : 0.382          
    ##     P-Value [Acc > NIR] : < 2.2e-16      
    ##                                          
    ##                   Kappa : 0.6436         
    ##                                          
    ##  Mcnemar's Test P-Value : < 2.2e-16      
    ## 
    ## Statistics by Class:
    ## 
    ##                      Class: 1 Class: 2 Class: 3
    ## Sensitivity            0.9490   0.6736   0.7185
    ## Specificity            0.8602   0.9742   0.8140
    ## Pos Pred Value         0.7045   0.9417   0.6831
    ## Neg Pred Value         0.9796   0.8285   0.8383
    ## Prevalence             0.2599   0.3820   0.3581
    ## Detection Rate         0.2467   0.2573   0.2573
    ## Detection Prevalence   0.3501   0.2732   0.3767
    ## Balanced Accuracy      0.9046   0.8239   0.7663

``` r
print("K=4")
```

    ## [1] "K=4"

``` r
conf_res4 <- get_conf_matrix("mayo.eset", "ClusterK4_NMF", new_df, old_df)
print(conf_res4)
```

    ## Confusion Matrix and Statistics
    ## 
    ##           Reference
    ## Prediction  1  2  3  4
    ##          1 64  0 37  0
    ##          2  1 78  9  0
    ##          3  0  1 40 79
    ##          4 26 17  0 25
    ## 
    ## Overall Statistics
    ##                                           
    ##                Accuracy : 0.5491          
    ##                  95% CI : (0.4973, 0.6001)
    ##     No Information Rate : 0.2759          
    ##     P-Value [Acc > NIR] : < 2.2e-16       
    ##                                           
    ##                   Kappa : 0.4016          
    ##                                           
    ##  Mcnemar's Test P-Value : < 2.2e-16       
    ## 
    ## Statistics by Class:
    ## 
    ##                      Class: 1 Class: 2 Class: 3 Class: 4
    ## Sensitivity            0.7033   0.8125   0.4651  0.24038
    ## Specificity            0.8706   0.9644   0.7251  0.84249
    ## Pos Pred Value         0.6337   0.8864   0.3333  0.36765
    ## Neg Pred Value         0.9022   0.9377   0.8210  0.74434
    ## Prevalence             0.2414   0.2546   0.2281  0.27586
    ## Detection Rate         0.1698   0.2069   0.1061  0.06631
    ## Detection Prevalence   0.2679   0.2334   0.3183  0.18037
    ## Balanced Accuracy      0.7870   0.8885   0.5951  0.54144

``` r
BA_mayo_nmf = c(conf_res2$byClass[11], conf_res3$byClass[,11], conf_res4$byClass[,11])
BA_mayo_nmf
```

    ## Balanced Accuracy          Class: 1          Class: 2          Class: 3 
    ##         0.9106145         0.9045973         0.8239300         0.7662841 
    ##          Class: 1          Class: 2          Class: 3          Class: 4 
    ##         0.7869630         0.8884564         0.5951011         0.5414377

### Tothill Kmeans

``` r
print("K=2")
```

    ## [1] "K=2"

``` r
conf_res2 <- get_conf_matrix("Tothill", "ClusterK2_kmeans", new_df, old_df)
print(conf_res2)
```

    ## Confusion Matrix and Statistics
    ## 
    ##           Reference
    ## Prediction   1   2
    ##          1 117   1
    ##          2  12 111
    ##                                          
    ##                Accuracy : 0.9461         
    ##                  95% CI : (0.9095, 0.971)
    ##     No Information Rate : 0.5353         
    ##     P-Value [Acc > NIR] : < 2.2e-16      
    ##                                          
    ##                   Kappa : 0.8923         
    ##                                          
    ##  Mcnemar's Test P-Value : 0.005546       
    ##                                          
    ##             Sensitivity : 0.9070         
    ##             Specificity : 0.9911         
    ##          Pos Pred Value : 0.9915         
    ##          Neg Pred Value : 0.9024         
    ##              Prevalence : 0.5353         
    ##          Detection Rate : 0.4855         
    ##    Detection Prevalence : 0.4896         
    ##       Balanced Accuracy : 0.9490         
    ##                                          
    ##        'Positive' Class : 1              
    ## 

``` r
print("K=3")
```

    ## [1] "K=3"

``` r
conf_res3 <- get_conf_matrix("Tothill", "ClusterK3_kmeans", new_df, old_df)
print(conf_res3)
```

    ## Confusion Matrix and Statistics
    ## 
    ##           Reference
    ## Prediction  1  2  3
    ##          1 93  0  3
    ##          2  0 46  2
    ##          3 12  0 85
    ## 
    ## Overall Statistics
    ##                                           
    ##                Accuracy : 0.9295          
    ##                  95% CI : (0.8895, 0.9584)
    ##     No Information Rate : 0.4357          
    ##     P-Value [Acc > NIR] : < 2.2e-16       
    ##                                           
    ##                   Kappa : 0.8895          
    ##                                           
    ##  Mcnemar's Test P-Value : NA              
    ## 
    ## Statistics by Class:
    ## 
    ##                      Class: 1 Class: 2 Class: 3
    ## Sensitivity            0.8857   1.0000   0.9444
    ## Specificity            0.9779   0.9897   0.9205
    ## Pos Pred Value         0.9688   0.9583   0.8763
    ## Neg Pred Value         0.9172   1.0000   0.9653
    ## Prevalence             0.4357   0.1909   0.3734
    ## Detection Rate         0.3859   0.1909   0.3527
    ## Detection Prevalence   0.3983   0.1992   0.4025
    ## Balanced Accuracy      0.9318   0.9949   0.9325

``` r
print("K=4")
```

    ## [1] "K=4"

``` r
conf_res4 <- get_conf_matrix("Tothill", "ClusterK4_kmeans", new_df, old_df)
print(conf_res4)
```

    ## Confusion Matrix and Statistics
    ## 
    ##           Reference
    ## Prediction  1  2  3  4
    ##          1 73  0  0  2
    ##          2  0 43  2  0
    ##          3  0  0 64  3
    ##          4  0  0  1 53
    ## 
    ## Overall Statistics
    ##                                           
    ##                Accuracy : 0.9668          
    ##                  95% CI : (0.9356, 0.9856)
    ##     No Information Rate : 0.3029          
    ##     P-Value [Acc > NIR] : < 2.2e-16       
    ##                                           
    ##                   Kappa : 0.9552          
    ##                                           
    ##  Mcnemar's Test P-Value : NA              
    ## 
    ## Statistics by Class:
    ## 
    ##                      Class: 1 Class: 2 Class: 3 Class: 4
    ## Sensitivity            1.0000   1.0000   0.9552   0.9138
    ## Specificity            0.9881   0.9899   0.9828   0.9945
    ## Pos Pred Value         0.9733   0.9556   0.9552   0.9815
    ## Neg Pred Value         1.0000   1.0000   0.9828   0.9733
    ## Prevalence             0.3029   0.1784   0.2780   0.2407
    ## Detection Rate         0.3029   0.1784   0.2656   0.2199
    ## Detection Prevalence   0.3112   0.1867   0.2780   0.2241
    ## Balanced Accuracy      0.9940   0.9949   0.9690   0.9542

``` r
BA_tothill_kmeans = c(conf_res2$byClass[11], conf_res3$byClass[,11], conf_res4$byClass[,11])
BA_tothill_kmeans
```

    ## Balanced Accuracy          Class: 1          Class: 2          Class: 3 
    ##         0.9490241         0.9318277         0.9948718         0.9324871 
    ##          Class: 1          Class: 2          Class: 3          Class: 4 
    ##         0.9940476         0.9949495         0.9689913         0.9541643

### Tothill NMF

``` r
print("K=2")
```

    ## [1] "K=2"

``` r
conf_res2 <- get_conf_matrix("Tothill", "ClusterK2_NMF", new_df, old_df)
print(conf_res2)
```

    ## Confusion Matrix and Statistics
    ## 
    ##           Reference
    ## Prediction   1   2
    ##          1 115  12
    ##          2   0 114
    ##                                          
    ##                Accuracy : 0.9502         
    ##                  95% CI : (0.9146, 0.974)
    ##     No Information Rate : 0.5228         
    ##     P-Value [Acc > NIR] : < 2.2e-16      
    ##                                          
    ##                   Kappa : 0.9007         
    ##                                          
    ##  Mcnemar's Test P-Value : 0.001496       
    ##                                          
    ##             Sensitivity : 1.0000         
    ##             Specificity : 0.9048         
    ##          Pos Pred Value : 0.9055         
    ##          Neg Pred Value : 1.0000         
    ##              Prevalence : 0.4772         
    ##          Detection Rate : 0.4772         
    ##    Detection Prevalence : 0.5270         
    ##       Balanced Accuracy : 0.9524         
    ##                                          
    ##        'Positive' Class : 1              
    ## 

``` r
print("K=3")
```

    ## [1] "K=3"

``` r
conf_res3 <- get_conf_matrix("Tothill", "ClusterK3_NMF", new_df, old_df)
print(conf_res3)
```

    ## Confusion Matrix and Statistics
    ## 
    ##           Reference
    ## Prediction  1  2  3
    ##          1 63  1 18
    ##          2  6 60  1
    ##          3  0 29 63
    ## 
    ## Overall Statistics
    ##                                           
    ##                Accuracy : 0.7718          
    ##                  95% CI : (0.7135, 0.8232)
    ##     No Information Rate : 0.3734          
    ##     P-Value [Acc > NIR] : < 2.2e-16       
    ##                                           
    ##                   Kappa : 0.6588          
    ##                                           
    ##  Mcnemar's Test P-Value : 2.461e-10       
    ## 
    ## Statistics by Class:
    ## 
    ##                      Class: 1 Class: 2 Class: 3
    ## Sensitivity            0.9130   0.6667   0.7683
    ## Specificity            0.8895   0.9536   0.8176
    ## Pos Pred Value         0.7683   0.8955   0.6848
    ## Neg Pred Value         0.9623   0.8276   0.8725
    ## Prevalence             0.2863   0.3734   0.3402
    ## Detection Rate         0.2614   0.2490   0.2614
    ## Detection Prevalence   0.3402   0.2780   0.3817
    ## Balanced Accuracy      0.9013   0.8102   0.7930

``` r
print("K=4")
```

    ## [1] "K=4"

``` r
conf_res4 <- get_conf_matrix("Tothill", "ClusterK4_NMF", new_df, old_df)
print(conf_res4)
```

    ## Confusion Matrix and Statistics
    ## 
    ##           Reference
    ## Prediction  1  2  3  4
    ##          1 62  0  1  5
    ##          2  3 45  1  0
    ##          3  0  0 50  1
    ##          4  0 12 15 46
    ## 
    ## Overall Statistics
    ##                                           
    ##                Accuracy : 0.8423          
    ##                  95% CI : (0.7901, 0.8859)
    ##     No Information Rate : 0.278           
    ##     P-Value [Acc > NIR] : < 2.2e-16       
    ##                                           
    ##                   Kappa : 0.7902          
    ##                                           
    ##  Mcnemar's Test P-Value : 6.019e-06       
    ## 
    ## Statistics by Class:
    ## 
    ##                      Class: 1 Class: 2 Class: 3 Class: 4
    ## Sensitivity            0.9538   0.7895   0.7463   0.8846
    ## Specificity            0.9659   0.9783   0.9943   0.8571
    ## Pos Pred Value         0.9118   0.9184   0.9804   0.6301
    ## Neg Pred Value         0.9827   0.9375   0.9105   0.9643
    ## Prevalence             0.2697   0.2365   0.2780   0.2158
    ## Detection Rate         0.2573   0.1867   0.2075   0.1909
    ## Detection Prevalence   0.2822   0.2033   0.2116   0.3029
    ## Balanced Accuracy      0.9599   0.8839   0.8703   0.8709

``` r
BA_tothill_nmf = c(conf_res2$byClass[11], conf_res3$byClass[,11], conf_res4$byClass[,11])
BA_tothill_nmf
```

    ## Balanced Accuracy          Class: 1          Class: 2          Class: 3 
    ##         0.9523810         0.9012892         0.8101545         0.7929514 
    ##          Class: 1          Class: 2          Class: 3          Class: 4 
    ##         0.9598776         0.8838673         0.8702608         0.8708791

### Yoshihara Kmeans

``` r
print("K=2")
```

    ## [1] "K=2"

``` r
conf_res2 <- get_conf_matrix("Yoshihara", "ClusterK2_kmeans", new_df, old_df)
print(conf_res2)
```

    ## Confusion Matrix and Statistics
    ## 
    ##           Reference
    ## Prediction   1   2
    ##          1  95   3
    ##          2   1 156
    ##                                           
    ##                Accuracy : 0.9843          
    ##                  95% CI : (0.9603, 0.9957)
    ##     No Information Rate : 0.6235          
    ##     P-Value [Acc > NIR] : <2e-16          
    ##                                           
    ##                   Kappa : 0.9667          
    ##                                           
    ##  Mcnemar's Test P-Value : 0.6171          
    ##                                           
    ##             Sensitivity : 0.9896          
    ##             Specificity : 0.9811          
    ##          Pos Pred Value : 0.9694          
    ##          Neg Pred Value : 0.9936          
    ##              Prevalence : 0.3765          
    ##          Detection Rate : 0.3725          
    ##    Detection Prevalence : 0.3843          
    ##       Balanced Accuracy : 0.9854          
    ##                                           
    ##        'Positive' Class : 1               
    ## 

``` r
print("K=3")
```

    ## [1] "K=3"

``` r
conf_res3 <- get_conf_matrix("Yoshihara", "ClusterK3_kmeans", new_df, old_df)
print(conf_res3)
```

    ## Confusion Matrix and Statistics
    ## 
    ##           Reference
    ## Prediction   1   2   3
    ##          1  90   1   3
    ##          2   1  37   0
    ##          3   0   6 117
    ## 
    ## Overall Statistics
    ##                                           
    ##                Accuracy : 0.9569          
    ##                  95% CI : (0.9241, 0.9783)
    ##     No Information Rate : 0.4706          
    ##     P-Value [Acc > NIR] : < 2e-16         
    ##                                           
    ##                   Kappa : 0.9299          
    ##                                           
    ##  Mcnemar's Test P-Value : 0.02929         
    ## 
    ## Statistics by Class:
    ## 
    ##                      Class: 1 Class: 2 Class: 3
    ## Sensitivity            0.9890   0.8409   0.9750
    ## Specificity            0.9756   0.9953   0.9556
    ## Pos Pred Value         0.9574   0.9737   0.9512
    ## Neg Pred Value         0.9938   0.9677   0.9773
    ## Prevalence             0.3569   0.1725   0.4706
    ## Detection Rate         0.3529   0.1451   0.4588
    ## Detection Prevalence   0.3686   0.1490   0.4824
    ## Balanced Accuracy      0.9823   0.9181   0.9653

``` r
print("K=4")
```

    ## [1] "K=4"

``` r
conf_res4 <- get_conf_matrix("Yoshihara", "ClusterK4_kmeans", new_df, old_df)
print(conf_res4)
```

    ## Confusion Matrix and Statistics
    ## 
    ##           Reference
    ## Prediction  1  2  3  4
    ##          1 87  0  1  0
    ##          2  0 29  0  0
    ##          3  0  1 71  1
    ##          4  0  0  2 63
    ## 
    ## Overall Statistics
    ##                                           
    ##                Accuracy : 0.9804          
    ##                  95% CI : (0.9548, 0.9936)
    ##     No Information Rate : 0.3412          
    ##     P-Value [Acc > NIR] : < 2.2e-16       
    ##                                           
    ##                   Kappa : 0.9728          
    ##                                           
    ##  Mcnemar's Test P-Value : NA              
    ## 
    ## Statistics by Class:
    ## 
    ##                      Class: 1 Class: 2 Class: 3 Class: 4
    ## Sensitivity            1.0000   0.9667   0.9595   0.9844
    ## Specificity            0.9940   1.0000   0.9890   0.9895
    ## Pos Pred Value         0.9886   1.0000   0.9726   0.9692
    ## Neg Pred Value         1.0000   0.9956   0.9835   0.9947
    ## Prevalence             0.3412   0.1176   0.2902   0.2510
    ## Detection Rate         0.3412   0.1137   0.2784   0.2471
    ## Detection Prevalence   0.3451   0.1137   0.2863   0.2549
    ## Balanced Accuracy      0.9970   0.9833   0.9742   0.9870

``` r
BA_yoshihara_kmeans = c(conf_res2$byClass[11], conf_res3$byClass[,11], conf_res4$byClass[,11])
BA_yoshihara_kmeans
```

    ## Balanced Accuracy          Class: 1          Class: 2          Class: 3 
    ##         0.9853577         0.9823104         0.9180849         0.9652778 
    ##          Class: 1          Class: 2          Class: 3          Class: 4 
    ##         0.9970238         0.9833333         0.9742049         0.9869519

### Yoshihara NMF

``` r
print("K=2")
```

    ## [1] "K=2"

``` r
conf_res2 <- get_conf_matrix("Yoshihara", "ClusterK2_NMF", new_df, old_df)
print(conf_res2)
```

    ## Confusion Matrix and Statistics
    ## 
    ##           Reference
    ## Prediction   1   2
    ##          1  97   5
    ##          2   2 151
    ##                                           
    ##                Accuracy : 0.9725          
    ##                  95% CI : (0.9443, 0.9889)
    ##     No Information Rate : 0.6118          
    ##     P-Value [Acc > NIR] : <2e-16          
    ##                                           
    ##                   Kappa : 0.9425          
    ##                                           
    ##  Mcnemar's Test P-Value : 0.4497          
    ##                                           
    ##             Sensitivity : 0.9798          
    ##             Specificity : 0.9679          
    ##          Pos Pred Value : 0.9510          
    ##          Neg Pred Value : 0.9869          
    ##              Prevalence : 0.3882          
    ##          Detection Rate : 0.3804          
    ##    Detection Prevalence : 0.4000          
    ##       Balanced Accuracy : 0.9739          
    ##                                           
    ##        'Positive' Class : 1               
    ## 

``` r
print("K=3")
```

    ## [1] "K=3"

``` r
conf_res3 <- get_conf_matrix("Yoshihara", "ClusterK3_NMF", new_df, old_df)
print(conf_res3)
```

    ## Confusion Matrix and Statistics
    ## 
    ##           Reference
    ## Prediction   1   2   3
    ##          1  63   0   1
    ##          2   0  68   0
    ##          3   4  13 106
    ## 
    ## Overall Statistics
    ##                                           
    ##                Accuracy : 0.9294          
    ##                  95% CI : (0.8907, 0.9576)
    ##     No Information Rate : 0.4196          
    ##     P-Value [Acc > NIR] : < 2.2e-16       
    ##                                           
    ##                   Kappa : 0.8909          
    ##                                           
    ##  Mcnemar's Test P-Value : NA              
    ## 
    ## Statistics by Class:
    ## 
    ##                      Class: 1 Class: 2 Class: 3
    ## Sensitivity            0.9403   0.8395   0.9907
    ## Specificity            0.9947   1.0000   0.8851
    ## Pos Pred Value         0.9844   1.0000   0.8618
    ## Neg Pred Value         0.9791   0.9305   0.9924
    ## Prevalence             0.2627   0.3176   0.4196
    ## Detection Rate         0.2471   0.2667   0.4157
    ## Detection Prevalence   0.2510   0.2667   0.4824
    ## Balanced Accuracy      0.9675   0.9198   0.9379

``` r
print("K=4")
```

    ## [1] "K=4"

``` r
conf_res4 <- get_conf_matrix("Yoshihara", "ClusterK4_NMF", new_df, old_df)
print(conf_res4)
```

    ## Confusion Matrix and Statistics
    ## 
    ##           Reference
    ## Prediction  1  2  3  4
    ##          1 47  7 11  0
    ##          2  0 36  8  0
    ##          3  0  0 90 15
    ##          4  1  7  1 32
    ## 
    ## Overall Statistics
    ##                                           
    ##                Accuracy : 0.8039          
    ##                  95% CI : (0.7498, 0.8508)
    ##     No Information Rate : 0.4314          
    ##     P-Value [Acc > NIR] : < 2.2e-16       
    ##                                           
    ##                   Kappa : 0.7242          
    ##                                           
    ##  Mcnemar's Test P-Value : 2.64e-08        
    ## 
    ## Statistics by Class:
    ## 
    ##                      Class: 1 Class: 2 Class: 3 Class: 4
    ## Sensitivity            0.9792   0.7200   0.8182   0.6809
    ## Specificity            0.9130   0.9610   0.8966   0.9567
    ## Pos Pred Value         0.7231   0.8182   0.8571   0.7805
    ## Neg Pred Value         0.9947   0.9336   0.8667   0.9299
    ## Prevalence             0.1882   0.1961   0.4314   0.1843
    ## Detection Rate         0.1843   0.1412   0.3529   0.1255
    ## Detection Prevalence   0.2549   0.1725   0.4118   0.1608
    ## Balanced Accuracy      0.9461   0.8405   0.8574   0.8188

``` r
BA_yoshihara_nmf = c(conf_res2$byClass[11], conf_res3$byClass[,11], conf_res4$byClass[,11])
BA_yoshihara_nmf
```

    ## Balanced Accuracy          Class: 1          Class: 2          Class: 3 
    ##         0.9738733         0.9674897         0.9197531         0.9378947 
    ##          Class: 1          Class: 2          Class: 3          Class: 4 
    ##         0.9461051         0.8404878         0.8573668         0.8187909

# Balanced Accuracy Summary

``` r
class_ids = c("K=2", "K=3, 1", "K=3, 2", "K=3, 3", "K=4, 1", "K=4, 2", "K=4, 3", "K=4, 4")


# Balanced Accuracy for Kmeans
ba_df = data.frame(class_ids,
                    BA_tcga_kmeans,
                    BA_mayo_kmeans,
                    BA_tothill_kmeans,
                    BA_yoshihara_kmeans)
ba_df
```

    ##   class_ids BA_tcga_kmeans BA_mayo_kmeans BA_tothill_kmeans BA_yoshihara_kmeans
    ## 1       K=2      0.9304779      0.9720654         0.9490241           0.9853577
    ## 2    K=3, 1      0.9857947      0.9774953         0.9318277           0.9823104
    ## 3    K=3, 2      0.9734313      0.9085323         0.9948718           0.9180849
    ## 4    K=3, 3      0.9783502      0.9446197         0.9324871           0.9652778
    ## 5    K=4, 1      0.9522723      0.9668800         0.9940476           0.9970238
    ## 6    K=4, 2      0.9722475      0.9679721         0.9949495           0.9833333
    ## 7    K=4, 3      0.9561779      0.9750682         0.9689913           0.9742049
    ## 8    K=4, 4      0.9349008      0.9772576         0.9541643           0.9869519

``` r
# Balanced Accuracy for NMF
ba_df = data.frame(class_ids,
                    BA_tcga_nmf,
                    BA_mayo_nmf,
                    BA_tothill_nmf,
                    BA_yoshihara_nmf)
ba_df
```

    ##   class_ids BA_tcga_nmf BA_mayo_nmf BA_tothill_nmf BA_yoshihara_nmf
    ## 1       K=2   0.9669697   0.9106145      0.9523810        0.9738733
    ## 2    K=3, 1   0.9646960   0.9045973      0.9012892        0.9674897
    ## 3    K=3, 2   0.9475958   0.8239300      0.8101545        0.9197531
    ## 4    K=3, 3   0.9674671   0.7662841      0.7929514        0.9378947
    ## 5    K=4, 1   0.8736466   0.7869630      0.9598776        0.9461051
    ## 6    K=4, 2   0.8093724   0.8884564      0.8838673        0.8404878
    ## 7    K=4, 3   0.7145083   0.5951011      0.8702608        0.8573668
    ## 8    K=4, 4   0.5847732   0.5414377      0.8708791        0.8187909

# Conclusion

We find that overall the clusters align between the previously published
Way et. al. run and my run of the pipeline. NOTE: There exist slight
data and possibly code differences between our two runs, so we do not
expect 100% matching between runs. However, we do expect high coherence
between runs. We find that the clusters are more consistent in K-Means
than NMF. We also find that K=4 is the most inconsistent between the two
runs of the pipeline. We will consider K=4 and KMeans results for
remaining main analyses.
