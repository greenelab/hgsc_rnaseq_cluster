compare_clusters
================
Natalie Davidson
9/08/2021

# Intro

In this notebook I compare the clusters from the two runs of the
pipeline. The runs differ in the fact that they either have or don’t
have samples with neoadjuvant therapy. Since the clusters are mapped to
TCGA clusters, we can use the cluster IDs to compare them. To do this, I
treated the cluster IDs as expected (with neoadj) and predicted (withOUT
neoadj) labels. This enables me to calculate balanced accuracy between
the published results and my pipeline results. We do this across all
datasets and in Kmeans and NMF.

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
                    "/data/way_pipeline_results_10removed/2.Clustering_DiffExprs-Tables-ClusterMembership/FullClusterMembership.csv")

old_df = data.frame(fread(old_file))
colnames(old_df)[1] = "ID"

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

### aaces.eset Kmeans

``` r
print("K=2")
```

    ## [1] "K=2"

``` r
conf_res2 <- get_conf_matrix("aaces.eset", "ClusterK2_kmeans", new_df, old_df)
print(conf_res2)
```

    ## Confusion Matrix and Statistics
    ## 
    ##           Reference
    ## Prediction   1   2
    ##          1 107  31
    ##          2   0 124
    ##                                           
    ##                Accuracy : 0.8817          
    ##                  95% CI : (0.8363, 0.9182)
    ##     No Information Rate : 0.5916          
    ##     P-Value [Acc > NIR] : < 2.2e-16       
    ##                                           
    ##                   Kappa : 0.7657          
    ##                                           
    ##  Mcnemar's Test P-Value : 7.118e-08       
    ##                                           
    ##             Sensitivity : 1.0000          
    ##             Specificity : 0.8000          
    ##          Pos Pred Value : 0.7754          
    ##          Neg Pred Value : 1.0000          
    ##              Prevalence : 0.4084          
    ##          Detection Rate : 0.4084          
    ##    Detection Prevalence : 0.5267          
    ##       Balanced Accuracy : 0.9000          
    ##                                           
    ##        'Positive' Class : 1               
    ## 

``` r
print("K=3")
```

    ## [1] "K=3"

``` r
conf_res3 <- get_conf_matrix("aaces.eset", "ClusterK3_kmeans", new_df, old_df)
print(conf_res3)
```

    ## Confusion Matrix and Statistics
    ## 
    ##           Reference
    ## Prediction  1  2  3
    ##          1 70  6 23
    ##          2  0 71  0
    ##          3  1  1 90
    ## 
    ## Overall Statistics
    ##                                           
    ##                Accuracy : 0.8817          
    ##                  95% CI : (0.8363, 0.9182)
    ##     No Information Rate : 0.4313          
    ##     P-Value [Acc > NIR] : < 2.2e-16       
    ##                                           
    ##                   Kappa : 0.8222          
    ##                                           
    ##  Mcnemar's Test P-Value : 5.432e-06       
    ## 
    ## Statistics by Class:
    ## 
    ##                      Class: 1 Class: 2 Class: 3
    ## Sensitivity            0.9859   0.9103   0.7965
    ## Specificity            0.8482   1.0000   0.9866
    ## Pos Pred Value         0.7071   1.0000   0.9783
    ## Neg Pred Value         0.9939   0.9634   0.8647
    ## Prevalence             0.2710   0.2977   0.4313
    ## Detection Rate         0.2672   0.2710   0.3435
    ## Detection Prevalence   0.3779   0.2710   0.3511
    ## Balanced Accuracy      0.9170   0.9551   0.8915

``` r
print("K=4")
```

    ## [1] "K=4"

``` r
conf_res4 <- get_conf_matrix("aaces.eset", "ClusterK4_kmeans", new_df, old_df)
print(conf_res4)
```

    ## Confusion Matrix and Statistics
    ## 
    ##           Reference
    ## Prediction  1  2  3  4
    ##          1 66  7  8  0
    ##          2  0 64  0  0
    ##          3  0  1 97  0
    ##          4  0  0  0 19
    ## 
    ## Overall Statistics
    ##                                           
    ##                Accuracy : 0.9389          
    ##                  95% CI : (0.9027, 0.9647)
    ##     No Information Rate : 0.4008          
    ##     P-Value [Acc > NIR] : < 2.2e-16       
    ##                                           
    ##                   Kappa : 0.9127          
    ##                                           
    ##  Mcnemar's Test P-Value : NA              
    ## 
    ## Statistics by Class:
    ## 
    ##                      Class: 1 Class: 2 Class: 3 Class: 4
    ## Sensitivity            1.0000   0.8889   0.9238  1.00000
    ## Specificity            0.9235   1.0000   0.9936  1.00000
    ## Pos Pred Value         0.8148   1.0000   0.9898  1.00000
    ## Neg Pred Value         1.0000   0.9596   0.9512  1.00000
    ## Prevalence             0.2519   0.2748   0.4008  0.07252
    ## Detection Rate         0.2519   0.2443   0.3702  0.07252
    ## Detection Prevalence   0.3092   0.2443   0.3740  0.07252
    ## Balanced Accuracy      0.9617   0.9444   0.9587  1.00000

``` r
BA_aaces_kmeans = c(conf_res2$byClass[11], conf_res3$byClass[,11], conf_res4$byClass[,11])
BA_aaces_kmeans
```

    ## Balanced Accuracy          Class: 1          Class: 2          Class: 3 
    ##         0.9000000         0.9170415         0.9551282         0.8915187 
    ##          Class: 1          Class: 2          Class: 3          Class: 4 
    ##         0.9617347         0.9444444         0.9587200         1.0000000

### AACES NMF

``` r
print("K=2")
```

    ## [1] "K=2"

``` r
conf_res2 <- get_conf_matrix("aaces.eset", "ClusterK2_NMF", new_df, old_df)
print(conf_res2)
```

    ## Confusion Matrix and Statistics
    ## 
    ##           Reference
    ## Prediction   1   2
    ##          1  78  28
    ##          2   1 155
    ##                                           
    ##                Accuracy : 0.8893          
    ##                  95% CI : (0.8449, 0.9246)
    ##     No Information Rate : 0.6985          
    ##     P-Value [Acc > NIR] : 1.649e-13       
    ##                                           
    ##                   Kappa : 0.7605          
    ##                                           
    ##  Mcnemar's Test P-Value : 1.379e-06       
    ##                                           
    ##             Sensitivity : 0.9873          
    ##             Specificity : 0.8470          
    ##          Pos Pred Value : 0.7358          
    ##          Neg Pred Value : 0.9936          
    ##              Prevalence : 0.3015          
    ##          Detection Rate : 0.2977          
    ##    Detection Prevalence : 0.4046          
    ##       Balanced Accuracy : 0.9172          
    ##                                           
    ##        'Positive' Class : 1               
    ## 

``` r
print("K=3")
```

    ## [1] "K=3"

``` r
conf_res3 <- get_conf_matrix("aaces.eset", "ClusterK3_NMF", new_df, old_df)
print(conf_res3)
```

    ## Confusion Matrix and Statistics
    ## 
    ##           Reference
    ## Prediction   1   2   3
    ##          1  16   0   0
    ##          2   3  78  18
    ##          3  16   0 131
    ## 
    ## Overall Statistics
    ##                                           
    ##                Accuracy : 0.8588          
    ##                  95% CI : (0.8106, 0.8986)
    ##     No Information Rate : 0.5687          
    ##     P-Value [Acc > NIR] : < 2.2e-16       
    ##                                           
    ##                   Kappa : 0.7479          
    ##                                           
    ##  Mcnemar's Test P-Value : 4.601e-08       
    ## 
    ## Statistics by Class:
    ## 
    ##                      Class: 1 Class: 2 Class: 3
    ## Sensitivity           0.45714   1.0000   0.8792
    ## Specificity           1.00000   0.8859   0.8584
    ## Pos Pred Value        1.00000   0.7879   0.8912
    ## Neg Pred Value        0.92276   1.0000   0.8435
    ## Prevalence            0.13359   0.2977   0.5687
    ## Detection Rate        0.06107   0.2977   0.5000
    ## Detection Prevalence  0.06107   0.3779   0.5611
    ## Balanced Accuracy     0.72857   0.9429   0.8688

``` r
print("K=4")
```

    ## [1] "K=4"

``` r
conf_res4 <- get_conf_matrix("aaces.eset", "ClusterK4_NMF", new_df, old_df)
print(conf_res4)
```

    ## Confusion Matrix and Statistics
    ## 
    ##           Reference
    ## Prediction   1   2   3   4
    ##          1  13   2  15   4
    ##          2   0  33   0   0
    ##          3   0   1 122   0
    ##          4   0   8   5  59
    ## 
    ## Overall Statistics
    ##                                           
    ##                Accuracy : 0.8664          
    ##                  95% CI : (0.8191, 0.9052)
    ##     No Information Rate : 0.542           
    ##     P-Value [Acc > NIR] : < 2.2e-16       
    ##                                           
    ##                   Kappa : 0.7951          
    ##                                           
    ##  Mcnemar's Test P-Value : 4.31e-06        
    ## 
    ## Statistics by Class:
    ## 
    ##                      Class: 1 Class: 2 Class: 3 Class: 4
    ## Sensitivity           1.00000   0.7500   0.8592   0.9365
    ## Specificity           0.91566   1.0000   0.9917   0.9347
    ## Pos Pred Value        0.38235   1.0000   0.9919   0.8194
    ## Neg Pred Value        1.00000   0.9520   0.8561   0.9789
    ## Prevalence            0.04962   0.1679   0.5420   0.2405
    ## Detection Rate        0.04962   0.1260   0.4656   0.2252
    ## Detection Prevalence  0.12977   0.1260   0.4695   0.2748
    ## Balanced Accuracy     0.95783   0.8750   0.9254   0.9356

``` r
BA_aaces_nmf = c(conf_res2$byClass[11], conf_res3$byClass[,11], conf_res4$byClass[,11])
BA_aaces_nmf
```

    ## Balanced Accuracy          Class: 1          Class: 2          Class: 3 
    ##         0.9171682         0.7285714         0.9429348         0.8688009 
    ##          Class: 1          Class: 2          Class: 3          Class: 4 
    ##         0.9578313         0.8750000         0.9254108         0.9355907

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
    ##          1 313   0
    ##          2   0 186
    ##                                      
    ##                Accuracy : 1          
    ##                  95% CI : (0.9926, 1)
    ##     No Information Rate : 0.6273     
    ##     P-Value [Acc > NIR] : < 2.2e-16  
    ##                                      
    ##                   Kappa : 1          
    ##                                      
    ##  Mcnemar's Test P-Value : NA         
    ##                                      
    ##             Sensitivity : 1.0000     
    ##             Specificity : 1.0000     
    ##          Pos Pred Value : 1.0000     
    ##          Neg Pred Value : 1.0000     
    ##              Prevalence : 0.6273     
    ##          Detection Rate : 0.6273     
    ##    Detection Prevalence : 0.6273     
    ##       Balanced Accuracy : 1.0000     
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
    ##          1 140   0   0
    ##          2   0 133   0
    ##          3   0   0 226
    ## 
    ## Overall Statistics
    ##                                      
    ##                Accuracy : 1          
    ##                  95% CI : (0.9926, 1)
    ##     No Information Rate : 0.4529     
    ##     P-Value [Acc > NIR] : < 2.2e-16  
    ##                                      
    ##                   Kappa : 1          
    ##                                      
    ##  Mcnemar's Test P-Value : NA         
    ## 
    ## Statistics by Class:
    ## 
    ##                      Class: 1 Class: 2 Class: 3
    ## Sensitivity            1.0000   1.0000   1.0000
    ## Specificity            1.0000   1.0000   1.0000
    ## Pos Pred Value         1.0000   1.0000   1.0000
    ## Neg Pred Value         1.0000   1.0000   1.0000
    ## Prevalence             0.2806   0.2665   0.4529
    ## Detection Rate         0.2806   0.2665   0.4529
    ## Detection Prevalence   0.2806   0.2665   0.4529
    ## Balanced Accuracy      1.0000   1.0000   1.0000

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
    ##          1 123   0   0   0
    ##          2   0  96   0   0
    ##          3   0   0 120   0
    ##          4   0   0   0 160
    ## 
    ## Overall Statistics
    ##                                      
    ##                Accuracy : 1          
    ##                  95% CI : (0.9926, 1)
    ##     No Information Rate : 0.3206     
    ##     P-Value [Acc > NIR] : < 2.2e-16  
    ##                                      
    ##                   Kappa : 1          
    ##                                      
    ##  Mcnemar's Test P-Value : NA         
    ## 
    ## Statistics by Class:
    ## 
    ##                      Class: 1 Class: 2 Class: 3 Class: 4
    ## Sensitivity            1.0000   1.0000   1.0000   1.0000
    ## Specificity            1.0000   1.0000   1.0000   1.0000
    ## Pos Pred Value         1.0000   1.0000   1.0000   1.0000
    ## Neg Pred Value         1.0000   1.0000   1.0000   1.0000
    ## Prevalence             0.2465   0.1924   0.2405   0.3206
    ## Detection Rate         0.2465   0.1924   0.2405   0.3206
    ## Detection Prevalence   0.2465   0.1924   0.2405   0.3206
    ## Balanced Accuracy      1.0000   1.0000   1.0000   1.0000

``` r
BA_tcga_kmeans = c(conf_res2$byClass[11], conf_res3$byClass[,11], conf_res4$byClass[,11])
BA_tcga_kmeans
```

    ## Balanced Accuracy          Class: 1          Class: 2          Class: 3 
    ##                 1                 1                 1                 1 
    ##          Class: 1          Class: 2          Class: 3          Class: 4 
    ##                 1                 1                 1                 1

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
    ##          1 247   0
    ##          2  14 238
    ##                                           
    ##                Accuracy : 0.9719          
    ##                  95% CI : (0.9534, 0.9846)
    ##     No Information Rate : 0.523           
    ##     P-Value [Acc > NIR] : < 2.2e-16       
    ##                                           
    ##                   Kappa : 0.9439          
    ##                                           
    ##  Mcnemar's Test P-Value : 0.000512        
    ##                                           
    ##             Sensitivity : 0.9464          
    ##             Specificity : 1.0000          
    ##          Pos Pred Value : 1.0000          
    ##          Neg Pred Value : 0.9444          
    ##              Prevalence : 0.5230          
    ##          Detection Rate : 0.4950          
    ##    Detection Prevalence : 0.4950          
    ##       Balanced Accuracy : 0.9732          
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
    ##          1 131   3  13
    ##          2   1 150   0
    ##          3   0  16 185
    ## 
    ## Overall Statistics
    ##                                          
    ##                Accuracy : 0.9339         
    ##                  95% CI : (0.9084, 0.954)
    ##     No Information Rate : 0.3968         
    ##     P-Value [Acc > NIR] : < 2.2e-16      
    ##                                          
    ##                   Kappa : 0.8998         
    ##                                          
    ##  Mcnemar's Test P-Value : 1.38e-06       
    ## 
    ## Statistics by Class:
    ## 
    ##                      Class: 1 Class: 2 Class: 3
    ## Sensitivity            0.9924   0.8876   0.9343
    ## Specificity            0.9564   0.9970   0.9468
    ## Pos Pred Value         0.8912   0.9934   0.9204
    ## Neg Pred Value         0.9972   0.9454   0.9564
    ## Prevalence             0.2645   0.3387   0.3968
    ## Detection Rate         0.2625   0.3006   0.3707
    ## Detection Prevalence   0.2946   0.3026   0.4028
    ## Balanced Accuracy      0.9744   0.9423   0.9406

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
    ## Prediction   1   2   3   4
    ##          1 109   0   6   1
    ##          2   4  87   6   0
    ##          3   0   5 122  49
    ##          4   7  32   0  71
    ## 
    ## Overall Statistics
    ##                                           
    ##                Accuracy : 0.7796          
    ##                  95% CI : (0.7406, 0.8152)
    ##     No Information Rate : 0.2685          
    ##     P-Value [Acc > NIR] : < 2.2e-16       
    ##                                           
    ##                   Kappa : 0.7051          
    ##                                           
    ##  Mcnemar's Test P-Value : < 2.2e-16       
    ## 
    ## Statistics by Class:
    ## 
    ##                      Class: 1 Class: 2 Class: 3 Class: 4
    ## Sensitivity            0.9083   0.7016   0.9104   0.5868
    ## Specificity            0.9815   0.9733   0.8521   0.8968
    ## Pos Pred Value         0.9397   0.8969   0.6932   0.6455
    ## Neg Pred Value         0.9713   0.9080   0.9628   0.8715
    ## Prevalence             0.2405   0.2485   0.2685   0.2425
    ## Detection Rate         0.2184   0.1743   0.2445   0.1423
    ## Detection Prevalence   0.2325   0.1944   0.3527   0.2204
    ## Balanced Accuracy      0.9449   0.8375   0.8813   0.7418

``` r
BA_tcga_nmf = c(conf_res2$byClass[11], conf_res3$byClass[,11], conf_res4$byClass[,11])
BA_tcga_nmf
```

    ## Balanced Accuracy          Class: 1          Class: 2          Class: 3 
    ##         0.9731801         0.9744138         0.9422718         0.9405936 
    ##          Class: 1          Class: 2          Class: 3          Class: 4 
    ##         0.9449318         0.8374731         0.8812513         0.7418011

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
    ##          1 190   0
    ##          2   0 187
    ##                                      
    ##                Accuracy : 1          
    ##                  95% CI : (0.9903, 1)
    ##     No Information Rate : 0.504      
    ##     P-Value [Acc > NIR] : < 2.2e-16  
    ##                                      
    ##                   Kappa : 1          
    ##                                      
    ##  Mcnemar's Test P-Value : NA         
    ##                                      
    ##             Sensitivity : 1.000      
    ##             Specificity : 1.000      
    ##          Pos Pred Value : 1.000      
    ##          Neg Pred Value : 1.000      
    ##              Prevalence : 0.504      
    ##          Detection Rate : 0.504      
    ##    Detection Prevalence : 0.504      
    ##       Balanced Accuracy : 1.000      
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
    ##          1 147   0   0
    ##          2   0  88   0
    ##          3   0   0 142
    ## 
    ## Overall Statistics
    ##                                      
    ##                Accuracy : 1          
    ##                  95% CI : (0.9903, 1)
    ##     No Information Rate : 0.3899     
    ##     P-Value [Acc > NIR] : < 2.2e-16  
    ##                                      
    ##                   Kappa : 1          
    ##                                      
    ##  Mcnemar's Test P-Value : NA         
    ## 
    ## Statistics by Class:
    ## 
    ##                      Class: 1 Class: 2 Class: 3
    ## Sensitivity            1.0000   1.0000   1.0000
    ## Specificity            1.0000   1.0000   1.0000
    ## Pos Pred Value         1.0000   1.0000   1.0000
    ## Neg Pred Value         1.0000   1.0000   1.0000
    ## Prevalence             0.3899   0.2334   0.3767
    ## Detection Rate         0.3899   0.2334   0.3767
    ## Detection Prevalence   0.3899   0.2334   0.3767
    ## Balanced Accuracy      1.0000   1.0000   1.0000

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
    ##          1 106   0   0   0
    ##          2   0  79   0   0
    ##          3   0   0  95   0
    ##          4   0   0   0  97
    ## 
    ## Overall Statistics
    ##                                      
    ##                Accuracy : 1          
    ##                  95% CI : (0.9903, 1)
    ##     No Information Rate : 0.2812     
    ##     P-Value [Acc > NIR] : < 2.2e-16  
    ##                                      
    ##                   Kappa : 1          
    ##                                      
    ##  Mcnemar's Test P-Value : NA         
    ## 
    ## Statistics by Class:
    ## 
    ##                      Class: 1 Class: 2 Class: 3 Class: 4
    ## Sensitivity            1.0000   1.0000    1.000   1.0000
    ## Specificity            1.0000   1.0000    1.000   1.0000
    ## Pos Pred Value         1.0000   1.0000    1.000   1.0000
    ## Neg Pred Value         1.0000   1.0000    1.000   1.0000
    ## Prevalence             0.2812   0.2095    0.252   0.2573
    ## Detection Rate         0.2812   0.2095    0.252   0.2573
    ## Detection Prevalence   0.2812   0.2095    0.252   0.2573
    ## Balanced Accuracy      1.0000   1.0000    1.000   1.0000

``` r
BA_mayo_kmeans = c(conf_res2$byClass[11], conf_res3$byClass[,11], conf_res4$byClass[,11])
BA_mayo_kmeans
```

    ## Balanced Accuracy          Class: 1          Class: 2          Class: 3 
    ##                 1                 1                 1                 1 
    ##          Class: 1          Class: 2          Class: 3          Class: 4 
    ##                 1                 1                 1                 1

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
    ##          2   3 227
    ##                                           
    ##                Accuracy : 0.992           
    ##                  95% CI : (0.9769, 0.9984)
    ##     No Information Rate : 0.6021          
    ##     P-Value [Acc > NIR] : <2e-16          
    ##                                           
    ##                   Kappa : 0.9833          
    ##                                           
    ##  Mcnemar's Test P-Value : 0.2482          
    ##                                           
    ##             Sensitivity : 0.9800          
    ##             Specificity : 1.0000          
    ##          Pos Pred Value : 1.0000          
    ##          Neg Pred Value : 0.9870          
    ##              Prevalence : 0.3979          
    ##          Detection Rate : 0.3899          
    ##    Detection Prevalence : 0.3899          
    ##       Balanced Accuracy : 0.9900          
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
    ## Prediction   1   2   3
    ##          1 129   1   2
    ##          2   0  98   5
    ##          3   0   0 142
    ## 
    ## Overall Statistics
    ##                                           
    ##                Accuracy : 0.9788          
    ##                  95% CI : (0.9586, 0.9908)
    ##     No Information Rate : 0.3952          
    ##     P-Value [Acc > NIR] : < 2e-16         
    ##                                           
    ##                   Kappa : 0.9678          
    ##                                           
    ##  Mcnemar's Test P-Value : 0.04601         
    ## 
    ## Statistics by Class:
    ## 
    ##                      Class: 1 Class: 2 Class: 3
    ## Sensitivity            1.0000   0.9899   0.9530
    ## Specificity            0.9879   0.9820   1.0000
    ## Pos Pred Value         0.9773   0.9515   1.0000
    ## Neg Pred Value         1.0000   0.9964   0.9702
    ## Prevalence             0.3422   0.2626   0.3952
    ## Detection Rate         0.3422   0.2599   0.3767
    ## Detection Prevalence   0.3501   0.2732   0.3767
    ## Balanced Accuracy      0.9940   0.9860   0.9765

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
    ## Prediction   1   2   3   4
    ##          1  73   0  13  15
    ##          2  12  68   6   2
    ##          3   0  11 105   4
    ##          4   0  13   0  55
    ## 
    ## Overall Statistics
    ##                                           
    ##                Accuracy : 0.7984          
    ##                  95% CI : (0.7543, 0.8378)
    ##     No Information Rate : 0.3289          
    ##     P-Value [Acc > NIR] : < 2.2e-16       
    ##                                           
    ##                   Kappa : 0.7282          
    ##                                           
    ##  Mcnemar's Test P-Value : 9.145e-10       
    ## 
    ## Statistics by Class:
    ## 
    ##                      Class: 1 Class: 2 Class: 3 Class: 4
    ## Sensitivity            0.8588   0.7391   0.8468   0.7237
    ## Specificity            0.9041   0.9298   0.9407   0.9568
    ## Pos Pred Value         0.7228   0.7727   0.8750   0.8088
    ## Neg Pred Value         0.9565   0.9170   0.9261   0.9320
    ## Prevalence             0.2255   0.2440   0.3289   0.2016
    ## Detection Rate         0.1936   0.1804   0.2785   0.1459
    ## Detection Prevalence   0.2679   0.2334   0.3183   0.1804
    ## Balanced Accuracy      0.8815   0.8345   0.8937   0.8402

``` r
BA_mayo_nmf = c(conf_res2$byClass[11], conf_res3$byClass[,11], conf_res4$byClass[,11])
BA_mayo_nmf
```

    ## Balanced Accuracy          Class: 1          Class: 2          Class: 3 
    ##         0.9900000         0.9939516         0.9859567         0.9765101 
    ##          Class: 1          Class: 2          Class: 3          Class: 4 
    ##         0.8814666         0.8344775         0.8937428         0.8402474

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
    ##          1 118   0
    ##          2   0 123
    ##                                      
    ##                Accuracy : 1          
    ##                  95% CI : (0.9848, 1)
    ##     No Information Rate : 0.5104     
    ##     P-Value [Acc > NIR] : < 2.2e-16  
    ##                                      
    ##                   Kappa : 1          
    ##                                      
    ##  Mcnemar's Test P-Value : NA         
    ##                                      
    ##             Sensitivity : 1.0000     
    ##             Specificity : 1.0000     
    ##          Pos Pred Value : 1.0000     
    ##          Neg Pred Value : 1.0000     
    ##              Prevalence : 0.4896     
    ##          Detection Rate : 0.4896     
    ##    Detection Prevalence : 0.4896     
    ##       Balanced Accuracy : 1.0000     
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
    ##          1 96  0  0
    ##          2  0 48  0
    ##          3  0  0 97
    ## 
    ## Overall Statistics
    ##                                      
    ##                Accuracy : 1          
    ##                  95% CI : (0.9848, 1)
    ##     No Information Rate : 0.4025     
    ##     P-Value [Acc > NIR] : < 2.2e-16  
    ##                                      
    ##                   Kappa : 1          
    ##                                      
    ##  Mcnemar's Test P-Value : NA         
    ## 
    ## Statistics by Class:
    ## 
    ##                      Class: 1 Class: 2 Class: 3
    ## Sensitivity            1.0000   1.0000   1.0000
    ## Specificity            1.0000   1.0000   1.0000
    ## Pos Pred Value         1.0000   1.0000   1.0000
    ## Neg Pred Value         1.0000   1.0000   1.0000
    ## Prevalence             0.3983   0.1992   0.4025
    ## Detection Rate         0.3983   0.1992   0.4025
    ## Detection Prevalence   0.3983   0.1992   0.4025
    ## Balanced Accuracy      1.0000   1.0000   1.0000

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
    ##          1 75  0  0  0
    ##          2  0 45  0  0
    ##          3  0  0 67  0
    ##          4  0  0  0 54
    ## 
    ## Overall Statistics
    ##                                      
    ##                Accuracy : 1          
    ##                  95% CI : (0.9848, 1)
    ##     No Information Rate : 0.3112     
    ##     P-Value [Acc > NIR] : < 2.2e-16  
    ##                                      
    ##                   Kappa : 1          
    ##                                      
    ##  Mcnemar's Test P-Value : NA         
    ## 
    ## Statistics by Class:
    ## 
    ##                      Class: 1 Class: 2 Class: 3 Class: 4
    ## Sensitivity            1.0000   1.0000    1.000   1.0000
    ## Specificity            1.0000   1.0000    1.000   1.0000
    ## Pos Pred Value         1.0000   1.0000    1.000   1.0000
    ## Neg Pred Value         1.0000   1.0000    1.000   1.0000
    ## Prevalence             0.3112   0.1867    0.278   0.2241
    ## Detection Rate         0.3112   0.1867    0.278   0.2241
    ## Detection Prevalence   0.3112   0.1867    0.278   0.2241
    ## Balanced Accuracy      1.0000   1.0000    1.000   1.0000

``` r
BA_tothill_kmeans = c(conf_res2$byClass[11], conf_res3$byClass[,11], conf_res4$byClass[,11])
BA_tothill_kmeans
```

    ## Balanced Accuracy          Class: 1          Class: 2          Class: 3 
    ##                 1                 1                 1                 1 
    ##          Class: 1          Class: 2          Class: 3          Class: 4 
    ##                 1                 1                 1                 1

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
    ##          1 121   6
    ##          2   0 114
    ##                                           
    ##                Accuracy : 0.9751          
    ##                  95% CI : (0.9466, 0.9908)
    ##     No Information Rate : 0.5021          
    ##     P-Value [Acc > NIR] : < 2e-16         
    ##                                           
    ##                   Kappa : 0.9502          
    ##                                           
    ##  Mcnemar's Test P-Value : 0.04123         
    ##                                           
    ##             Sensitivity : 1.0000          
    ##             Specificity : 0.9500          
    ##          Pos Pred Value : 0.9528          
    ##          Neg Pred Value : 1.0000          
    ##              Prevalence : 0.5021          
    ##          Detection Rate : 0.5021          
    ##    Detection Prevalence : 0.5270          
    ##       Balanced Accuracy : 0.9750          
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
    ##          1 78  0  4
    ##          2  2 65  0
    ##          3  0  4 88
    ## 
    ## Overall Statistics
    ##                                          
    ##                Accuracy : 0.9585         
    ##                  95% CI : (0.925, 0.9799)
    ##     No Information Rate : 0.3817         
    ##     P-Value [Acc > NIR] : < 2e-16        
    ##                                          
    ##                   Kappa : 0.9373         
    ##                                          
    ##  Mcnemar's Test P-Value : 0.01857        
    ## 
    ## Statistics by Class:
    ## 
    ##                      Class: 1 Class: 2 Class: 3
    ## Sensitivity            0.9750   0.9420   0.9565
    ## Specificity            0.9752   0.9884   0.9732
    ## Pos Pred Value         0.9512   0.9701   0.9565
    ## Neg Pred Value         0.9874   0.9770   0.9732
    ## Prevalence             0.3320   0.2863   0.3817
    ## Detection Rate         0.3237   0.2697   0.3651
    ## Detection Prevalence   0.3402   0.2780   0.3817
    ## Balanced Accuracy      0.9751   0.9652   0.9648

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
    ##          1 37  0 27  4
    ##          2  9 39  1  0
    ##          3  0  3 48  0
    ##          4  0 13  0 60
    ## 
    ## Overall Statistics
    ##                                           
    ##                Accuracy : 0.7635          
    ##                  95% CI : (0.7047, 0.8157)
    ##     No Information Rate : 0.3154          
    ##     P-Value [Acc > NIR] : < 2.2e-16       
    ##                                           
    ##                   Kappa : 0.6857          
    ##                                           
    ##  Mcnemar's Test P-Value : NA              
    ## 
    ## Statistics by Class:
    ## 
    ##                      Class: 1 Class: 2 Class: 3 Class: 4
    ## Sensitivity            0.8043   0.7091   0.6316   0.9375
    ## Specificity            0.8410   0.9462   0.9818   0.9266
    ## Pos Pred Value         0.5441   0.7959   0.9412   0.8219
    ## Neg Pred Value         0.9480   0.9167   0.8526   0.9762
    ## Prevalence             0.1909   0.2282   0.3154   0.2656
    ## Detection Rate         0.1535   0.1618   0.1992   0.2490
    ## Detection Prevalence   0.2822   0.2033   0.2116   0.3029
    ## Balanced Accuracy      0.8227   0.8277   0.8067   0.9320

``` r
BA_tothill_nmf = c(conf_res2$byClass[11], conf_res3$byClass[,11], conf_res4$byClass[,11])
BA_tothill_nmf
```

    ## Balanced Accuracy          Class: 1          Class: 2          Class: 3 
    ##         0.9750000         0.9750776         0.9652005         0.9648381 
    ##          Class: 1          Class: 2          Class: 3          Class: 4 
    ##         0.8226867         0.8276637         0.8066986         0.9320268

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
    ##          1  96   2
    ##          2   0 157
    ##                                         
    ##                Accuracy : 0.9922        
    ##                  95% CI : (0.972, 0.999)
    ##     No Information Rate : 0.6235        
    ##     P-Value [Acc > NIR] : <2e-16        
    ##                                         
    ##                   Kappa : 0.9834        
    ##                                         
    ##  Mcnemar's Test P-Value : 0.4795        
    ##                                         
    ##             Sensitivity : 1.0000        
    ##             Specificity : 0.9874        
    ##          Pos Pred Value : 0.9796        
    ##          Neg Pred Value : 1.0000        
    ##              Prevalence : 0.3765        
    ##          Detection Rate : 0.3765        
    ##    Detection Prevalence : 0.3843        
    ##       Balanced Accuracy : 0.9937        
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
    ##          1  93   0   1
    ##          2   0  38   0
    ##          3   0   0 123
    ## 
    ## Overall Statistics
    ##                                           
    ##                Accuracy : 0.9961          
    ##                  95% CI : (0.9783, 0.9999)
    ##     No Information Rate : 0.4863          
    ##     P-Value [Acc > NIR] : < 2.2e-16       
    ##                                           
    ##                   Kappa : 0.9936          
    ##                                           
    ##  Mcnemar's Test P-Value : NA              
    ## 
    ## Statistics by Class:
    ## 
    ##                      Class: 1 Class: 2 Class: 3
    ## Sensitivity            1.0000    1.000   0.9919
    ## Specificity            0.9938    1.000   1.0000
    ## Pos Pred Value         0.9894    1.000   1.0000
    ## Neg Pred Value         1.0000    1.000   0.9924
    ## Prevalence             0.3647    0.149   0.4863
    ## Detection Rate         0.3647    0.149   0.4824
    ## Detection Prevalence   0.3686    0.149   0.4824
    ## Balanced Accuracy      0.9969    1.000   0.9960

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
    ##          1 88  0  0  0
    ##          2  0 29  0  0
    ##          3  0  0 73  0
    ##          4  0  0  0 65
    ## 
    ## Overall Statistics
    ##                                      
    ##                Accuracy : 1          
    ##                  95% CI : (0.9856, 1)
    ##     No Information Rate : 0.3451     
    ##     P-Value [Acc > NIR] : < 2.2e-16  
    ##                                      
    ##                   Kappa : 1          
    ##                                      
    ##  Mcnemar's Test P-Value : NA         
    ## 
    ## Statistics by Class:
    ## 
    ##                      Class: 1 Class: 2 Class: 3 Class: 4
    ## Sensitivity            1.0000   1.0000   1.0000   1.0000
    ## Specificity            1.0000   1.0000   1.0000   1.0000
    ## Pos Pred Value         1.0000   1.0000   1.0000   1.0000
    ## Neg Pred Value         1.0000   1.0000   1.0000   1.0000
    ## Prevalence             0.3451   0.1137   0.2863   0.2549
    ## Detection Rate         0.3451   0.1137   0.2863   0.2549
    ## Detection Prevalence   0.3451   0.1137   0.2863   0.2549
    ## Balanced Accuracy      1.0000   1.0000   1.0000   1.0000

``` r
BA_yoshihara_kmeans = c(conf_res2$byClass[11], conf_res3$byClass[,11], conf_res4$byClass[,11])
BA_yoshihara_kmeans
```

    ## Balanced Accuracy          Class: 1          Class: 2          Class: 3 
    ##         0.9937107         0.9969136         1.0000000         0.9959677 
    ##          Class: 1          Class: 2          Class: 3          Class: 4 
    ##         1.0000000         1.0000000         1.0000000         1.0000000

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
    ##          1 100   2
    ##          2   0 153
    ##                                         
    ##                Accuracy : 0.9922        
    ##                  95% CI : (0.972, 0.999)
    ##     No Information Rate : 0.6078        
    ##     P-Value [Acc > NIR] : <2e-16        
    ##                                         
    ##                   Kappa : 0.9836        
    ##                                         
    ##  Mcnemar's Test P-Value : 0.4795        
    ##                                         
    ##             Sensitivity : 1.0000        
    ##             Specificity : 0.9871        
    ##          Pos Pred Value : 0.9804        
    ##          Neg Pred Value : 1.0000        
    ##              Prevalence : 0.3922        
    ##          Detection Rate : 0.3922        
    ##    Detection Prevalence : 0.4000        
    ##       Balanced Accuracy : 0.9935        
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
    ##          1  62   2   0
    ##          2   0  68   0
    ##          3   1   7 115
    ## 
    ## Overall Statistics
    ##                                          
    ##                Accuracy : 0.9608         
    ##                  95% CI : (0.9291, 0.981)
    ##     No Information Rate : 0.451          
    ##     P-Value [Acc > NIR] : < 2e-16        
    ##                                          
    ##                   Kappa : 0.9387         
    ##                                          
    ##  Mcnemar's Test P-Value : 0.01857        
    ## 
    ## Statistics by Class:
    ## 
    ##                      Class: 1 Class: 2 Class: 3
    ## Sensitivity            0.9841   0.8831   1.0000
    ## Specificity            0.9896   1.0000   0.9429
    ## Pos Pred Value         0.9688   1.0000   0.9350
    ## Neg Pred Value         0.9948   0.9519   1.0000
    ## Prevalence             0.2471   0.3020   0.4510
    ## Detection Rate         0.2431   0.2667   0.4510
    ## Detection Prevalence   0.2510   0.2667   0.4824
    ## Balanced Accuracy      0.9869   0.9416   0.9714

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
    ## Prediction   1   2   3   4
    ##          1  57   6   2   0
    ##          2   0  43   1   0
    ##          3   0   3 102   0
    ##          4   0   1   9  31
    ## 
    ## Overall Statistics
    ##                                           
    ##                Accuracy : 0.9137          
    ##                  95% CI : (0.8723, 0.9451)
    ##     No Information Rate : 0.4471          
    ##     P-Value [Acc > NIR] : < 2.2e-16       
    ##                                           
    ##                   Kappa : 0.8774          
    ##                                           
    ##  Mcnemar's Test P-Value : NA              
    ## 
    ## Statistics by Class:
    ## 
    ##                      Class: 1 Class: 2 Class: 3 Class: 4
    ## Sensitivity            1.0000   0.8113   0.8947   1.0000
    ## Specificity            0.9596   0.9950   0.9787   0.9554
    ## Pos Pred Value         0.8769   0.9773   0.9714   0.7561
    ## Neg Pred Value         1.0000   0.9526   0.9200   1.0000
    ## Prevalence             0.2235   0.2078   0.4471   0.1216
    ## Detection Rate         0.2235   0.1686   0.4000   0.1216
    ## Detection Prevalence   0.2549   0.1725   0.4118   0.1608
    ## Balanced Accuracy      0.9798   0.9032   0.9367   0.9777

``` r
BA_yoshihara_nmf = c(conf_res2$byClass[11], conf_res3$byClass[,11], conf_res4$byClass[,11])
BA_yoshihara_nmf
```

    ## Balanced Accuracy          Class: 1          Class: 2          Class: 3 
    ##         0.9935484         0.9868552         0.9415584         0.9714286 
    ##          Class: 1          Class: 2          Class: 3          Class: 4 
    ##         0.9797980         0.9031851         0.9367301         0.9776786

# Balanced Accuracy Summary

``` r
class_ids = c("K=2", "K=3, 1", "K=3, 2", "K=3, 3", "K=4, 1", "K=4, 2", "K=4, 3", "K=4, 4")


# Balanced Accuracy for Kmeans
ba_df = data.frame(class_ids,
                    BA_aaces_kmeans,
                    BA_tcga_kmeans,
                    BA_mayo_kmeans,
                    BA_tothill_kmeans,
                    BA_yoshihara_kmeans)
ba_df
```

    ##   class_ids BA_aaces_kmeans BA_tcga_kmeans BA_mayo_kmeans BA_tothill_kmeans
    ## 1       K=2       0.9000000              1              1                 1
    ## 2    K=3, 1       0.9170415              1              1                 1
    ## 3    K=3, 2       0.9551282              1              1                 1
    ## 4    K=3, 3       0.8915187              1              1                 1
    ## 5    K=4, 1       0.9617347              1              1                 1
    ## 6    K=4, 2       0.9444444              1              1                 1
    ## 7    K=4, 3       0.9587200              1              1                 1
    ## 8    K=4, 4       1.0000000              1              1                 1
    ##   BA_yoshihara_kmeans
    ## 1           0.9937107
    ## 2           0.9969136
    ## 3           1.0000000
    ## 4           0.9959677
    ## 5           1.0000000
    ## 6           1.0000000
    ## 7           1.0000000
    ## 8           1.0000000

``` r
# Balanced Accuracy for NMF
ba_df = data.frame(class_ids,
                    BA_aaces_nmf,
                    BA_tcga_nmf,
                    BA_mayo_nmf,
                    BA_tothill_nmf,
                    BA_yoshihara_nmf)
ba_df
```

    ##   class_ids BA_aaces_nmf BA_tcga_nmf BA_mayo_nmf BA_tothill_nmf
    ## 1       K=2    0.9171682   0.9731801   0.9900000      0.9750000
    ## 2    K=3, 1    0.7285714   0.9744138   0.9939516      0.9750776
    ## 3    K=3, 2    0.9429348   0.9422718   0.9859567      0.9652005
    ## 4    K=3, 3    0.8688009   0.9405936   0.9765101      0.9648381
    ## 5    K=4, 1    0.9578313   0.9449318   0.8814666      0.8226867
    ## 6    K=4, 2    0.8750000   0.8374731   0.8344775      0.8276637
    ## 7    K=4, 3    0.9254108   0.8812513   0.8937428      0.8066986
    ## 8    K=4, 4    0.9355907   0.7418011   0.8402474      0.9320268
    ##   BA_yoshihara_nmf
    ## 1        0.9935484
    ## 2        0.9868552
    ## 3        0.9415584
    ## 4        0.9714286
    ## 5        0.9797980
    ## 6        0.9031851
    ## 7        0.9367301
    ## 8        0.9776786

# Conclusion

We find that overall the clusters align between the previously published
Way et. al. run and my run of the pipeline. NOTE: There exist slight
data and possibly code differences between our two runs, so we do not
expect 100% matching between runs. However, we do expect high coherence
between runs. We find that the clusters are more consistent in K-Means
than NMF. We also find that K=4 is the most inconsistent between the two
runs of the pipeline. We will consider K=4 and KMeans results for
remaining main analyses.
