### Scenario 1 variation in event fraction
Constant:
-	AUC = 0.75
-	Dimensionality = 10
Varying:
-	Sample size: exactly at and below minimal required sample size
-	Event-fraction: 0.05 & 0.2
-	Prediction models: Elastic net regression + SVM (Gaussian kernel? RBF kernel?)
7 internal validation methods

Gives: 2 * 2 * (2 + 1) = 8 different scenarios
8 * 7 (i.v. approaches) = 56

### Scenario 2 variation in dimensionality
Constant:
-	AUC @ 0.75
-	Event fraction @ 0.2
Varying:
-	Sample size: exactly at and below minimal required sample size
-	Dimensionality: 6 - 33 - 60
-	Prediction models: Elastic net regression (or OLS?) + SVM (+ gaussian kernel/rbf?)
7 internal validation methods

Gives: 2 * 3 * 2 = 12 different scenarios
12 * 7 (i.v. approaches) = 84
 
### Scenario 3: variation in prediction models
Constant:
-	AUC @ 0.75
-	Event fraction @ 0.2
-	Dimensionality @ 10 (?)
Varying:
-	Sample size: exactly at and below minimal required sample size
-	Prediction models:
*  OLS
	Ridge
	Lasso
	Elastic net
	SVM
	ANN
	RF
7 internal validation methods

Gives: 2 * 7  = 14 different scenarios
14 * 7 (i.v. approaches) = 98

### Scenario 4:
Overfitting  calibratie op .6 of .3
 
