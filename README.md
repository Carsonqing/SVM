# Goal:  
Write python code to solve **LS-SVM**, **LS-OCSVM** and **LS-SVDD** problems.  

# Dataset:
* [`charlie.csv`](https://drive.google.com/file/d/11gKBiIGYTEWPNmXZZU8L8EXT0dNdl2GO/view?usp=sharing)
* The response has two classes: "original" and "new"
* x1 to x4 are predictors

# Note:
* In one class case, use "original" for training and "new" for testing
* No svm packages from _python_ are used. Matrix inversion solutions for alpha and b are utilized to solve the problems
