# SubGesture library

### Author:
###### Víctor Ponce López.
###### University of Barcelona, Open University of Catalonia, Computer Vision Center.
######© 2015. All rights reserved ©


##### Please, cite the following publication if you make any usage of this code

V. Ponce-López, H.J. Escalante, S. Escalera, X. Baró. Gesture and Action Recognition by Evolved Dynamic Subgestures. Proceedings of the British Machine Vision Conference (BMVC), pp. 129.1-129.13, 2015.

```BIB
@inproceedings{BMVC2015_129,
	title={Gesture and Action Recognition by Evolved Dynamic Subgestures},
	author={Víctor Ponce-López and Hugo Jair Escalante and Sergio Escalera and Xavier Baró},
	year={2015},
	month={September},
	pages={129.1-129.13},
	articleno={129},
	numpages={13},
	booktitle={Proceedings of the British Machine Vision Conference (BMVC)},
	publisher={BMVA Press},
	editor={Xianghua Xie, Mark W. Jones, and Gary K. L. Tam},
	isbn={1-901725-53-7},
}
```
##### Thanks for using this code.

### Main functions:

##### testDTWalign:

Simple tests with the implementations both for the DTW computation and for the mean aligment.

##### testHMM:

Test the HMM setting to discretize the gesture samples, learn a HMM model for each class and evaluate the models.

##### runGenDTW:

Main file to run the whole evolutionary algorithm with temporal segmentation over training data, using the parameter settings from the file `varload.m`.