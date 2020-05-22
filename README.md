# Reliable Domain Adaptation(RDA Demo)
This is the source code of the RDA method, which has been published as “Reliable Domain Adaptation with Classiﬁers Competition for Image Classiﬁcation“ IEEE Trans.CAS II Express Briefs, 2020. 

## In the file, there are multiple .m files.

`RDA.m : Core codes of RDA algorithm.`

`demo_RDA.m : Evaluate RDA/NRDA on an example task (C1-C2 in COIL20 dataset).`

## You can run the "demo_RDA.m" code for your reference. Then you will get the results:

`RDA: Choose 'primal' as 'options.kernel', the final accrucy will be ‘93.75%’.`

`NRDA: Choose 'linear'(linear kernel in our paper,for example) as 'options.kernel' , the final accrucy will be ‘96.25%’.`

## When you use our code, there are some parameters need to be adjusted according to different tasks. We recommend that you adjust these parameters when using this code on a new dataset.

## Once you run the code, please correctly set the path of the data and liblinear toolbox

## Please contact leizhang@cqu.edu.cn or jrfu@cqu.edu.cn if there is any problem

## Enjoy it!
