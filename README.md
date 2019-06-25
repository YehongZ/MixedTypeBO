# MixedTypeBO
This is the code for the following paper:

"Bayesian Optimization with Binary Auxiliary Information." Yehong Zhang, Zhongxiang Dai and Bryan Kian Hsiang Low. In Proceedings of the 35th Conference on Uncertainty in Artificial Intelligence (UAI), 2019

## Dependencies

This code depends on the following toolboxes:

* [gp](https://github.com/SheffieldML/GPmat/tree/master/gp)
* [mltools](https://github.com/SheffieldML/GPmat/tree/master/mltools)
* [ndlutil](https://github.com/SheffieldML/GPmat/tree/master/ndlutil)
* [optimi](https://github.com/SheffieldML/GPmat/tree/master/optimi)
* [multigp](https://github.com/SheffieldML/multigp)
* [netlab](https://github.com/sods/netlab)

## Examples

### Mixed-type PES using synthetic functions

This example shows how to do BO with mixed-type PES using the synthetic functions. The performance of the tested algorithms are evaluated using ten groups (i.e., one target function and one auxiliary function) of synthetic functions whose images are shown in [images](https://github.com/YehongZ/MixedTypeBO/tree/master/examples/synthetic/images). An averaged immediate regret (IR) is obtained by optimizing the target function in each of them with 10 different initializations.

```
>> run_syn
```

Compute the averaged IR and plot the result:

```
>> plot_syn
```

### Hyper-parameters tuning of a convolutional neural network (CNN) with CIFAR-10 dataset

Download the Bayesian optimal stopping (bos) function "bos_function.py" from [BO-BOS](https://github.com/daizhongxiang/Bayesian-Optimization-Meets-Bayesian-Optimal-Stopping) and put it in the CNN folder. Set the parameters in "bos_function.py" as: fs_sample_number, K1, K2, C = 10000, 100, 100, 1.

Prepare the dataset:

```
>> prepare_cifar10
```

Run the CNN example over 5 different initializations:

```
>> run_cnn
```
