# Recreating Results From The Paper

This document provides steps to recreate the results for the two case studies, namely Anesthesia and ACC.

## Anesthesia

### Recreating Figure 3

![Fig3](Fig3.png)

Note that Figure 3 contains four figures within it. Let the bottom-left be figure number 1, bottom-right be figure number 2, top-left be figure number 3, and top-right be figure number 4. Following are the steps to recreate these figures.

#### Step 1: Go to the corresponding folder.

```shell
cd /my/location/MoULDyS/src/recreate_results_from_paper/
```

#### Step 2: Recreate the results.

To recreate results as per figure number `f`, use `python python Anesthesia.py -offline f`. For example, to recreate figure number 1:

```shell
python python Anesthesia.py -offline f
```

### Recreating Figure 4 

![Fig4](Fig4.png)

* To recreate figure 4 left:

```shell
python Anesthesia.py -online
```

* To recreate figure 4 right:

```shell
python Anesthesia.py -compare
```

## Adaptive Cruise Control (ACC)

ToDo

