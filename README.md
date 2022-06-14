# `MoULDyS`: Monitoring of scattered uncertain logs using uncertain linear dynamical systems

Monitoring the correctness of distributed cyber-physical systems is essential. We address the analysis of the log of a black-box cyber-physical system. Detecting possible safety violations can be hard when some samples are uncertain or missing. In this work, the log is made of values known with some uncertainty; in addition, we make use of an over-approximated yet expressive model, given by a non-linear extension of dynamical systems. Given an offline log, our approach is able to monitor the log against safety specifications with a limited number of false alarms. As a second contribution, we show that our approach can be used online to minimize the number of sample triggers, with the aim at energetic efficiency. 

![mouldys](mouldys.png)

## Installation

A detailed installation guide is provided in `/documentation/installation_guide.md`.

### Dependencies

- [`Python 3.9.x`](https://www.python.org/)
- [`NumPy`](https://numpy.org/)
- [`SciPy`](https://scipy.org/)
- [`mpmath`](https://mpmath.org/)
- Gurobi Python Interface:
  - Please obtain appropriate Gurobi License from [here](http://www.gurobi.com/downloads/licenses/license-center). Please refer to this [link](https://www.gurobi.com/documentation/8.1/quickstart_windows/academic_validation.html) for details. After the license is installed properly, Gurobi can be used from home network.
  - Install Gurobi. Please note that we will need Gurobi Python Interface: 
    - On-line documentation on installation can be found [here](http://www.gurobi.com/documentation/).
    - **[Recommend]** Gurobi Python Interface can also be installed through [Anaconda](https://www.anaconda.com/). Details on installing Gurobi Python Interface through `conda` can be found [here](https://www.gurobi.com/documentation/8.1/quickstart_mac/installing_the_anaconda_py.html#section:Anaconda).

### Downloading the tool

1. Download the repository to your desired location `/my/location/`:

2. Once the repository is downloaded, please open `~/.bashrc`, and add the line `export MNTR_ROOT_DIR=/my/location/MoULDyS/`, mentioned in the following steps:

   1. ```shell
      vi ~/.baschrc
      ```

   2. Once `.bashrc` is opened, please add the location, where the tool was downloaded, to a path variable `MNTR_ROOT_DIR` (This step is crucial to run the tool):

      1. ```shell
         export MNTR_ROOT_DIR=/my/location/MoULDyS/
         ```

## Running the tool

Once the dependencies are installed properly, and the path variable is set, following steps should run without any error.

Following are some of the crucial functionalities offered by this prototype tool:

1. Given uncertain logs, perform offline monitoring.
2. Online monitoring

### Case studies

We offer to case studies:

1. [Anesthesia](https://cps-vo.org/node/12111).
2. [Adaptive Cruise Control (ACC)](https://ieeexplore.ieee.org/document/7349170).

Here, we illustrate the Anesthesia case study, as the other one can be run in similar fashion.

1. ```shell
   cd src/
   ```

2. ```shell
   python Anasthesia.py
   ```

```shell
>>STATUS: Time Taken:  0.13864684104919434
>>STATUS: Log generated!
>> STATUS: Performing online monitoring using reachable sets . . .
	>> SUBSTATUS: Triggering a log at time  49 	 Time Taken (Diff):  4.118176698684692
	>> SUBSTATUS: Triggering a log at time  84 	 Time Taken (Diff):  2.3078575134277344
	.
	.
	.
	>> SUBSTATUS: Triggering a log at time  1978 	 Time Taken (Diff):  0.8975484371185303
	>> SUBSTATUS: Triggering a log at time  1997 	 Time Taken (Diff):  0.8754169940948486
	>> Safety: Safe
	>> Number of Logs:  84
	>> Time Taken:  108.90439081192017
>> STATUS: Performed online monitoring using reachable sets!
>> STATUS: Computing reachable sets as per monitors . . .
	>> SUBSTATUS: Log Step:  0 / 94 	Time Diff:  5
	>> SUBSTATUS: Log Step:  1 / 94 	Time Diff:  45
		>> SUBSTATUS: Computing refinement at at time step  41
		.
		.
		.
		>> SUBSTATUS: Computing refinement at at time step  215
	>> SUBSTATUS: Log Step:  9 / 94 	Time Diff:  25
	.
	.
	.
	>> SUBSTATUS: Log Step:  21 / 94 	Time Diff:  133
		>> SUBSTATUS: Computing refinement at at time step  439
	>>Safety: Unsafe at log step  21
	>>Time Taken:  41.25862264633179
>>STATUS: Visualizing monitors . . .
```



One should see the following plot in `/my/location/Monitoring/output/`:

![viz_compare_monitors_cp_4761-72](viz_compare_monitors_cp_4761-72.png)

Note: Visualization takes quite some time.

The results provided in the draft can be found in `/my/location/Monitoring/output/Draft`.
