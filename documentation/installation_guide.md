# `MoULDyS` Installation Guide

**The tool requires a Linux environment**. 

`MoULDyS` can be used in the following two ways:

1. **Virtual Machine Image (Recommend)**. This is the simplest way to use `MoULDyS`, which does not necessitate the installation of any dependencies or code downloading. Nevertheless, it is required to acquire and install the [Gurobi](https://www.gurobi.com/solutions/gurobi-optimizer/?campaignid=193283256&adgroupid=138872523040&creative=596136082776&keyword=gurobi&matchtype=e&gclid=CjwKCAjw6IiiBhAOEiwALNqncXIGRe-OYdzuBIwq3Waarc4fe6rP6DRYPh1xTWfA86OQSH_oX5zbdRoC7IUQAvD_BwE) license. Users can recreate the results easily using this method. This also requires [VirtualBox](https://www.virtualbox.org/) installed on the user's machine.
2. **Install `MoULDyS` on Local Machine**. This option requires installation of the tool from scratch. 

## **Virtual Machine Image (Recommend)**

1. If [VirtualBox](https://www.virtualbox.org/) is not already installed, please install the appropriate platform package from [here](https://www.virtualbox.org/wiki/Downloads).

2. Once the VirtualBox is installed, the VM image of `MoULDyS` can be downloaded from [here](..). Once downloaded, the image must be loaded in the VirtualBox. It should startup a Ubuntu virtual machine (VM), which comes preloaded with `MoULDyS`.

3. Once the Ubuntu VM is started, please login to the `admin` account using password `mouldys123`.

4. Please be advised that we are presently offering a license, thus omitting this step is feasible. Nonetheless, this arrangement might be terminated in the near future, and we strongly recommend users to proceed to the next step and verify if the license is activated. In case it is not, kindly ensure to carry out this step.

   1. Once logged in, please obtain appropriate Gurobi License from [here](http://www.gurobi.com/downloads/licenses/license-center). After the license is installed properly, Gurobi can be used from home network.

   2. Though one should use a license that is appropriate for them and their organization, it is worth point out that Gurobi offers free academic licenses. Here, we mention the steps to obtain a free a license. Obtain your free license by following the instructions [here](https://www.gurobi.com/academia/academic-program-and-licenses/) (please select `Individual Academic Licenses`). The license can be installed as follows (**note: must be in your university network, or in VPN**):

      * ```shell
        grbgetkey <your-license-key>
        ```

      * **Note:** Gurobi doesn't allow the same license to be used on two different computers. Please see the details [here](https://www.gurobi.com/downloads/end-user-license-agreement-academic/) (especially if you want to use it on two different computers).

5. We have provided a testing script, [`testInstall.py`](https://github.com/bineet-coderep/MoULDyS/blob/main/env_test/testInstall.py), in the folder [`/my/location/MoULDyS/env_test/`](https://github.com/bineet-coderep/MoULDyS/tree/main/env_test) to check if the environment is ready. To perform the testing, please follow the following steps.

   1. One can simply test if their environment is ready by running the following script:

      * ```shell
        python /home/MoULDyS/env_test/testInstall.py
        ```

   2. If the following output message displays in the console (in cyan color), the environment is ready:

      * ```shell
        =======================
        Environment is Ready!
        =======================
        ```

      * Following is a screenshot displaying a ready environment: Look for "Environment is Ready!" in cyan.

        * ![test_env_op](test_env_op.png)

   3. If any other error message pops up, the environment is most likely not ready.

6. One can now recreate the results following the instructions in [`documentation/recreate_results.md`](https://github.com/bineet-coderep/MoULDyS/blob/main/documentation/recreate_results.md).

## Install `MoULDyS` on Local Machine

`MoULDyS` can be installed on a local machine by performing the following steps

### Install Dependencies

One needs to install the following dependencies first. The `debian` package names are provided in brackets.

- [`Python 3.7.x`](https://www.python.org/).

- [`NumPy`](https://numpy.org/)([`python-numpy`](https://packages.debian.org/search?keywords=python-numpy)).

- [`SciPy`](https://scipy.org/)([`python-scipy`](https://packages.debian.org/search?keywords=python-scipy)).

- [`mpmath`](https://mpmath.org/)([`python3-mpmath`](https://packages.debian.org/search?keywords=python3-mpmath)).

- [`pandas`](https://pandas.pydata.org/)([`python-pandas`](https://packages.debian.org/search?suite=default&section=all&arch=any&searchon=names&keywords=python-pandas)).

- Gurobi Python Interface:
  1. Install Gurobi. Please note that we will need Gurobi Python Interface. On-line documentation on installation can be found [here](https://www.gurobi.com/documentation/9.5/quickstart_linux/cs_using_pip_to_install_gr.html). 

     1. **[If you are NOT a `conda` user]** One can use `pip` to install Gurobi Python interface as follows:

        * ```shell
          python -m pip install gurobipy
          ```

        * Make sure to have upgraded `numpy`. If not already done, one can issue the following command:

           * ```shell
             pip install --upgrade numpy
             ```

        * if using `pip`, one has to install `grbgetkey` manually (see [here](https://support.gurobi.com/hc/en-us/articles/360059842732)).

     2. **[If you are a `conda` user]** Gurobi Python Interface can also be installed through [Anaconda](https://www.anaconda.com/). Details on installing Gurobi Python Interface through `conda` can be found [here](https://www.gurobi.com/documentation/9.5/quickstart_mac/cs_anaconda_and_grb_conda_.html). One can use the following steps to install:

        1. ```shell
           conda config --add channels https://conda.anaconda.org/gurobi
           ```

        2. ```shell
           conda install gurobi
           ```

  2. Please obtain appropriate Gurobi License from [here](http://www.gurobi.com/downloads/licenses/license-center). After the license is installed properly, Gurobi can be used from home network.

     * Though one should use a license that is appropriate for them and their organization, it is worth point out that Gurobi offers free academic licenses. Here, we mention the steps to obtain a free a license. Obtain your free license by following the instructions [here](https://www.gurobi.com/academia/academic-program-and-licenses/) (please select `Individual Academic Licenses`). The license can be installed as follows (**note: must be in your university network, or in VPN**):

       * ```shell
         grbgetkey <your-license-key>
         ```

       * **Note:** Gurobi doesn't allow the same license to be used on two different computers. Please see the details [here](https://www.gurobi.com/downloads/end-user-license-agreement-academic/) (especially if you want to use it on two different computers).

### Downloading and Setting-Up `MoULDyS`

1. Download the repository to your desired location `/my/location/`:

   * ```shell
     git clone https://github.com/bineet-coderep/MoULDyS.git
     ```

2. Once the repository is downloaded, please perform the following steps:

   1. Open `~/.bashrc` using your choice of editor (say, `vi`):

     * ```shell
       vi ~/.baschrc
       ```

   2. Once `.bashrc` is opened, please add the location of the tool to a path variable `MNTR_ROOT_DIR` (This step is crucial to run the tool):

     * ```shell
       export MNTR_ROOT_DIR=/my/location/MoULDyS/
       ```

### Testing The Environment (Optional)

We have provided a testing script, [`testInstall.py`](https://github.com/bineet-coderep/MoULDyS/blob/main/env_test/testInstall.py), in the folder [`/my/location/MoULDyS/env_test/`](https://github.com/bineet-coderep/MoULDyS/tree/main/env_test) to check if the environment is ready. To perform the testing, please follow the following steps.

1. One can simply test if their environment is ready by running the following script:

   * ```shell
     python /my/location/MoULDyS/env_test/testInstall.py
     ```

2. If the following output message displays in the console (in cyan color), the environment is ready:

   * ```shell
     =======================
     Environment is Ready!
     =======================
     ```

   * Following is a screenshot displaying a ready environment: Look for "Environment is Ready!" in cyan.

     * ![test_env_op](test_env_op.png)

3. If any other error message pops up, the environment is most likely not ready.

## Development Environment Details (Optional)

Following environment has been used for development and testing `MoULDyS`.

### System Details of The Development Environment 

* OS Name: `Ubuntu 20.04.4 LTS`.
* OS Type: `64 bit`.
* GNOME Version: `3.36.8`.
* Windowing System: `X11`.

### Hardware Details of The Development Environment 

* Model: Alienware Area 51m R2.
* Processor: `Intel® Core™ i7-10700 CPU @ 2.90GHz × 16 `.
* Memory: `31.1 GiB`.
* Graphics: `NVIDIA RTX 2070 Super 8GB GDDR6 Dual `.
* Disk Capacity: `4.0 TB`.

Note: This is just the development platform, not a requirement for the tool to work.
