# Configuration profile for the Dardel HPC system

This directory contains a configuration profile for the Dardel HPC system at
PDC.

Read more about the [available partitions](https://www.pdc.kth.se/support/documents/run_jobs/job_scheduling.html#dardel-partitions) 
and [compute nodes](https://www.pdc.kth.se/support/documents/run_jobs/job_scheduling.html#dardel-compute-nodes).

## Using the profile

To use this profile, first update the `config.yaml` file with your compute
account ID under `default-resources:`. For example, if your compute account is
`naiss2024-1-100` the file should look like this:
    
```yaml
default-resources: 
  slurm_account: naiss2024-1-100
  slurm_partition: shared
  runtime: 240
  tasks: f"{threads}"
  cpus_per_task: 1
```

To use this profile you run Snakemake from the repository root and add
`--profile dardel` to the command.

```bash
snakemake --configfile <path-to-your-configfile.yml> --profile dardel <additional-arguments>
```

## Additional tips for running on Dardel

### Using Apptainer

The `happ` workflow can be run with Apptainer to handle software dependencies.
To use Apptainer on Dardel you first need to load the Apptainer module with:

```bash
module load PDC apptainer
```

Also, make sure to set the `APPTAINER_CACHE` environment variable to a directory
somewhere in your scratch space to avoid Disk quota errors. For example, you can create a folder `$TMPDIR/.cache/apptainer` and set the environment variable like this:

```bash
mkdir -p $TMPDIR/.cache/apptainer
export APPTAINER_CACHE=$TMPDIR/.cache/apptainer
```

You can also add the `export APPTAINER_CACHE=$TMPDIR/.cache/apptainer` line to your shell configuration file (e.g. `~/.bashrc` if using bash or `~/.zshrc` if using zsh) to make the change permanent.

Then add the `--software-deployment-method apptainer` (or `--sdm apptainer`) flag to your Snakemake command:

```bash
snakemake --configfile <path-to-your-configfile.yml> --sdm apptainer --profile dardel  <additional-arguments>
```

### Using Conda

The Dardel system has a module for Conda that you can load with the following command:

```bash
module load bioinfo-tools conda
```

If this is your first time loading conda, also initialize your shell with:

```bash
conda init
```

then close and reopen your terminal.

To avoid Disk quota errors create a `pkgs` directory in your scratch path and
configure Conda to use this for storing package downloads. Run the following:

```bash
mkdir -p $TMPDIR/pkgs
conda config --add pkgs_dirs $TMPDIR/pkgs
```

It is also recommended to create the main `happ` environment in the root of the
git repository, for instance if you cloned `happ` into a directory called
`/cfs/klemming/projects/snic/snic2024-1-100/happ` you first make a subdirectory
`envs/` there and then create the environment with:

```bash
conda env create -f environment.yml -p envs/happ
```

You can then update your Conda config to make this environment findable:

```bash
conda config --append envs_dirs /cfs/klemming/projects/snic/snic2024-1-100/happ/envs
```

Once you have completed the above, on subsequent logins you load conda then activate the environment with:

```bash
module load bioinfo-tools conda
conda activate happ
```

### Using pixi

If you want to use [pixi](https://pixi.sh/) to handle software environments you
will need to download and install pixi on Dardel. You can do this by running

```bash
curl -fsSL https://pixi.sh/install.sh | bash
```

You may have to close and reopen your shell after this to get the `pixi`
command to work.

To make sure you don't run into Disk quota errors you should configure pixi to
use your scratch path for storing its cache. Add the following line to your
shell configuration file (e.g. `~/.bashrc` if using bash or `~/.zshrc` if using zsh):

```bash
export PIXI_CACHE_DIR="$TMPDIR/.cache/pixi"
````

After adding this line you need to source the file to apply the changes:

```bash
source ~/.bashrc
```

Also, create the cache directory path by running:

```bash
mkdir -p $PIXI_CACHE_DIR
```

### Using terminal multiplexers

It is highly recommended to use a terminal multiplexer such as
[Screen](https://www.gnu.org/software/screen/) or
[tmux](https://github.com/tmux/tmux) when working on Dardel. This allows you to
keep your session alive even if you get disconnected from the system.

Both these programs are available on Dardel. To use `tmux` first run `module
load tmux`. `screen` is available by default.

To start a new session with `tmux` simply run:

```bash
tmux
```

then work as normal (see this [cheat sheet](https://tmuxcheatsheet.com/) for
some useful commands). If you get disconnected from the system you can simply
reconnect and then run `tmux attach` to reattach to your session without losing
any work.

To start a new session with `screen` you can run:

```bash
screen -S mysession
```

which will start a session called `mysession`. See the
[documentation](https://www.gnu.org/software/screen/) for more information on
how to use `screen`. If you get disconnected from the system, reconnect and
re-attach to your session with:

```bash
screen -x mysession
```