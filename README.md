# StatInference

## How to run binning optimisation on lxplus

### Server side
Better to run on screen.

Example:
```sh
source env.sh lcg
python tools/optimize_channel.py --input /afs/cern.ch/user/f/fbrivio/public/datacards_5July2021_1000bins_DYmerged --channel tauTau_2018 --output output/binning_v2
```

### Worker side
Running on screen is not required, jobs are submitted on HTCondor.

Example:
```sh
source env.sh hh
python tools/submitLimitWorkers.py --output output/binning_v2 --channel tauTau_2018 --n-workers 8
```

If runtime is longer than one day, both server and workers should be restarted.
Workers directory should be removed before restarting server/workers.

E.g.
```sh
rm -r output/binning_v2/tauTau_2018/workers
```
